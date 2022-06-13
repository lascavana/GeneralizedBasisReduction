#include <math.h>
#include <limits>
#include <vector>
#include <cassert>
#include <sstream>
#include <fstream>
#include <iostream>


#include "gurobi_c++.h"


using namespace std;
using matrix = vector<vector<double>>;

string itos(int i) {stringstream s; s << i; return s.str(); }

void read_constraints
(
  const char *filename,
  matrix &A,
  vector<double> &b,
  vector<double> &lb,
  vector<double> &ub
)
{
  int n, m; double coeff; string str;

  /* open file for reading */
  ifstream input_file(filename, ifstream::in);

  /* read data */
  if (input_file.is_open())
  {
    input_file >> n >> m;

    /* read A matrix */
    for (int i = 0; i < m; i++)
  	{
        vector<double> row(n);
  	    for (int j = 0; j < n; j++)
          input_file >> row[j];
        A.push_back( row );
    }

    /* read rhs vector b */
    for (int i = 0; i < m; i++)
    {
      input_file >> coeff;
      b.push_back( coeff );
    }

    /* read variable lower bounds */
    for (int j = 0; j < n; j++)
    {
      input_file >> str;
      if (str == "-inf")
        lb.push_back( -numeric_limits<double>::max() );
      else
        lb.push_back( stod(str) );
    }

    /* read variable upper bounds */
    for (int j = 0; j < n; j++)
    {
      input_file >> str;
      if (str == "inf")
        ub.push_back( numeric_limits<double>::max() );
      else
        ub.push_back( stod(str) );
    }

    input_file.close();
  }
  else cerr << "Unable to open file" << endl;
}

class Polytope
{
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);

  vector<GRBVar> xvars;
  vector<GRBVar> yvars;

  void initialize(const matrix &A, const vector<double> &b,
                  const vector<double> &lb, const vector<double> &ub)
  {
    int m = b.size();
    int n = lb.size();
    assert(lb.size()==ub.size());

    /* disable console output */
    model.set(GRB_IntParam_OutputFlag, 0);

    /* create problem variables */
    for( int j = 0; j < n; ++j )
    {
      GRBVar xvar, yvar;
      string xname = "X" + to_string(j);
      string yname = "Y" + to_string(j);
      if (ub[j] < numeric_limits<double>::max())
      {
        xvar = model.addVar(lb[j], ub[j], 0.0, GRB_CONTINUOUS, xname);
        yvar = model.addVar(lb[j], ub[j], 0.0, GRB_CONTINUOUS, yname);
      }
      else
      {
        xvar = model.addVar(lb[j], GRB_INFINITY, 0.0, GRB_CONTINUOUS, xname);
        yvar = model.addVar(lb[j], GRB_INFINITY, 0.0, GRB_CONTINUOUS, yname);
      }
      xvars.push_back( xvar );
      yvars.push_back( yvar );
    }

    /* create constraints Ax<=b */
    for( int i = 0; i < m; ++i )
    {
      GRBLinExpr expr = 0;
      for( int j = 0; j < n; ++j )
        expr += A[i][j]*xvars[j];
      string name = "CX_" + itos(i);
      model.addConstr(expr, GRB_LESS_EQUAL, b[i], name);
    }

    /* create constraints Ay<=b */
    for( int i = 0; i < m; ++i )
    {
      GRBLinExpr expr = 0;
      for( int j = 0; j < n; ++j )
        expr += A[i][j]*yvars[j];
      string name = "CY_" + itos(i);
      model.addConstr(expr, GRB_LESS_EQUAL, b[i], name);
    }

  }

public:
  matrix basis;

  Polytope(const matrix &A, const vector<double> &b,
           const vector<double> &lb, const vector<double> &ub)
  {
    /* add basic variables and constraints */
    initialize(A, b, lb, ub);

    /* initialize basis to standard Zn basis */
    int n = lb.size();
    vector<double> zeros(n);
    for( int j = 0; j < n; ++j )
    {
      basis.push_back(zeros);
      basis[j][j] = 1;
    }

  }

  ~Polytope()
  {
    xvars.clear();
    yvars.clear();
  }

  double distance(vector<double> w, int k, double* alpha, double* beta)
  {
    int n = xvars.size();
    assert(w.size() == n);

    /* change objective function */
    GRBLinExpr obj = 0;
    for( int j = 0; j < n; ++j )
      obj+=w[j]*xvars[j] - w[j]*yvars[j];
    model.setObjective(obj, GRB_MAXIMIZE);

    /* add new rows */
    vector<GRBConstr> added_conss;
    for( int i = 1; i < k; ++i )
    {
       GRBLinExpr expr;
       for( int j = 0; j < n; ++j )
         expr += basis[i-1][j]*xvars[j] - basis[i-1][j]*yvars[j];
       string name = "N_" + itos(i);
       GRBConstr cons = model.addConstr(expr, GRB_EQUAL, 0.0, name);
       added_conss.push_back(cons);
    }

    /* solve */
    model.optimize();

    /* get solution */
    double bestsol = model.get(GRB_DoubleAttr_ObjVal);

    /* get dual variables for last two constraints */
    for(int i = 0; i < added_conss.size(); ++i)
    if ( k > 1 && alpha != nullptr)
      *alpha = -added_conss[k-2].get(GRB_DoubleAttr_Pi);
    if ( k > 2 && beta != nullptr)
      *beta = -added_conss[k-3].get(GRB_DoubleAttr_Pi);

    /* release transformed problem and added constraints */
    for(int i = 0; i < added_conss.size(); ++i)
      model.remove(added_conss[i]);
    added_conss.clear();

    return bestsol;
  }

  double distance(int k, int p, double* alpha, double* beta)
  /* calculates F_{k}(b^p) for k,p \in {1,...,n}. If k>=2, returns
  the dual vairable associated with b^{p-1} (alpha), and furthermore,
  if k>=3, returns the dual vairable associated with b^{p-2} (beta). */
  {
    int n = xvars.size();

    /* change objective function */
    GRBLinExpr obj;
    for( int j = 0; j < n; ++j )
      obj += basis[p-1][j]*xvars[j] - basis[p-1][j]*yvars[j];
    model.setObjective(obj, GRB_MAXIMIZE);

    /* add new rows */
    vector<GRBConstr> added_conss;
    for( int i = 1; i < k; ++i )
    {
       GRBLinExpr expr;
       for( int j = 0; j < n; ++j )
         expr += basis[i-1][j]*xvars[j] -basis[i-1][j]*yvars[j];
       string name = "N_" + itos(i);
       GRBConstr cons = model.addConstr(expr, GRB_EQUAL, 0.0, name);
       added_conss.push_back(cons);
    }

    /* solve */
    model.optimize();

    /* get solution */
    double bestsol = model.get(GRB_DoubleAttr_ObjVal);

    /* get dual variables for last two constraints */
    if ( k > 1 && alpha != nullptr)
      *alpha = -added_conss[k-2].get(GRB_DoubleAttr_Pi);
    if ( k > 2 && beta != nullptr)
      *beta = -added_conss[k-3].get(GRB_DoubleAttr_Pi);

    /* release transformed problem and added constraints */
    for(int i = 0; i < added_conss.size(); ++i)
      model.remove(added_conss[i]);
    added_conss.clear();

    return bestsol;
  }

};



/** main function for queens example */
int
main(
     int argc,
     char ** argv
     )
{
  if (argc < 2)
  {
    cout << "No file provided. Usage: reduce <inputfile>";
    return 0;
  }

  /* get problem data */
  matrix A;
  vector<double> b;
  vector<double> lb;
  vector<double> ub;
  read_constraints(argv[1], A, b, lb, ub);
  int m = A.size();
  int n = lb.size();

  /* parameters */
  double eps = 1/4;

  /* initialize polytope */
  Polytope P(A, b, lb, ub);

  double h, alpha, beta;
  vector<double> F(n+1, -1);
  double* prev_alpha = nullptr;

  /* reduction loop */
  int i = 1;
  while (i < n)
  {
    cout << "*** i = " << i << endl;
    /* get F[i] */
    if ( F[i] < 0 )
      F[i] = P.distance(i, i, nullptr, nullptr);

    /* get F[i+1] */
    if ( F[i+1] < 0 || prev_alpha == nullptr)
      h = P.distance(i+1, i+1, &alpha, &beta);
    else
    {
      h = F[i+1]; alpha = *prev_alpha;
      prev_alpha = nullptr;
    }
    cout << "   F" << i << "(b" << i << ") = " << F[i] << endl;
    cout << "   F" << i+1 << "(b" << i+1 << ") = " << h << endl;
    cout << "   alpha = " << alpha << endl;

    /* get mu */
    double mu, fpp;
    if ( trunc(alpha) == alpha )
    {
      mu = alpha; fpp = h;
    }
    else
    {
      double beta1, beta2;
      vector<double> vec1(n), vec2(n);
      for( int j = 0; j < n; ++j )
      {
        vec1[j] = P.basis[i][j] + ceil(alpha)*P.basis[i-1][j];
        vec2[j] = P.basis[i][j] + floor(alpha)*P.basis[i-1][j];
      }
      double h1 = P.distance(vec1, i, nullptr, &beta1);
      double h2 = P.distance(vec2, i, nullptr, &beta2);
      // cout << "F" << i << "(b_{i+1}+ceil(alpha)*b_{i}) = " << h1 <<endl;
      // cout << "F" << i << "(b_{i+1}+floor(alpha)*b_{i}) = " << h2 <<endl;
      if (h1<h2)
      {
        mu = ceil(alpha); fpp = h1; beta = beta1;
      }
      else
      {
        mu = floor(alpha); fpp = h2; beta = beta2;
      }

    }

    cout << "   mu = " << mu << endl;

    /* update b^{i+1} */
    for( int j = 0; j < n; ++j )
      P.basis[i][j] = P.basis[i][j] + mu*P.basis[i-1][j];

    /* do basis check */
    if ( fpp < (1-eps)*F[i] )
    {
      cout << "   condition not satisfied. i--" << endl;
      /* interchange b^i and b^{i+1} */
      vector<double> tempcopy = P.basis[i];
      for( int j = 0; j < n; ++j )
      {
        P.basis[i][j] = P.basis[i-1][j];
        P.basis[i-1][j] = tempcopy[j];
      }

      /* update F[i] */
      F[i] = fpp;

      /* forget F[i+1] if saved */
      F[i+1] = -1;
      assert(F[i+2]==-1);

      /* save correct dual variable */
      prev_alpha = &beta;

      /* go back a step */
      i = max(1,i-1);

    }
    else
    {
      cout << "   condition satisfied. i++" << endl;
      F[i+1] = h;
      i++;
    }

  }

  /* print basis */
  cout << endl;
  cout << "Reduction finished. Reduced basis is: " << endl;
  for( int i = 0; i < n; ++i )
  {
    cout << "b" << i+1 << " = [";
    for( int j = 0; j < n; ++j )
      cout << P.basis[i][j] << " ";
    cout << "]" << endl;
  }

  /* check basis */
  cout << endl;
  for( int i = 1; i < n; ++i )
  {
    double Fii = P.distance(i, i, nullptr, nullptr);
    double Fiipp = P.distance(i, i+1, nullptr, nullptr);
    cout << "Condition ";
    cout << Fiipp << " >= " << (1.0-eps)*Fii;
    if ((1.0-eps)*Fii <= Fiipp)
      cout << " satisfied " << endl;
    else
      cout << " not satisfied " << endl;
  }

}
