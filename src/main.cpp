#include <math.h>
#include <chrono>
#include <limits>
#include <vector>
#include <cassert>
#include <sstream>
#include <fstream>
#include <iostream>

// gurobi
#include "gurobi_c++.h"


using namespace std;
using matrix = vector<vector<double>>;

struct element
{
  int idx;
  double value;
};

using SparseVec = vector<element>;
using Basis = vector<SparseVec>;

string itos(int i) {stringstream s; s << i; return s.str(); }

bool is_integer(double x)
{
  if (abs(x-trunc(x)) < 1e-9)
  {
    return true;
  }
  else
  {
    return false;
  }
}

class Body
{
  GRBEnv env = GRBEnv();
  GRBModel* model;

  vector<GRBVar> xvars;
  vector<GRBVar> yvars;
  vector<char> vartypes;
  vector<string> varnames;

  void initialize(const string filename)
  {
    /* disable console output */
    env.set(GRB_IntParam_OutputFlag, 0);

    /* read problem */
    model = new GRBModel(env, filename);

    /* set tolerance */
    model->set(GRB_DoubleParam_OptimalityTol, 1e-9);

    /* get problem size */
    int nvars  = model->get(GRB_IntAttr_NumVars);
    int nconss  = model->get(GRB_IntAttr_NumConstrs);

    /* duplicate variables */
    GRBVar* origvars = model->getVars();
    for (int j=0; j<nvars; j++)
    {
      double lb = origvars[j].get(GRB_DoubleAttr_LB);
      double ub = origvars[j].get(GRB_DoubleAttr_UB);
      char type = origvars[j].get(GRB_CharAttr_VType);
      string name = origvars[j].get(GRB_StringAttr_VarName);

      /* register the original variable type */
      vartypes.push_back(type);

      /* register the original variable name */
      varnames.push_back(name);

      /* change type to continuous */
      if (type=='N') // semi-integer -> semi-continuos
      {
        origvars[j].set(GRB_CharAttr_VType, 'S');
      }
      else
      {
        origvars[j].set(GRB_CharAttr_VType, 'C');
      }

      /* change objective function */
      origvars[j].set(GRB_DoubleAttr_Obj, 0.0);

      /* create new variable and save */
      yvars.push_back( model->addVar(lb, ub, 0.0, 'C', "Y"+itos(j+1)) );
      xvars.push_back( origvars[j] );
    }

    /* duplicate constraints */
    GRBConstr* origconss = model->getConstrs();
    for (int i=0; i<nconss; i++)
    {
      GRBLinExpr expr = 0;
      for (int j=0; j<nvars; j++)
      {
        double coeff = model->getCoeff(origconss[i], origvars[j]);
        if (coeff != 0.0) expr += coeff * yvars[j];
      }
      double rhs = origconss[i].get(GRB_DoubleAttr_RHS);
      char sense = origconss[i].get(GRB_CharAttr_Sense);

      model->addConstr(expr, sense, rhs, "CY_"+itos(i+1));
    }

    /* delete arrays */
    delete[] origvars;
    delete[] origconss;
  }

public:
  int nvars;    // number of variables (= dimension of the lattice)
  int latrank;  // rank of the lattice
  Basis basis; // lattice basis

  Body(const string filename)
  {
    /* initialize body */
    initialize(filename);

    /* initialize basis */
    nvars = xvars.size();
    vector<double> zeros(nvars, 0.0);
    for( int j = 0; j < nvars; ++j )
    {
      /* add basis vector for each integer/binary variable */
      if (vartypes[j] == 'C' || vartypes[j] == 'S') continue;
      SparseVec newvec(1, {j, 1.0});
      basis.push_back(newvec);
    }
    latrank = basis.size();
  }

  ~Body()
  {
    xvars.clear();
    yvars.clear();
  }

  double distance(int k, const vector<double> &w)
  /* calculates F_{k}(w) for k \in {1,...,n}. */
  {
    vector<double> dummy;
    return distance(k, w, dummy);
  }

  double distance(int k, const vector<double> &w, vector<double> &alpha)
  /* calculates F_{k}(w) for k \in {1,...,n}. Also returns the dual
  variables alpha_1, ..., alpha_{k-1}. */
  {
    int n = xvars.size();
    assert(w.size() == n);

    /* change objective function */
    GRBLinExpr obj = 0;
    for( int j = 0; j < n; ++j )
      obj+=w[j]*xvars[j] - w[j]*yvars[j];
    model->setObjective(obj, GRB_MAXIMIZE);

    /* add new rows */
    vector<GRBConstr> added_conss;
    for( int i = 1; i < k; ++i )
    {
       GRBLinExpr expr;
       SparseVec bi = basis[i-1];
       for( element e: bi )
         expr += e.value * xvars[e.idx] - e.value * yvars[e.idx];
       string name = "N_" + itos(i);
       GRBConstr cons = model->addConstr(expr, GRB_EQUAL, 0.0, name);
       added_conss.push_back(cons);
    }

    /* solve */
    model->optimize();

    /* get status */
    int status = model->get(GRB_IntAttr_Status);
    assert(status == 2); // solution should be optimal

    /* get solution */
    double bestsol = model->get(GRB_DoubleAttr_ObjVal);

    /* get dual variables for the added constraints */
    alpha.clear();
    for(int i = 0; i < added_conss.size(); ++i)
    {
      alpha.push_back( -added_conss[i].get(GRB_DoubleAttr_Pi) );
    }

    /* release transformed problem and added constraints */
    for(int i = 0; i < added_conss.size(); ++i)
      model->remove(added_conss[i]);
    added_conss.clear();

    return max(0.0, bestsol);
  }

  double distance(int k, int p)
  /* calculates F_{k}(b^p) for k,p \in {1,...,n}. */
  {
    vector<double> dummy;
    return distance(k, p, dummy);
  }

  double distance(int k, int p, vector<double> &alpha)
  /* calculates F_{k}(b^p) for k,p \in {1,...,n}. Also returns the dual
  variables alpha_1, ..., alpha_{k-1}. */
  {
    int n = xvars.size();

    /* change objective function */
    GRBLinExpr obj;
    SparseVec bp = basis[p-1];
    for( element e: bp )
      obj += e.value * xvars[e.idx] - e.value * yvars[e.idx];
    model->setObjective(obj, GRB_MAXIMIZE);

    /* add new rows */
    vector<GRBConstr> added_conss;
    for( int i = 1; i < k; ++i )
    {
        GRBLinExpr expr;
        SparseVec bi = basis[i-1];
        for( element e: bi )
          expr += e.value * xvars[e.idx] - e.value * yvars[e.idx];
        string name = "N_" + itos(i);
        GRBConstr cons = model->addConstr(expr, GRB_EQUAL, 0.0, name);
        added_conss.push_back(cons);
    }

    /* solve */
    model->optimize();

    /* get status */
    int status = model->get(GRB_IntAttr_Status);
    assert(status == 2); // solution should be optimal

    /* get solution */
    double bestsol = model->get(GRB_DoubleAttr_ObjVal);

    /* get dual variables for the added constraints */
    alpha.clear();
    for(int i = 0; i < added_conss.size(); ++i)
    {
      alpha.push_back( -added_conss[i].get(GRB_DoubleAttr_Pi) );
    }

    /* release transformed problem and added constraints */
    for(int i = 0; i < added_conss.size(); ++i)
      model->remove(added_conss[i]);
    added_conss.clear();

    return max(0.0, bestsol);
  }

  void print_best(unsigned long k)
  /* prints best k branching directions into file */
  {
    ofstream outputfile;
    outputfile.open("hyperplanes.txt");

    k = std::min(k, basis.size());

    for (int j=0; j<k; j++)
    {
      for (element e: basis[j])
      {
        string name = varnames[e.idx];
        double value = e.value;
        if (value > 0.0)
          outputfile << "+" << value << " " << name << " ";
        if (value < 0.0)
          outputfile << value << " " << name << " ";
      }
      outputfile << "\n";
    }

    outputfile.close();
  }

};

void print_basis
(
  Body &P
)
{
  for( int i = 0; i < P.latrank; ++i )
  {
    vector<double> bi(P.nvars, 0.0);
    for (element e: P.basis[i])
      bi[e.idx] =  e.value;
    cout << "b" << i+1 << " = [";
    for( int j = 0; j < P.nvars; ++j )
      cout << bi[j] << " ";
    cout << "]" << endl;
  }
}

void reduce
(
  Body &P,
  double eps,
  double tol
)
{
  int n = P.nvars;
  int latrank = P.latrank;

  double mu;
  double f, fp;
  vector<double> alphas;

  double fpp = -1;
  vector<double> f_stored(latrank+1, -1);

  int i = 1;
  while (i < latrank)
  {
    cout << "*** i = " << i << " " << endl;

    /* get f = F_i(b^i) */
    f = (f_stored[i] >= 0.0) ? f_stored[i] : P.distance(i,i);

    /* get fpp = F_{i+1}(b^{i+1}) */
    fpp = (fpp >= 0.0) ? fpp : P.distance(i+1, i+1, alphas);

    cout << "   F" << i << "(b" << i << ") = " << f << endl;
    cout << "   F" << i+1 << "(b" << i+1 << ") = " << fpp << endl;

    /* get mu and fp */  
    double alpha = alphas.back();  
    if ( is_integer(alpha) )
    {
      mu = trunc(alpha); fp = fpp;
    }
    else
    {
      /* create vectors 
        vec1 = b^{i+1}+ ceil(alpha)b^i
        vec2 = b^{i+1}+ floor(alpha)b^i */
      vector<double> vec1(n, 0.0), vec2(n, 0.0);
      for (element e: P.basis[i])
      {
        vec1[e.idx] += e.value;
        vec2[e.idx] += e.value;
      }
      for (element e: P.basis[i-1])
      {
        vec1[e.idx] += ceil(alpha)*e.value;
        vec2[e.idx] += floor(alpha)*e.value;
      }

      vector<double> alphas1, alphas2;
      double z_ceil = P.distance(i, vec1, alphas1);
      double z_floor = P.distance(i, vec2, alphas2);
      if (z_ceil < z_floor)
      {
        mu = ceil(alpha); fp = z_ceil; alphas = alphas1;
      }
      else
      {
        mu = floor(alpha); fp = z_floor; alphas = alphas2;
      }

    }

    cout << "   alpha = " << alpha << endl;
    cout << "   mu = " << mu << endl;

    /* update b^{i+1} */
    if (mu != 0)
    {
      // b^{i+1} \gets  b^{i+1} + mu * b^{i}
      for( element e: P.basis[i-1] )
      {
        bool present = false;
        for( int j=0; j<P.basis[i].size(); j++ )
        {
          if (P.basis[i][j].idx == e.idx)
          {
            P.basis[i][j].value += mu * e.value;
            present = true;
          }
        }

        if (not present)
          P.basis[i].push_back( {e.idx, mu * e.value} );
      }
    }

    /* do basis check */
    if ( fp + tol < (1-eps)*f )
    {
      cout << "   condition not satisfied. i--" << endl;

      /* interchange b^i and b^{i+1} */
      SparseVec tempcopy = P.basis[i];
      P.basis[i] = P.basis[i-1];
      P.basis[i-1] = tempcopy;

      /* clear obsolete saved values */
      f_stored[i] = -1;
      f_stored[i+1] = -1;

      /* save values for next iter */
      fpp = fp;
      alphas.resize(i-1);
      if( i == 1 ) fpp = -1; // reusing not valid at i=1

      /* go back a step */
      i = max(1,i-1);

    }
    else
    {
      cout << "   condition satisfied. i++" << endl;
      /* store useful values */
      f_stored[i] = f;
      f_stored[i+1] = fpp;

      /* clear fpp */
      fpp = -1;

      i++;
    }

  }

}


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

  /* parameters */
  double eps = 0.25;
  double tol = 1e-6;

  /* initialize body */
  string filename(argv[1]);
  Body P(filename);

  /* track reduction time: start clock */
  auto start = chrono::high_resolution_clock::now();

  /* reduce basis */
  reduce(P, eps, tol);

  /* track reduction time: stop clock */
  auto end = chrono::high_resolution_clock::now();

  /* print basis */
  cout << endl;
  cout << "Reduction finished. Reduced basis is: " << endl;
  print_basis(P);

  /* check basis */
  cout << endl;
  for( int i = 1; i < P.latrank; ++i )
  {
    double Fii = P.distance(i, i);
    double Fiipp = P.distance(i, i+1);
    if ((1.0-eps)*Fii > Fiipp + tol)
    {
      cout << "Condition " << i << " " << Fiipp << " >= " << (1.0-eps)*Fii << " not satisfied " << endl;
      assert(0);
    }
    
  }

  /* print hyperplane */
  P.print_best(2);

  auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
  cout << "~ Total time: " << duration.count() << " milliseconds" << endl;

}
