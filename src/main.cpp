#include <vector>
#include <sstream>
#include <iostream>

#include <scip/scip.h>
#include <scip/scipdefplugins.h>


using namespace std;
using matrix = vector<vector<double>>;

void CATCH_SCIP(SCIP_RETCODE retcode)
{
  if (retcode!=SCIP_OKAY)
  {
    throw "SCIP Error: recode ";
  }
}

class Polytope
{
  SCIP* scip;
  vector<SCIP_VAR*> xvars;
  vector<SCIP_VAR*> yvars;
  vector<SCIP_CONS*> conss;

  void initialize(const matrix &A, const vector<double> &b,
                  const vector<double> &lb, const vector<double> &ub)
  {
    SCIP_RETCODE retcode;
    int m = b.size();
    int n = lb.size();

    /* change objective sense to maximize */
    CATCH_SCIP( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

    /* create problem variables */
    ostringstream namebuf;
    for( int j = 0; j < n; ++j )
    {
      SCIP_VAR *xvar, *yvar;

      /* create the x SCIP_VAR object  and add to the scip problem */
      namebuf.str("");
      namebuf << "x#" << j;
      retcode = SCIPcreateVar(scip, &xvar, namebuf.str().c_str(),
                                   lb[j], SCIPinfinity(scip), 0.0,  // lb, ub, obj
                                   SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL);
      retcode = SCIPaddVar(scip, xvar);
      xvars.push_back(xvar);

      /* create the y SCIP_VAR object  and add to the scip problem */
      namebuf.str("");
      namebuf << "y#" << j;
      retcode = SCIPcreateVar(scip, &yvar, namebuf.str().c_str(),
                                   lb[j], SCIPinfinity(scip), 0.0,  // lb, ub, obj
                                   SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL);
      retcode = SCIPaddVar(scip, yvar);
      yvars.push_back(yvar);

    }

    /* create constraints Ax<=b */
    for( int i = 0; i < m; ++i )
    {
       SCIP_CONS* cons;
       namebuf.str("");
       namebuf<<"CX_"<<i;

       /* create SCIP_CONS object */
       retcode = SCIPcreateConsBasicLinear(scip, & cons, namebuf.str().c_str(),
                                                0, NULL, NULL, -SCIPinfinity(scip), b[i]); // nvars, vars, vals, lhs, rhs

       /* add the vars belonging to field in this row to the constraint */
       for( int j = 0; j < n; ++j )
          retcode = SCIPaddCoefLinear(scip, cons, xvars[j], A[i][j]);

       /* add the constraint to scip */
       retcode = SCIPaddCons(scip, cons);

       conss.push_back(cons);
    }

    /* create constraints Ay<=b */
    for( int i = 0; i < m; ++i )
    {
       SCIP_CONS * cons;
       namebuf.str("");
       namebuf<<"CY_"<<i;

       /* create SCIP_CONS object */
       retcode = SCIPcreateConsBasicLinear(scip, & cons, namebuf.str().c_str(),
                                                0, NULL, NULL, -SCIPinfinity(scip), b[i]); // nvars, vars, vals, lhs, rhs

       /* add the vars belonging to field in this row to the constraint */
       for( int j = 0; j < n; ++j )
          retcode = SCIPaddCoefLinear(scip, cons, yvars[j], A[i][j]);

       /* add the constraint to scip */
       retcode = SCIPaddCons(scip, cons);

       conss.push_back(cons);
    }

  }

public:
  matrix basis;

  Polytope(const matrix &A, const vector<double> &b,
           const vector<double> &lb, const vector<double> &ub)
  {
    /* initialize scip */
    CATCH_SCIP( SCIPcreate(&scip) );

    /* load default plugins */
    CATCH_SCIP( SCIPincludeDefaultPlugins(scip) );

    /* load solving settings */
    CATCH_SCIP( SCIPreadParams(scip, "settingsfile.set") );

    /* disable scip output to stdout */
    SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scip), TRUE);

    /* disable heuristics */
    CATCH_SCIP( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );

    /* create an empty problem */
    CATCH_SCIP( SCIPcreateProb(scip, "width", NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

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
    try
    {
      /* release variables */
       for(int i = 0; i < xvars.size(); ++i)
         CATCH_SCIP( SCIPreleaseVar(scip, &xvars[i]) );
       for(int i = 0; i < yvars.size(); ++i)
         CATCH_SCIP( SCIPreleaseVar(scip, &yvars[i]) );
       xvars.clear();
       yvars.clear();

       /* release constraints */
       for(int i = 0; i < conss.size(); ++i)
         CATCH_SCIP( SCIPreleaseCons(scip, &conss[i]) );
       conss.clear();

       /* free the scip problem */
       CATCH_SCIP( SCIPfree(&scip) );
    }
    catch ( const char* msg )
    {
       cerr << msg << endl;
    }
  }

  double distance(vector<double> w, int k, double* alpha, double* beta)
  {
    int n = xvars.size();

    /* change objective function */
    for( int j = 0; j < n; ++j )
    {
      CATCH_SCIP( SCIPchgVarObj(scip, xvars[j], w[j]) );
      CATCH_SCIP( SCIPchgVarObj(scip, yvars[j], -w[j]) );
    }

    /* add new rows */
    ostringstream namebuf;
    vector<SCIP_CONS*> added_conss;
    for( int i = 1; i < k; ++i )
    {
       SCIP_CONS* cons;
       namebuf.str("");
       namebuf<<"b_"<<i;

       /* create SCIP_CONS object */
       CATCH_SCIP( SCIPcreateConsBasicLinear(scip, & cons, namebuf.str().c_str(),
                                                0, NULL, NULL, 0.0, 0.0) ); // nvars, vars, vals, lhs, rhs

       /* add the vars belonging to field in this row to the constraint */
       for( int j = 0; j < n; ++j )
       {
         CATCH_SCIP( SCIPaddCoefLinear(scip, cons, xvars[j], basis[i-1][j]) );
         CATCH_SCIP( SCIPaddCoefLinear(scip, cons, yvars[j], -basis[i-1][j]) );
       }

       /* add the constraint to scip */
       CATCH_SCIP( SCIPaddCons(scip, cons) );

       added_conss.push_back(cons);
    }

    /* solve */
    CATCH_SCIP( SCIPsolve(scip) );

    /* get solution */
    SCIP_SOL *sol = SCIPgetBestSol(scip);
    double bestsol = SCIPgetSolOrigObj(scip, sol);

    /* get dual variables for last two constraints */
    SCIP_Bool success = FALSE;
    if ( k > 1 && alpha != nullptr)
    {
      SCIPconsGetDualsol(scip, added_conss[k-2], alpha, &success);
      assert(success);
    }
    if ( k > 2 && beta != nullptr)
    {
      SCIPconsGetDualsol(scip, added_conss[k-3], beta, &success);
      assert(success);
    }

    /* release transformed problem and added constraints */
    CATCH_SCIP( SCIPfreeTransform(scip) );
    for(int i = 0; i < added_conss.size(); ++i)
    {
      CATCH_SCIP( SCIPdelCons(scip, added_conss[i]) );
      CATCH_SCIP( SCIPreleaseCons(scip, &added_conss[i]) );
    }
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
    for( int j = 0; j < n; ++j )
    {
      CATCH_SCIP( SCIPchgVarObj(scip, xvars[j], basis[p-1][j]) );
      CATCH_SCIP( SCIPchgVarObj(scip, yvars[j], -basis[p-1][j]) );
    }

    /* add new rows */
    ostringstream namebuf;
    vector<SCIP_CONS*> added_conss;
    for( int i = 1; i < k; ++i )
    {
       SCIP_CONS* cons;
       namebuf.str("");
       namebuf<<"b_"<<i;

       /* create SCIP_CONS object */
       CATCH_SCIP( SCIPcreateConsLinear(scip, &cons, namebuf.str().c_str(),
                                        0, NULL, NULL, 0.0, 1.0, // nvars, vars, vals, lhs, rhs
 					                              TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );
       /* add the vars belonging to field in this row to the constraint */
       for( int j = 0; j < n; ++j )
       {
         CATCH_SCIP( SCIPaddCoefLinear(scip, cons, xvars[j], basis[i-1][j]) );
         CATCH_SCIP( SCIPaddCoefLinear(scip, cons, yvars[j], -basis[i-1][j]) );
       }

       /* add the constraint to scip */
       CATCH_SCIP( SCIPaddCons(scip, cons) );

       added_conss.push_back(cons);
    }

    /* solve */
    CATCH_SCIP( SCIPsolve(scip) );

    /* get solution */
    SCIP_SOL *sol = SCIPgetBestSol(scip);
    double bestsol = SCIPgetSolOrigObj(scip, sol);

    /* get dual variables for last two constraints */
    SCIP_Bool success = FALSE;
    double sss;
    cout << "nrows " << SCIPgetNLPRows(scip) <<endl;

    for(int i = 0; i < added_conss.size(); ++i)
    {

      SCIP_ROW* row = SCIPconsGetRow(scip, added_conss[i]);
      //assert(row!=nullptr);
      //sss = SCIProwGetDualsol(row);
      SCIPconsGetDualsol(scip, added_conss[i], &sss, &success);
      cout << sss <<endl;
    }
    if ( k > 1 && alpha != nullptr)
    {
      SCIPconsGetDualsol(scip, added_conss[k-2], alpha, &success);
      assert(success);
    }
    if ( k > 2 && beta != nullptr)
    {
      SCIPconsGetDualsol(scip, added_conss[k-3], beta, &success);
      assert(success);
    }

    /* release transformed problem and added constraints */
    CATCH_SCIP( SCIPfreeTransform(scip) );
    for(int i = 0; i < added_conss.size(); ++i)
    {
      CATCH_SCIP( SCIPdelCons(scip, added_conss[i]) );
      CATCH_SCIP( SCIPreleaseCons(scip, &added_conss[i]) );
    }
    added_conss.clear();

    return bestsol;
  }

};



/** main function for queens example */
int
main(
     int args,
     char ** argv
     )
{
  SCIP_RETCODE retcode;

  /* get problem data */
  int n = 2;
  int m = 3;
  double eps = 1/4;
  matrix A(m, vector<double>(n));
  vector<double> b(m);
  vector<double> lb(n);
  vector<double> ub(n);

  A = { {-1, -7}, {2, 7}, {-5, 4} };
  b = { -7, 14, 4 };
  lb = {0,0};
  ub = {1000, 1000};

  /* initialize polytope */
  Polytope P(A, b, lb, ub);


  double h, alpha, beta;
  vector<double> F(n+1, -1);
  double* prev_alpha = nullptr;

  int i = 1;
  while (i < n)
  {
    cout << i << endl;
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
    cout << "F" << i << "(b" << i << ")= " << F[i] << endl;
    cout << "F" << i+1 << "(b" << i+1 << ")= " << h << endl;
    cout << "alpha= " << alpha << endl;

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
      double h1 = P.distance(vec1, i+1, &alpha, &beta1);
      double h2 = P.distance(vec2, i+1, &alpha, &beta2);
      if (h1<h2)
      {
        mu = ceil(alpha); fpp = h1; beta = beta1;
      }
      else
      {
        mu = floor(alpha); fpp = h2; beta = beta2;
      }

    }

    /* update b^{i+1} */
    for( int j = 0; j < n; ++j )
      P.basis[i][j] = P.basis[i][j] + mu*P.basis[i-1][j];

    /* do basis check */
    if ( fpp < (1-eps)*F[i] )
    {
      cout << "back" << endl;
      /* interchange b^i and b^{i+1} */
      vector<double> tempcopy = P.basis[i];
      for( int j = 0; j < n; ++j )
      {
        P.basis[i][j] = P.basis[i-1][j];
        P.basis[i-1][j] = tempcopy[j];
      }

      /* update F[i] */
      F[i] = fpp;

      /* save correct dual variable */
      prev_alpha = &beta;

      /* go back a step */
      i = max(1,i-1);

    }
    else
    {
      cout << "forward" << endl;
      F[i+1] = h;
      i++;
    }

  }

}
