#include<iostream>
#include<math.h>
#include<chrono>
#include<string>
#include<stdexcept>
#include"nwsqp.h"
#include"SQPhelp.h"

#ifndef N
#error N must be defined!
#endif
#ifndef M
#error M must be defined!
#endif
#ifndef dim
#error dim must be defined!
#endif
#ifdef dim
#if dim != M+N
#error dim must be M+N!
#endif
#endif

using namespace std;

int main(int argc, char* argv[]){
  cout.precision(10);
  double* x = new double[N];
  initx(x);
  double delta = 100.0;
  double nu = 1.0;
  double zeta = 0.8;

  double gtol = 1e-6;
  double ctol = 1e-6;
  double stol = 1e-8;

  double** jac = new double*[M];
  for(int i = 0; i < M; i++){
    jac[i] = new double[N];
  }
  setJac(x, jac);

  if(jac[0][0] != 6.0 || jac[0][1] != 4.0 || jac[0][2] != 4.0 ||
     jac[0][3] != 2.0 || jac[0][4] != 2.0){
    cerr << "Error in defining first row of jacobian!" << endl;
    return 0;
  }
  if(jac[1][0] != 0.0 || jac[1][1] != 2.0 || jac[1][2] != 2.0 ||
     jac[1][3] != -5.0 || jac[1][4] != -5.0){
    cerr << "Error in defining second row of jacobian!" << endl;
    return 0;
  }
  if(jac[2][0] != 27.0 || jac[2][1] != 12.0 || jac[2][2] != 0.0 ||
     jac[2][3] != 0.0 || jac[2][4] != 0.0){
    cerr << "Error in defining third row of jacobian!" << endl;
    return 0;
  }

  double f = fval(x);

  if(f < 1.621e5 || f > 1.622e5){
    cerr << "Error in calculating function value" << endl;
    return 0;
  }

  double* c = new double[M];
  cval(x, c);

  if(c[0] != 9.0 || c[1] != -1.0 || c[2] != 36.0){
    cerr << "Error in defining constraint values" << endl;
    return 0;
  }

  double* fgx = new double[N];
  fgrad(x, fgx);

  if(fgx[0] < 6.500e5 || fgx[0] > 6.501e5 ||
     fgx[1] < 9.760e5 || fgx[1] > 9.761e5 ||
     fgx[2] < 9.765e5 || fgx[2] > 9.766e5 ||
     fgx[3] < 1.953e6 || fgx[3] > 1.954e6 ||
     fgx[4] < 1.953e6 || fgx[4] > 1.954e6){
    cerr << "Error in calculating gradient values" << endl;
    return 0;
  }

  double* lambda = new double[M];
  augsys(fgx, true, jac, lambda, false);

  if(lambda[0] <  4.065e5 || lambda[0] >  4.066e5 ||
     lambda[1] < -2.173e5 || lambda[1] > -2.172e5 ||
     lambda[2] < -5.830e4 || lambda[2] > -5.829e4){
    cerr << "Error in calculating lagrange multipliers" << endl;
    return 0;
  }

  double* gl = new double[N];
  double* cj = new double[N];
  cjac(jac, lambda, cj, true); 
  setGl(fgx, cj, gl);

  if(gl[0] < -2.152e5 || gl[0] > -2.151e5 ||
     gl[1] <  4.840e5 || gl[1] >  4.841e5 ||
     gl[2] < -2.152e5 || gl[2] > -2.151e5 ||
     gl[3] <  5.377e4 || gl[3] >  5.378e4 ||
     gl[4] <  5.377e4 || gl[4] >  5.378e4){
    cerr << "Error in calculating gl" << endl;
    return 0;
  }

  double ngl = norm(gl, N);

  if(ngl < 5.767e5 || ngl > 5.768e5){
    cerr << "error in calculating norm" << endl;
    return 0;
  }

  double nc = norm(c, M);
  double* v = new double[N];
  double* w = new double[N];

  bool acc = false;
  while(!acc && delta >= stol){
    qns(c, delta, v, jac); 
    int cgitmax = 200;
    double cgtol = 1e-2;
    tgs(x, gl, v, lambda, jac, delta, cgitmax, cgtol, w);
    acceptStep(x, lambda, f, fgx, c, v, w, jac, delta, nu, acc);
    if(!acc){
      setJac(x, jac);
    }
  }
  
  if(v[0] < -1.140    || v[0] > -1.139    ||
     v[1] < -4.370e-1 || v[1] > -4.369e-1 ||
     v[2] <  6.920e-2 || v[2] >  6.921e-2 ||
     v[3] < -1.736e-1 || v[3] > -1.735e-1 ||
     v[4] < -1.736e-1 || v[4] > -1.735e-1){ 
    cerr << "final quasinormal step incorrect" << endl;
    return 0;
  }

  if(w[0] < -1.541    || w[0] > -1.540    ||
     w[1] <  3.465    || w[1] >  3.466    ||
     w[2] < -1.541    || w[2] > -1.540    ||
     w[3] <  3.851e-1 || w[3] >  3.852e-1 ||
     w[4] <  3.851e-1 || w[4] >  3.852e-1){ 
    cerr << "final tangential step incorrect" << endl;
    return 0;
  }

  if(delta != 100){
    cerr << "final delta incorrect" << endl;
    return 0;
  }

  if(nu != 1){
    cerr << "final nu incorrect" << endl;
    return 0;
  }

  if(!acc){
    cerr << "final step not accepted" << endl;
    return 0;
  }

  delete[] w;
  delete[] v;
  delete[] gl;
  delete[] cj;
  delete[] lambda;
  delete[] fgx;
  delete[] c;
  for(int i = 0; i < M; i++){
    delete[] jac[i];
  }
  delete[] jac;
  delete[] x;
  return 1;
}
