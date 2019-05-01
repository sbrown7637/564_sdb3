#include<iostream>
#include<math.h>
#include<chrono>
#include<string>
#include<stdexcept>
#include<omp.h>
#include"ecacp.h"
#include"SQPhpar.h"

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

void solveNLP(double* x);

int main(int argc, char* argv[]){
  omp_set_num_threads(4);
  auto start = std::chrono::system_clock::now();
  double* x = new double[N];
  initx(x);
  solveNLP(x);
  delete[] x;
  auto end = std::chrono::system_clock::now();
  cout << "time: " <<
    std::chrono::duration_cast<std::chrono::microseconds>(end-start).count() <<
    " us" << endl;
  return 1;
}

void solveNLP(double* x){
  double delta = 100.0;
  double nu = 1.0;
  double zeta = 0.8;

  double gtol = 1e-6;
  double ctol = 1e-6;
  double stol = 1e-8;

  int maxiter = 50;
  int curiter = 0;

  double** jac = new double*[M];
  for(int i = 0; i < M; i++){
    jac[i] = new double[N];
  }
  setJac(x, jac);
  double f = fval(x);
  double* c = new double[M];
  cval(x, c);
  double* fgx = new double[N];
  fgrad(x, fgx);
  double* lambda = new double[M];
  augsys(fgx, true, jac, lambda, false);
  double* gl = new double[N];
  double* cj = new double[N];
  cjac(jac, lambda, cj, true);
  setGl(fgx, cj, gl);
  double ngl = norm(gl, N);
  double nc = norm(c, M);
  double* v = new double[N];
  double* w = new double[N];
  while(curiter < maxiter &&
        ((ngl >= gtol) || (nc >= ctol)) &&
	delta >= stol){
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
    if(acc){
      curiter++;
      f = fval(x);
      cval(x, c);
      fgrad(x, fgx);
      augsys(fgx, true, jac, lambda, false);
      cjac(jac, lambda, cj, true);
      setGl(fgx, cj, gl);
      ngl = norm(gl, N);
      nc = norm(c, M);
    }
  }
  for(int i = 0; i < N; i++){
    cout << "x[" << i << "] = " << x[i] <<endl;
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
}
