#include<iostream>
#include<math.h>
#include<chrono>
#include<string>
#include<stdexcept>
#include"nwsqp.h"

#define M 3
#define N 5
#define dim 8

using namespace std;

void initx(double* x){
  x[0] = 3.0;
  x[1] = 2.0;
  x[2] = 2.0;
  x[3] = 1.0;
  x[4] = 1.0;
}

void setHV(double* x, double* lambda, double* s, double* hv){
  double prod = 1.0;
  for(int i = 0; i < N; i++){
    prod *= x[i];
  }
  double ep = exp(prod);

  double ae = prod*prod;

  double** cae = new double*[N];
  double** pae = new double*[N];
  for(int i = 0; i < N; i++){
    cae[i] = new double[N];
    pae[i] = new double[N];
  }
  double curae;
  double purae;
  for(int i = 0; i < N; i++){
    if (x[i] != 0.0){
      curae = ae/x[i];
      purae = prod/x[i];
    } else{
      curae = 0.0;
      purae = 1.0;
      for(int k = 0; k < N; k++){
        if (i != k){
	  purae *= x[k];
	}
      }
    }
    for(int j = i; j < N; j++){
      double fc = curae;
      double fp = purae;
      if (x[j] != 0.0){
        fc /= x[j];
	fp /= x[j];
      } else if(i == j){
	fc = 1.0;
	fp = 0.0;
	for(int k = 0; k < N; k++){
          if(k != i){
            fc *= x[k]*x[k];
	  }
	}
      } else{
	fc = 0.0;
        fp = 1.0;
	for(int k = 0; k < N; k++){
	  if (k != i && k != j){
            fp *= x[k];
	  }
	}
      }
      cae[i][j] = fc;
      cae[j][i] = fc;
      pae[i][j] = fp;
      pae[j][i] = fp;
    }
  }
  double x02 = x[0]*x[0];
  double x12 = x[1]*x[1];
  double hvc;
  double x023 = x02*x[0]+x12*x[1]+1.0;

  for(int i = 0; i < N; i++){
    hvc = 0.0;
    for(int j = 0; j < N; j++){
      if(i != j){
	hvc += s[j]*pae[i][j]*ep;
	if ((i == 0 || i == 1) && (j == 0 || j == 1)){
          hvc -= s[j]*9.0*x02*x12;
	}
	if ((i == 1 || i == 2) && (j == 1 || j == 2)){
	  hvc -= s[j]*lambda[1];
	}
	if ((i == 3 || i == 4) && (j == 3 || j == 4)){
          hvc += 5.0*s[j]*lambda[1];
	}
      } else{
	hvc -= s[j]*2.0*lambda[0];
	if (i == 0 || i == 1){
          hvc -= s[j]*9.0*x[i]*x[i]*x[i]*x[i];
	  hvc -= s[j]*6.0*x023*x[i];
	  hvc -= s[j]*6.0*lambda[2]*x[i];
	}
      }
      hvc+= s[j]*cae[i][j]*ep;
    }
    hv[i] = hvc;
  }
  delete[] cae;
  delete[] pae;
}

void setJac(double* x, double** jac){
  for(int i = 0; i < N; i++){
    jac[0][i] = 2.0*x[i];
  }

  jac[1][0] = 0.0;
  jac[1][1] = x[2];
  jac[1][2] = x[1];
  jac[1][3] = -5.0*x[4];
  jac[1][4] = -5.0*x[3];

  jac[2][0] = 3.0*x[0]*x[0];
  jac[2][1] = 3.0*x[1]*x[1];
  jac[2][2] = 0.0;
  jac[2][3] = 0.0;
  jac[2][4] = 0.0;
}

double fval(double* x){
  double prod = 1.0;
  for(int i = 0; i < N; i++){
    prod *= x[i];
  }
  double tt = x[0]*x[0]*x[0]+x[1]*x[1]*x[1]+1.0;
  return exp(prod) - 0.5*tt*tt;
}

void cval(double* x, double* c){
  double c0 = -10.0;
  for(int i = 0; i < N; i++){
    c0 += x[i]*x[i];
  }
  c[0] = c0;

  c[1] = x[1]*x[2]-5.0*x[3]*x[4];

  c[2] = x[0]*x[0]*x[0]+x[1]*x[1]*x[1]+1.0;
}

void fgrad(double* x, double* fgx){
  double prod = 1.0;
  for(int i = 0; i < N; i++){
    if (x[i] != 0.0){
      prod *= x[i];
    } else {
      prod = 0.0;
      break;
    }
  }
  double ex = prod*exp(prod);
  double tt = x[0]*x[0]*x[0]+x[1]*x[1]*x[1]+1.0;
  for(int i = 0; i < N; i++){
    double cur = 1.0;
    if (x[i] != 0.0){
      cur = ex/x[i];
    } else {
      for(int j = 0; j < N; j++){
	if(j != i){
	  cur *= x[j];
	}
      }
    }
    if(i == 0 || i == 1){
      cur -= 3.0*x[i]*x[i]*tt;
    }
    fgx[i] = cur;
  }
}
