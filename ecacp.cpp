#include<iostream>
#include<math.h>
#include<chrono>
#include<string>
#include<stdexcept>
#include"ecacp.h"
#include<omp.h>

#define N 10
#define M 4
#define dim 14

using namespace std;

void initx(double* x){
#pragma omp parallel for
  for(int i = 0; i < N; i++){
    x[i] = 1.0;
  }
}

void setHV(double* x, double* lambda, double* s, double* hv){
#pragma omp parallel for
  for(int i = 0; i < N; i++){
    hv[i] = s[i]/(x[i]*x[i]);
  }
}

void setJac(double* x, double** jac){
  jac[0][0] = -1.2141;
  jac[0][1] = -0.7697;
  jac[0][2] = -1.0891;
  jac[0][3] = 1.5442;
  jac[0][4] = 0.5377;
  jac[0][5] = 0.3188;
  jac[0][6] = 3.5784;
  jac[0][7] = 0.7254;
  jac[0][8] = -0.1241;
  jac[0][9] - 0.6715;

  jac[1][0] = -1.1135;
  jac[1][1] = 0.3714;
  jac[1][2] = 0.0326;
  jac[1][3] = 0.0859;
  jac[1][4] = 1.8339;
  jac[1][5] = -1.3077;
  jac[1][6] = 2.7694;
  jac[1][7] = -0.0631;
  jac[1][8] = 1.4897;
  jac[1][9] = -1.2075;

  jac[2][0] = -0.0068;
  jac[2][1] = -0.2256;
  jac[2][2] = 0.5525;
  jac[2][3] = -1.4916;
  jac[2][4] = -2.2588;
  jac[2][5] = -0.4336;
  jac[2][6] = -1.3499;
  jac[2][7] = 0.7147;
  jac[2][8] = 1.4090;
  jac[2][9] = 0.7172;

  jac[3][0] = 1.5326;
  jac[3][1] = 1.1174;
  jac[3][2] = 1.1006;
  jac[3][3] = -0.7423;
  jac[3][4] = 0.8622;
  jac[3][5] = 0.3426;
  jac[3][6] = 3.0349;
  jac[3][7] = -0.2050;
  jac[3][8] = 1.4172;
  jac[3][9] = 1.6302;
}

double fval(double* x){
  double fv = 0.0;
  for(int i = 0; i < N; i++){
    fv -= log(x[i]);
  }
  return fv;
}

void cval(double* x, double* c){
  c[0] = -1.2141*x[0]-0.7697*x[1]-1.0891*x[2]+1.5442*x[3]+0.5377*x[4]
	 +0.3188*x[5]+3.5784*x[6]+0.7254*x[7]-0.1241*x[8]+0.6715*x[9]-1.0933;

  c[1] = -1.1135*x[0]+0.3714*x[1]+0.0326*x[2]+0.0859*x[3]+1.8339*x[4]
	 -1.3077*x[5]+2.7694*x[6]-0.0631*x[7]+1.4897*x[8]-1.2075*x[9]-1.1093;

  c[2] = -0.0068*x[0]-0.2256*x[1]+0.5525*x[2]-1.4916*x[3]-2.2588*x[4]
	 -0.4336*x[5]-1.3499*x[6]+0.7147*x[7]+1.4090*x[8]-0.7172*x[9]+0.8637;

  c[3] =  1.5326*x[0]+1.1174*x[1]+1.1006*x[2]-0.7423*x[3]+0.8622*x[4]
	 +0.3426*x[5]+3.0349*x[6]-0.2050*x[7]+1.4172*x[8]+1.6302*x[9]-0.0774;
}

void fgrad(double* x, double* fgx){
#pragma omp parallel for
  for(int i = 0; i < N; i++){
    fgx[i] = -1.0/x[i];
  }
}
