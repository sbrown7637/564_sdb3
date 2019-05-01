#include<iostream>
#include<math.h>
#include<chrono>
#include<string>
#include<stdexcept>
#include"ecac.h"
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

void acceptStep(double* x, double* lambda, double f, double* fgx,
		double* c, double* v, double* w, double** jac,
		double &delta, double &nu, bool &acc){
  double eta = 1e-8;
  double beta = 1e-8;
  double cn2 = 0.0;
  for(int i = 0; i < M; i++){
    cn2 += c[i]*c[i];
  }
  double* s = new double[N];
  double* xs = new double[N];
  for (int i = 0; i < N; i++){
    s[i] = v[i]+w[i];
    xs[i] = x[i]+s[i];
  }
  double* jl = new double[N];
  cjac(jac, lambda, jl, true);
  double* jvc = new double[M];
  cjac(jac, v, jvc, false);
  double jvcn2 = 0.0;
  for(int i = 0; i < M; i++){
    jvc[i] += c[i];
    jvcn2 += jvc[i]*jvc[i];
  }
  setJac(xs, jac);
  double fs = fval(xs);
  double* fgxs = new double[N];
  double* cs = new double[M];
  fgrad(xs, fgxs);
  cval(xs, cs);
  double* lambdas = new double[M];
  augsys(fgxs, true, jac, lambdas, false);
  delete[] fgxs;
  double* hva = new double[N];
  double* hvb = new double[N];
  setHV(x, lambda, v, hva);
  setHV(x, lambda, w, hvb);
  double papr = 0.0;
  for(int i = 0; i < N; i++){
    papr -= (fgx[i]+jl[i]+hva[i])*w[i];
    papr -= (fgx[i]+jl[i])*v[i];
    papr -= 0.5*v[i]*hva[i];
    papr -= 0.5*w[i]*hvb[i];
    papr -= (lambdas[i]-lambda[i])*jvc[i];
  }
  delete[] hva;
  delete[] hvb;
  delete[] jvc;
  delete[] jl;
  if(papr < -0.5*nu*(cn2-jvcn2)){
    nu = -2.0*papr/(cn2-jvcn2)+beta;
  }
  double pred = papr+nu*(cn2-jvcn2);
  double ared = f-fs+nu*cn2;
  for(int i = 0; i < M; i++){
    ared += lambda[i]*c[i];
    ared -= lambdas[i]*cs[i];
    ared -= nu*cs[i]*cs[i];
  }
  delete[] lambdas;
  delete[] cs;
  double rat = ared/pred;
  if(rat >= eta){
    for(int i = 0; i < N; i++){
      x[i] = xs[i];
    }
    if(rat >= 0.9){
      delta = fmax(7.0*norm(s, N), delta);
    } else if(rat >= 0.3){
      delta = fmax(2.0*norm(s, N), delta);
    }
    acc = true;
  } else{
    delta = 0.5*fmax(norm(v, N), norm(w, N));
  }
  delete[] xs;
  delete[] s;
}

void tgs(double* x, double* gl, double* v, double* lambda, double** jac,
	 double delta, int cgitmax, double cgtol, double* s){
  double* hv = new double[N];
  setHV(x, lambda, v, hv);
  for(int i = 0; i < N; i++){
    s[i] = 0.0;
    hv[i] += gl[i];
  }
  double* z = new double[N];
  augsys(hv, true, jac, z, true);
  double* p = new double[N];
  double rz = 0.0;
  for(int i = 0; i < N; i++){
    p[i] = -1.0*z[i];
    rz += hv[i]*z[i];
  }
  double rz0 = rz;
  double* hp = new double[N];
  for(int iter = 0; iter < cgitmax; iter++){
    setHV(x, lambda, p, hp);
    double php = 0.0;
    for(int i = 0; i < N; i++){
      php += p[i]*hp[i];
    }
    if(php <= 0.0){
      double a = 0.0;
      double b = 0.0;
      double sn2 = 0.0;
      for(int i = 0; i < N; i++){
	a += p[i]*p[i];
	b += s[i]*p[i];
	sn2 += s[i]*s[i];
      }
      double theta = (-1.0*b+sqrt(b*b-a*(sn2-delta*delta)))/a;
      for(int i = 0; i < N; i++){
	s[i] += theta*p[i];
      }
      break;
    }
    double alpha = rz/php;
    double nrs = 0.0;
    for(int i = 0; i < N; i++){
      double tcr = s[i]+alpha*p[i];
      nrs += tcr*tcr;
    }
    nrs = sqrt(nrs);
    if(nrs > delta){
      double a = 0.0;
      double b = 0.0;
      double sn2 = 0.0;
      for(int i = 0; i < N; i++){
	a += p[i]*p[i];
	b += s[i]*p[i];
	sn2 += s[i]*s[i];
      }
      double theta = (-1.0*b+sqrt(b*b-a*(sn2-delta*delta)))/a;
      for(int i = 0; i < N; i++){
	s[i] += theta*p[i];
      }
      break;
    }
    for(int i = 0; i < N; i++){
      s[i] += alpha*p[i];
      hv[i] += alpha*hp[i];
    }
    augsys(hv, true, jac, z, true);
    double rzn = 0.0;
    for(int i = 0; i < N; i++){
      rzn += hv[i]*z[i];
    }
    if(rzn <= cgtol*cgtol*rz0){
      break;
    }
    double beta = rzn/rz;
    for(int i = 0; i < N; i++){
      p[i] *= beta;
      p[i] -= z[i];
    }
    rz = rzn;
  }
  delete[] hp;
  delete[] p;
  delete[] z;
  delete[] hv;
}

//if fgflag is true, uses fgx with 0s, else uses 0s with c
//if lgflag is true, will get projection, else will get lagrange multipliers
void augsys(double* fgxc, bool fgflag, double** jac,
	    double* mults, bool lgflag){
  double* luOut = new double[M+N];
  luSolve(jac, fgxc, fgflag, luOut);
  if(lgflag){
    for(int i = 0; i < N; i++){
      mults[i] = luOut[i];
    }
  } else{
    for(int i = N; i < dim; i++){
      mults[i-N] = luOut[i];
    }
  }
  delete[] luOut;
}

void luSolve(double** A, double* b, bool fgflag, double* x){
  //note - assumes A has rows, then columns
  double** U = new double*[dim];
  double** L = new double*[dim];
  for(int i = 0; i < dim; i++){
    U[i] = new double[dim];
    L[i] = new double[dim];
  }
  int* perm = new int[dim];
  //the ith entry is the row corresponding to column i for which P has a 1
  for(int col = 0; col < dim; col++){
    perm[col] = col; //initially, P is identity
    L[col][col] = 1.0; //so is L
    for(int row = 0; row < dim; row++){
      U[col][row] = getAugVal(A, col, row);
      if(row != col){
	L[row][col] = 0.0;
      }
      //initially, U is the augmented jacobian
    }
  }
  for(int k = 0; k < dim-1; k++){
    int curb = k;
    double curv = fabs(U[k][k]);
    for(int i = k; i < dim; i++){
      if(fabs(U[i][k]) > curv){
	curb = i;
	curv = fabs(U[i][k]);
      }
    }
    if(curb != k){
      for(int i = k; i < dim; i++){
	double temp = U[k][i];
	U[k][i] = U[curb][i];
	U[curb][i] = temp;
      }
      for(int i = 0; i < k; i++){
	double temp = L[k][i];
	L[k][i] = L[curb][i];
	L[curb][i] = temp;
      }
      int tp = perm[k];
      perm[k] = perm[curb];
      perm[curb] = tp;
    }
    for(int j = k+1; j < dim; j++){
      L[j][k] = U[j][k]/U[k][k];
      for(int m = k; m < dim; m++){
        U[j][m] -= L[j][k]*U[k][m];
      }
    }
  }

  //find permuted b
  double* bp = new double[dim];
  if(fgflag){
    for(int i = 0; i < dim; i++){
      if(perm[i] < N){
        bp[i] = b[perm[i]];
      } else{
        bp[i] = 0.0;
      }
    }
  } else{
    for(int i = 0; i < dim; i++){
      if(perm[i] >= N){
	bp[i] = b[perm[i]-N];
      } else{
	bp[i] = 0.0;
      }
    }
  }
  delete[] perm;
  //solve Ly = Pb with forward substitution
  double* y = new double[dim];
  for(int i = 0; i < dim; i++){
    double rem = bp[i];
    for(int j = 0; j < i; j++){
      rem -= y[j]*L[i][j];
    }
    y[i] = rem/L[i][i];
  }
  delete[] bp;
  for(int i = 0; i < dim; i++){
    delete[] L[i];
  }
  delete[] L;
  //solve Ux = y with backward substitution
  for(int i = dim-1; i >= 0; i--){
    double rem = y[i];
    for(int j = dim-1; j > i; j--){
      rem -= x[j]*U[i][j];
    }
    x[i] = rem/U[i][i];
  }
  delete[] U;
  delete[] y;
}

double getAugVal(double** A, int row, int col){
  if(row == col && col < N){
    return 1.0;
  } else if(row < N && col >= N){
    return A[col-N][row];
  } else if(row >= N && col < N){
    return A[row-N][col];
  }
  return 0;
}

void cjac(double** A, double* v, double* cv, bool adj){
  if(adj){
    for(int i = 0; i < N; i++){
      double tv = 0.0;
      for(int j = 0; j < M; j++){
	tv += v[j]*A[j][i];
      }
      cv[i] = tv;
    }
  } else{
    for(int i = 0; i < M; i++){
      double tv = 0.0;
      for(int j = 0; j < N; j++){
	tv += v[j]*A[i][j];
      }
      cv[i] = tv;
    }
  }
}

void setGl(double* fgx, double* cj, double* gl){
  for(int i = 0; i < N; i++){
    gl[i] = fgx[i]-cj[i];
  }
}

double norm(double* v, int s){
  double no = 0.0;
  for(int i = 0; i < s; i++){
    no += v[i]*v[i];
  }
  return sqrt(no);
}

void qns(double* c, double delta, double* v, double** A){
  cjac(A, c, v, true);
  double* vcp = new double[M];
  cjac(A, v, vcp, false); 
  double nv = norm(v, N);
  double tce = nv/norm(vcp, M);
  delete[] vcp;
  double tcc = tce*tce;
  double* scp = new double[N];
  for(int i = 0; i < N; i++){
    scp[i] = -1.0*tcc*v[i];
  }

  double nscp = tcc*nv;
  if(nscp >= delta){
    double coe = delta/nscp;
    for(int i = 0; i < N; i++){
      v[i] = coe*scp[i];
    }
    delete[] scp;
    return;
  }

  double* tmp = new double[M];
  augsys(c, false, A, tmp, false);
  double* sn = new double[N];
  cjac(A, tmp, sn, true);
  delete[] tmp;
  double nsn = norm(sn, N);
  if(nsn <= delta){
    for (int i = 0; i < N; i++){
      v[i] = sn[i];
    }
  } else{
    double a = 0.0;
    double b = 0.0;
    for (int i = 0; i < N; i++){
      v[i] = sn[i]-scp[i];
      a += v[i]*v[i];
      b += v[i]*scp[i];
    }
    double c = nscp*nscp-delta*delta;
    double tau = (-1.0*b+sqrt(b*b-a*c))/a;
    for (int i = 0; i < N; i++){
      v[i] *= tau;
      v[i] += scp[i];
    }
  }
  delete[] scp;
  delete[] sn;
}

