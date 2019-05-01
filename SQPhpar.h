#ifndef SQPHPAR_H
#define SQPHPAR_H

#include<iostream>
#include<math.h>
#include<chrono>
#include<string>
#include<stdexcept>
#include"ecacp.h"

void acceptStep(double* x, double* lambda, double f, double* fgx,
		double* c, double* v, double* w, double** jac,
		double &delta, double &nu, bool &acc);
void tgs(double* x, double* gl, double* v, double* lambda, double** jac,
	 double delta, int cgitmax, double cgtol, double* s);
void augsys(double* fgxc, bool fgflag, double** jac,
	    double* mults, bool lgflag);
void luSolve(double** A, double* b, bool fgflag, double* x);
double getAugVal(double** A, int row, int col);
void cjac(double** A, double* v, double* cv, bool adj);
void setGl(double* fgx, double* cj, double* gl);
double norm(double* v, int s);
void qns(double* c, double delta, double* v, double** A);

#endif
