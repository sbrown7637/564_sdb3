#ifndef ECACP_H
#define ECACP_H

void initx(double* x);
void setHV(double* x, double* lambda, double* s, double* hv);
void setJac(double* x, double** jac);
double fval(double* x);
void cval(double* x, double* c);
void fgrad(double* x, double* fgx);

#endif
