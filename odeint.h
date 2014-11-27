#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int kmax,kount;
double *xp,**yp,dxsav;

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

void rkqs(double *, double *, int, double *, double, double, double *,
		  double *, double *, void (*)(double, double *, double *));

void rkck(double *y, double *dydx, int n, double x, double h, double *yout,
		  double *yerr, void (*derivs)(double, double *, double *));

void odeint(double *ystart, int nvar, double x1, double x2, double eps, double h1,
			double hmin, int *nok, int *nbad,
			void (*derivs)(double, double *, double *),
			void (*rkqs)(double *, double *, int, double *, double, double, double *,
						 double *, double *, void (*)(double, double *, double *)));


#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
}
#endif
