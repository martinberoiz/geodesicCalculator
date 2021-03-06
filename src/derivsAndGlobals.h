/*
 *  derivsAndGlobals.h
 *  OutDirectionCalculator
 *
 *  Created by Martin Beroiz on 9/23/10.
 *  Copyright 2010 UTB. All rights reserved.
 *
 */


/********* GLOBAL VARS **********/
//a is the Kerr parameter J/M in units of M (a = J/M^2)
double a;

//b is L/E in units of M
double b;

//q is Q/E^2, Q is the cartan parameter in units of M
double q;

//for the theta -> chi transf
double thMin;

//for the r -> psi transf
double p, e;

//a flag that indicates I should integrate with increasing
//or decreasing theta
int isGoingUp;

/********* END OF GLOBAL VARS **********/

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

void derivsUsingRAndTheta(double t, double *v, double *dv);
void derivsUsingRAndChi(double t, double *v, double *dv);
void derivsUsingPsiAndTheta(double t, double *v, double *dv);
void derivsUsingPsiAndChi(double t, double *v, double *dv);
void findThetaMin(void);

#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif