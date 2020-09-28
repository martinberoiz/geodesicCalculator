/*
 *  derivsAndGlobals.c
 *  OutDirectionCalculator
 *
 *  Created by Martin Beroiz on 9/23/10.
 *  Copyright 2010 UTB. All rights reserved.
 *
 */

#include "derivsAndGlobals.h"
#include <math.h>
#include <stdio.h>

#define YES 1
#define NO 0

#define TINY 1E-6
void findThetaMin(void) {
	
	double d;
	if (fabs(a) < TINY) {
		thMin = acos(sqrt(q/(q + b*b) + b*b*q*a*a/pow(b*b + q, 3)));
	} else {
		d = q - a*a + b*b;
		if (q > TINY) {
			thMin = acos(sqrt(sqrt(d*d + 4*a*a*q) - d)/(M_SQRT2*a));
		} 
		else if (q < -TINY) {
			//There is a chance in this case to not have a solution.
			if (-sqrt(d*d + 4*a*a*q) - d < 0) printf("Error calculating thMin\n");
			thMin = acos(sqrt(-sqrt(d*d + 4*a*a*q) - d)/(M_SQRT2*a));
		} 
		else if (fabs(q) < TINY) {
			if (fabs(b) > a) {
				thMin = acos(sqrt(sqrt(d*d + 4*a*a*q) - d)/(M_SQRT2*a));
			} else {
				thMin = acos(sqrt(-sqrt(d*d + 4*a*a*q) - d)/(M_SQRT2*a));
			}
		}
	}
		
}
#undef TINY


void derivsUsingRAndTheta(double t, double *v, double *dv) {
	
	double x = v[0];
	double th = v[1];
	double radicand;
	
	radicand = -4*a*b*x - x*((b*b + q)*(x - 2) - pow(x,3)) - a*a*(q - x*(2 + x));
	
	if (radicand < 0) {
		if (radicand > -1E-6) radicand = -radicand;
	}
	dv[0] = sqrt(radicand)/(x*x + a*a*pow(cos(th),2));
	
	radicand = q + pow(cos(th),2)*(a*a - b*b/pow(sin(th),2));
	
	if (radicand < 0) {
		/*if (radicand < -1E-6) {
		 printf("Warning: Trying to square root a negative number: %g\n", radicand);
		 printf("Setting it to positive now.\n");
		 }*/
		if (radicand > -1E-6) radicand = -radicand;
	}
	
	dv[1] = sqrt(radicand)/(x*x + a*a*pow(cos(th),2));
	if (isGoingUp == YES) dv[1] = -dv[1];
	
	dv[2] = (a*(2*x - a*b) + b*(a*a + x*(x - 2))*pow(1./sin(th),2))/
	((a*a + x*(x - 2))*(x*x + a*a*pow(cos(th),2)));
	
	
}

void derivsUsingRAndChi(double t, double *v, double *dv) {
	
	double x = v[0];
	double th = acos(cos(thMin)*sin(v[1]));
	double radicand;
	
	radicand = -4*a*b*x - x*((b*b + q)*(x - 2) - pow(x,3)) - a*a*(q - x*(2 + x));
	
	if (radicand < 0) {
		if (radicand > -1E-6) radicand = -radicand;
	}
	dv[0] = sqrt(radicand)/(x*x + a*a*pow(cos(th),2));
	
	//This is the deriv of chi wrt lambda
	radicand = pow(b/sin(thMin), 2) - a*a + pow(a*cos(thMin)*sin(v[1]), 2);
	
	if (radicand < 0) {
		/*if (radicand < -1E-6) {
		 printf("Warning: Trying to square root a negative number: %g\n", radicand);
		 printf("Setting it to positive now.\n");
		 }*/
		if (radicand > -1E-6) radicand = -radicand;
	}
	
	dv[1] = sqrt(radicand)/(x*x + pow(a*cos(thMin)*sin(v[1]), 2));
	if (isGoingUp == NO) dv[1] = -dv[1];
	
	
	dv[2] = (a*(2*x - a*b) + b*(a*a + x*(x - 2))*pow(1./sin(th),2))/
	((a*a + x*(x - 2))*(x*x + a*a*pow(cos(th),2)));
	
}


void derivsUsingPsiAndChi(double t, double *v, double *dv) {
	
	double psi = v[0];
	double x = p/(1 + e*cos(psi));
	double th = acos(sin(v[1])*cos(thMin));
	
	dv[0] = sqrt(-e*(a*a - b*b - q)*(e + pow(1./cos(psi/2.),2)) - 
				 (2*(pow(a - b,2) + q)*(e*e*(3 + e*cos(psi)) + (e*(3 + e*e)*pow(1./cos(psi/2.),2))/2.))/p + 
				 (a*a*q*(2*e*(1 + e*e)*pow(1./cos(psi/2.),2) + 
						 e*e*(6 + e*(4*cos(psi) + e*(2 - pow(sin(psi),2))))))/(p*p))/(fabs(e)*(x*x + pow(a*cos(th), 2)));
	
	//This is the deriv of chi wrt lambda
	double radicand = pow(b/sin(thMin), 2) - a*a + pow(a*cos(thMin)*sin(v[1]), 2);
	
	if (radicand < 0) {
		/*if (radicand < -1E-6) {
		 printf("Warning: Trying to square root a negative number: %g\n", radicand);
		 printf("Setting it to positive now.\n");
		 } */
		if (radicand > -1E-6) radicand = -radicand;
	}
	dv[1] = sqrt(radicand)/(x*x + pow(a*cos(thMin)*sin(v[1]), 2));
	if (isGoingUp == NO) dv[1] = -dv[1];
	
	//deriv of phi wrt lambda
	dv[2] = (a*(2*x - a*b) + b*(a*a + x*(x - 2))*pow(1./sin(th),2))/((a*a + x*(x - 2))*(x*x + a*a*pow(cos(th),2)));
	
}


void derivsUsingPsiAndTheta(double t, double *v, double *dv) {
	
	double psi = v[0];
	double x = p/(1 + e*cos(psi));
	double th =v[1];
	
	dv[0] = sqrt(-e*(a*a - b*b - q)*(e + pow(1./cos(psi/2.),2)) - 
				 (2*(pow(a - b,2) + q)*(e*e*(3 + e*cos(psi)) + (e*(3 + e*e)*pow(1./cos(psi/2.),2))/2.))/p + 
				 (a*a*q*(2*e*(1 + e*e)*pow(1./cos(psi/2.),2) + 
						 e*e*(6 + e*(4*cos(psi) + e*(2 - pow(sin(psi),2))))))/(p*p))/(fabs(e)*(x*x + pow(a*cos(th), 2)));
	
	//This is the deriv of theta wrt lambda
	double radicand = q + pow(cos(th),2)*(a*a - b*b/pow(sin(th),2));
	
	if (radicand < 0) {
		/*if (radicand < -1E-6) {
		 printf("Warning: Trying to square root a negative number: %g\n", radicand);
		 printf("Setting it to positive now.\n");
		 }*/
		if (radicand > -1E-6) radicand = -radicand;
	}	
	dv[1] = sqrt(radicand)/(x*x + a*a*pow(cos(th),2));
	if (isGoingUp == YES) dv[1] = -dv[1];
	
	//This is deriv of phi wrt lambda
	dv[2] = (a*(2*x - a*b) + b*(a*a + x*(x - 2))*pow(1./sin(th),2))/
	((a*a + x*(x - 2))*(x*x + a*a*pow(cos(th),2)));
	
}


