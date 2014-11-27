//
//  RKIntegrator.h
//  OutDirectionCalculator
//
//  Created by Developer Bot on 11/23/14.
//
//

#ifndef __OutDirectionCalculator__RKIntegrator__
#define __OutDirectionCalculator__RKIntegrator__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef void (*derivFunc)(double, double *, double *);
typedef int (*terminatingCond)(int nvar, double x, double *y, double *y_prime);

class RKIntegrator {
public:
    RKIntegrator(int n_variables,
                 double* xvalues = NULL,
                 double** yvalues = NULL,
                 derivFunc derivatives = NULL);
    ~RKIntegrator();
    
    int integrateInterval(double x1, double x2, double min_record_step = 0., int max_number_points = 0);
    void integrateUntilCondition(terminatingCond condition, double min_record_step = 0., int max_number_points = 0);
    void setNumberOfVariables(int nvariables);
    void setAccuracy(double new_epsilon);
    void setDerivativeFunction(derivFunc aDeriv);
    void setInitialConditions(double* ystart);
    
private:
    derivFunc derivs;
    double eps;
    double h_firstGuess, h_min;
    int nok, nbad;
    double *ystart;
    double *xp, **yp; // This is where the values will be saved (should be vectors to avoid overflow)
    int nvar;
    
    //The scaling factor required in the computation of the error.
    double* yscal;
    //The current values of y for each iteration
    double* y;
    //The current values of y' for each iteration
    double* dydx;
    double *yerr, *ytemp;
    double *ak2,*ak3,*ak4,*ak5,*ak6;

    void decideQualityStepRK(double *y, double *dydx, double *x, double htry, double eps,
                                double *yscal, double *hdid, double *hnext);
    void integrateStep_CashKarpRK(double *y, double *dydx, double x, double h, double *yout,
              double *yerr);
};

#endif /* defined(__OutDirectionCalculator__RKIntegrator__) */
