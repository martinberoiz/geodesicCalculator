//
//  RKIntegrator.cpp
//  OutDirectionCalculator
//
//  Created by Developer Bot on 11/23/14.
//
//

#include "RKIntegrator.h"

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
(maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
(minarg1) : (minarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


RKIntegrator::RKIntegrator(int n_variables,
                           double* xvalues,
                           double** yvalues,
                           derivFunc derivatives) :
h_firstGuess(0.0001), h_min(0.), nok(0), nbad(0), eps(1E-4)
{
    setNumberOfVariables(n_variables);
    xp = xvalues;
    yp = yvalues;
    derivs = derivatives;
}

RKIntegrator::~RKIntegrator() {
    if (dydx != NULL) free(dydx);
    if (y != NULL) free(y);
    if (yscal != NULL) free(yscal);
    
    if (yerr != NULL) free(yerr);
    if (ytemp != NULL) free(ytemp);
    
    if (ak6 != NULL) free(ak6);
    if (ak5 != NULL) free(ak5);
    if (ak4 != NULL) free(ak4);
    if (ak3 != NULL) free(ak3);
    if (ak2 != NULL) free(ak2);

}

//Setters
void RKIntegrator::setNumberOfVariables(int nvariables) {
    nvar = nvariables;
    yscal = (double *)malloc(nvar*sizeof(*yscal));
    y = (double *)malloc(nvar*sizeof(*y));
    dydx = (double *)malloc(nvar*sizeof(*dydx));
    
    //To be used in decideQualityStepRK
    yerr = (double *)malloc(nvar*sizeof(*yerr));
    ytemp = (double *)malloc(nvar*sizeof(*ytemp));
    
    //To be used in integrateStep_CashKarpRK
    ak2 = (double *)malloc(nvar*sizeof(*ak2));
    ak3 = (double *)malloc(nvar*sizeof(*ak3));
    ak4 = (double *)malloc(nvar*sizeof(*ak4));
    ak5 = (double *)malloc(nvar*sizeof(*ak5));
    ak6 = (double *)malloc(nvar*sizeof(*ak6));

}
void RKIntegrator::setDerivativeFunction(derivFunc aDeriv) {derivs = aDeriv;}
void RKIntegrator::setAccuracy(double new_epsilon) {eps = new_epsilon;}
void RKIntegrator::setInitialConditions(double *new_ystart) {ystart = new_ystart;}


#define MAXSTP 1E6
#define TINY 1.0e-30
//My driver routine
int RKIntegrator::integrateInterval(double x1, double x2, double min_record_step, int max_number_points)
{
    //last saved value of x
    double xsav;
    //h is the current proposed step
    double h;
    //hdid is the step taken by the stepper
    double hdid;
    //hnext is the recommended next step by the stepper
    double hnext;
    
    int& kmax = max_number_points;
    double& dxsav = min_record_step;
    
    //Initialize x, the value of x on each iteration
    double x = x1;
	
    int kount = 0;
    h = SIGN(h_firstGuess, x2 - x1);
    for (int i = 0; i < nvar; i++) y[i] = ystart[i];
    if (kmax > 0) xsav = x - dxsav*2.0;
    for (int nstp = 0; nstp < MAXSTP; nstp++) {
        (*derivs)(x,y,dydx);
        for (int i = 0; i < nvar; i++)
            yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + TINY;
        if (kmax > 0 && kount < kmax - 1 && fabs(x - xsav) > fabs(dxsav)) {
            xp[kount] = x;
            for (int i = 0; i < nvar; i++) yp[i][kount] = y[i];
            xsav = x; kount++;
        }
        //If the stepsize h makes x to fall outside the interval (x1,x2),
        //adjust h so that x+h=x2
        if ((x + h - x2)*(x + h - x1) > 0.0) h = x2 - x;
        
        decideQualityStepRK(y, dydx, &x, h, eps, yscal, &hdid, &hnext);
        if (hdid == h) ++nok; else ++nbad;
        if ((x - x2)*(x2 - x1) >= 0.0) {
            for (int i = 0; i < nvar; i++) ystart[i] = y[i];
            if (kmax) {
                xp[kount] = x;
                for (int i = 0; i < nvar; i++) yp[i][kount] = y[i];
                kount++;
            }
            return kount;
        }
        if (fabs(hnext) <= h_min) printf("Step size too small in odeint\n");
        h = hnext;
    }
    printf("Too many steps in method integrateOnInterval\n");
    return kount;
    
}

void RKIntegrator::integrateUntilCondition(terminatingCond condition, double min_record_step, int max_number_points)
{
    //h is the current proposed step
    double h = h_firstGuess;
    //hdid is the step taken by the stepper
    double hdid;
    //hnext is the recommended next step by the stepper
    double hnext;
    
    int& kmax = max_number_points;
    double& dxsav = min_record_step;
    
    //Initialize x and y, the values on each iteration
    double x = 0;
    for (int i = 0; i < nvar; i++) y[i] = ystart[i];
    int kount = 0;
    
    //last saved value of x (set to a value that will save first x)
    double xsav = x - 2.0*dxsav;
    
    for (int nstp = 0; nstp < MAXSTP; nstp++) {
        (*derivs)(x,y,dydx);
        if (!(*condition)(nvar,x,y,dydx)) return;
        if (fabs(x - xsav) > fabs(dxsav)) {
            if (kmax > 0 && kount >= kmax) return;
            xp[kount] = x;
            for (int i = 0; i < nvar; i++) yp[i][kount] = y[i];
            xsav = x; kount++;
        }
        for (int i = 0; i < nvar; i++)
            yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + TINY;
        decideQualityStepRK(y, dydx, &x, h, eps, yscal, &hdid, &hnext);
        if (hdid == h) ++nok; else ++nbad;
        //if (fabs(hnext) <= h_min) printf("Step size too small in odeint\n");
        h = hnext;
    }
	printf("Too many steps in method integrateUntilCondition\n");
    
}
#undef MAXSTP
#undef TINY


#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
//That quality-control decision we encode in a stepper routine.
//The stepper routine calls the algorithm routine.
//It may reject the result, set a smaller stepsize,
// and call the algorithm routine again, until compatibility with a
// predetermined accuracy criterion has been achieved.
//The stepper’s fundamental task is to take the largest stepsize
// consistent with specified performance.
void RKIntegrator::decideQualityStepRK(double *y, double *dydx, double *x, double htry, double eps,
          double *yscal, double *hdid, double *hnext)
{
    
    int i;
    double errmax,h,htemp,xnew;
    
    h=htry;
    for (;;) {
        integrateStep_CashKarpRK(y,dydx,*x,h,ytemp,yerr);
        errmax=0.0;
        for (i=0; i < nvar; i++) errmax = FMAX(errmax, fabs(yerr[i]/yscal[i]));
        errmax /= eps;
        if (errmax > 1.0) {
            htemp = SAFETY*h*pow(errmax, PSHRNK);
            h = (h >= 0.0 ? FMAX(htemp, 0.1*h) : FMIN(htemp, 0.1*h));
            xnew = (*x) + h;
            //if (xnew == *x) printf("stepsize underflow in rkqs\n");
            continue;
        } else {
            if (errmax > ERRCON) *hnext = SAFETY*h*pow(errmax, PGROW);
            else *hnext = 5.0*h;
            *x += (*hdid = h);
            for (i = 0; i < nvar; i++) y[i] = ytemp[i];
            break;
        }
    }
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON


//The lowest or “nitty-gritty” level is the piece we call the algorithm routine.
//This implements the basic formulas of the method, starts with
//dependent variables y i at x, and calculates new values of the dependent variables at the value x + h.
//The algorithm routine also yields up some information about the quality of the solution after the step.
//The routine is dumb, however, and it is unable to make any adaptive decision
//about whether the solution is of acceptable quality or not.
void RKIntegrator::integrateStep_CashKarpRK(double *y, double *dydx, double x, double h, double *yout,
          double *yerr)
{
    static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
    double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
    
    for (int i = 0; i < nvar; i++)
        ytemp[i] = y[i] + b21*h*dydx[i];
    (*derivs)(x + a2*h, ytemp, ak2);
    for (int i = 0; i < nvar; i++)
        ytemp[i] = y[i] + h*(b31*dydx[i] + b32*ak2[i]);
    (*derivs)(x + a3*h, ytemp, ak3);
    for (int i = 0; i < nvar; i++)
        ytemp[i] = y[i] + h*(b41*dydx[i] + b42*ak2[i] + b43*ak3[i]);
    (*derivs)(x + a4*h, ytemp, ak4);
    for (int i = 0; i < nvar; i++)
        ytemp[i] = y[i] + h*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] + b54*ak4[i]);
    (*derivs)(x + a5*h, ytemp, ak5);
    for (int i = 0; i < nvar; i++)
        ytemp[i] = y[i] + h*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] + b64*ak4[i] + b65*ak5[i]);
    (*derivs)(x + a6*h, ytemp, ak6);
    for (int i = 0; i < nvar; i++)
        yout[i] = y[i] + h*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]);
    for (int i = 0;i < nvar; i++)
        yerr[i] = h*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] + dc5*ak5[i] + dc6*ak6[i]);
}


