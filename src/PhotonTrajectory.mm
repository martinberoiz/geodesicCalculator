//
//  PhotonTrajectory.m
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 8/29/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import "PhotonTrajectory.h"
#import "rpoly.h"
#import "derivsAndGlobals.h"
#include "RKIntegrator.h"

@interface PhotonTrajectory () {
	RKIntegrator* trajGenerator;
    double *xp, **yp;
    int kmax;
}
@end

@implementation PhotonTrajectory

@synthesize bhSpin, massToDistanceRatio;
@synthesize impactParam, carterConst;
@synthesize numberOfPoints;
@synthesize trajectoryPointsArray;

-(void) setBhMass: (double)newMass {
	bhMass = newMass;
	massToDistanceRatio = bhMass/pulsarDistance;
}

-(void) setPulsarDistance: (double)newDistance {
	pulsarDistance = newDistance;
	massToDistanceRatio = bhMass/pulsarDistance;
}

-(double) bhMass {
	return bhMass;	
}

-(double) pulsarDistance {
	return pulsarDistance;
}

-(NSNumber *) theta0 {
	return [NSNumber numberWithDouble:theta0];
}

-(void) setTheta0:(NSNumber *)value {
	theta0 = [value doubleValue];
}

-(NSNumber *) phi0 {
	return [NSNumber numberWithDouble:phi0];
}

-(void) setPhi0:(NSNumber *)value {
	phi0 = [value doubleValue];
}
-(NSNumber *) thetaIn {
	return [NSNumber numberWithDouble:thetaIn];
}

-(void) setThetaIn:(NSNumber *)value {
	thetaIn	= [value doubleValue];
}

-(NSNumber *) phiIn	{
	return [NSNumber numberWithDouble:phiIn];
}

-(void) setPhiIn:(NSNumber *)value {
	phiIn = [value doubleValue];
}

-(NSNumber *) thetaOut	{
	return [NSNumber numberWithDouble:thetaOut];
}

-(void) setThetaOut:(NSNumber *)value {
	thetaOut = [value doubleValue];
}

-(NSNumber *) phiOut {
	return [NSNumber numberWithDouble:phiOut];
}

-(void) setPhiOut:(NSNumber *)value {
	phiOut = [value doubleValue];
}


-(id)init {
	
	if ((self = [super init])) {
		phiIn = 0;
		thetaIn = M_PI/2.;
		bhMass = 1.0;
		[self setPulsarDistance:5];
		bhSpin = 0.0;
		theta0 = M_PI/2.;
		phi0 = 0.;
		thetaOut = 0.;
		phiOut = 0.;
		trajectoryPointsArray = NULL;
		
		/*This is to save intermediate values */
		int nvars = 3;
		kmax = 100;
		xp = (double *)malloc(kmax*sizeof(*xp));
		yp = (double **)malloc(nvars*sizeof(**yp));
		yp[0] = (double *)malloc(nvars*kmax*sizeof(**yp));
		for (int i = 1; i < nvars; i++) yp[i] = yp[i - 1] + kmax;
        
        trajGenerator = new RKIntegrator(nvars, xp, yp);

	}
		
	return self;
}


-(id)calculateTrajectoryParameters {
	
	double x = 1./[self massToDistanceRatio];
	double cosThSq = pow(cos(thetaIn), 2);
	double sinThSq = 1 - cosThSq;
	double sinPhiSq = pow(sin(phiIn), 2);
	double cosPhiSq = 1 - sinPhiSq;
	double th0 = theta0;
	a = [self bhSpin];
	double xxm2 = x*(x - 2);
	double a4 = pow(a, 4);
	
	double coeff2, coeff1, coeff0;
	
	
	coeff0 = -2*(-a4*a4*sinThSq*sinPhiSq - 2*sinThSq*sinPhiSq*pow(x - 2,2)*pow(x,6) - 
				 a4*a*a*sinThSq*sinPhiSq*x*(-2 + 5*x) + 
				 a*a*pow(x,3)*(8*cosPhiSq*(-2*cosThSq + x) + sinPhiSq*(x - 2)*(4 + 4*sinThSq*x - 7*sinThSq*x*x)) + 
				 a4*x*x*(8*cosPhiSq*cosThSq + sinPhiSq*(4 + 12*sinThSq*x - 9*sinThSq*x*x)) + 
				 a*a*sinPhiSq*(a*a + xxm2)*(-a4*sinThSq - 2*a*a*sinThSq*xxm2 + 
											x*x*(-4 + 4*sinThSq*x - sinThSq*x*x))*cos(2*th0));
	
	coeff1 = -8*a*x*(-(sinPhiSq*(a*a*(-2 + cosThSq) - 2*xxm2)*(a*a + xxm2)) + 
					 cosPhiSq*(a*a + 2*xxm2)*(a*a*cosThSq + x*(-2*cosThSq + x))*pow(1./sin(th0),2) + 
					 a*a*cos(2*th0)*(cosThSq*sinPhiSq*(a*a + xxm2) + 
									 cosPhiSq*(a*a*cosThSq + x*(-2*cosThSq + x))*pow(1./sin(th0),2)));
	
	coeff2 = -((a*a + 2*xxm2 + a*a*cos(2*th0))*
			   (4*a4*cosPhiSq*cosThSq + 4*a4*sinPhiSq - 3*a4*cosThSq*sinPhiSq - 24*a*a*cosPhiSq*cosThSq*x - 
				16*a*a*sinPhiSq*x + 6*a*a*cosThSq*sinPhiSq*x + 4*a*a*cosPhiSq*x*x + 32*cosPhiSq*cosThSq*x*x + 8*a*a*cosPhiSq*cosThSq*x*x + 
				16*sinPhiSq*x*x + 8*a*a*sinPhiSq*x*x - 3*a*a*cosThSq*sinPhiSq*x*x - 16*cosPhiSq*pow(x,3) - 
				16*cosPhiSq*cosThSq*pow(x,3) - 16*sinPhiSq*pow(x,3) + 8*cosPhiSq*pow(x,4) + 4*sinPhiSq*pow(x,4) + 
				4*(a4*(cosPhiSq*cosThSq + (-1 + cosThSq)*sinPhiSq) - sinPhiSq*xxm2*xxm2 + 
				   a*a*x*((-2 + cosThSq)*sinPhiSq*(x - 2) + cosPhiSq*(-2*cosThSq + x)))*cos(2*th0) - 
				a*a*cosThSq*sinPhiSq*(a*a + xxm2)*cos(4*th0))*pow(1./sin(th0),4))/4.;
	

	double b1,b2;
	quadratic_solver(coeff2, coeff1, coeff0, &b1, &b2);
	if (phiIn > 0) b = b2;
	else b = b1; 
	
	[self setImpactParam:b];
	
	q = (2*cosThSq*x*(-8*a*a4*b + 16*a*b*xxm2*x - 4*pow(a,3)*b*x*(-8 + 3*x) + a*a*a4*(4 + 3*x) + 
					  2*pow(x - 2,2)*x*x*(2*b*b + pow(x,3)) + a4*(b*b*(4 - 3*x) + 4*x*(-4 - x + 2*x*x)) + 
					  a*a*x*(-4*b*b*pow(x - 2,2) + x*(16 - 16*x*x + 7*pow(x,3))) + a*a*x*(a4 + 4*a*b*x +
																						  (x - 4)*pow(x,3) - a*a*(b*b - 2*xxm2))*cos(2*th0)) + 
		 2*x*x*pow(a*a + xxm2,2)*(-a*a + 2*b*b + a*a*cos(2*th0))*pow(1./tan(th0),2))
	/(4*pow(a*a + xxm2,2)*(a*a*cosThSq + x*(-2*cosThSq + x)));
	
	[self setCarterConst:q];
	
	return self;
}


#define TINY 1E-6
-(id)calculateTrajectoryPoints {
		
    int kount = 0;
	double coeff[5], zeror[4], zeroi[4];
	double vstart[3];
    //Schwarzchild radius
    double rs = 1./(2.*[self bhMass]);
    //Horizon radius
    double rh = (rs + sqrt(rs*rs-4.*pow([self bhSpin],2)))/2.;
	//Initial radial coord. of the pulsar
	double r0 = 1./massToDistanceRatio;
	
	//phi_0
	vstart[2] = phi0;
	
	
	if(thetaIn < M_PI && thetaIn > M_PI/2.) isGoingUp = NO;
	else isGoingUp = YES;
	
    double dxsav = 0.;
    double lambda1 = 0.;
    double lambda2 = 1E6;

	/*shouldUsePsi being positive or negative tells wether r dot is >0 or <0.
	 If r dot is negative there will be a turning point => use psi
	 This is the derivative of r along r_in at r0 (nabla r times n_in | r0)*/
	double shouldUsePsi = sin(theta0)*sin(thetaIn)*cos(phi0 - phiIn) + cos(theta0)*cos(thetaIn);	
	
	if (shouldUsePsi >= 0) {
		
		//r0
		vstart[0] = r0;
		
		if (fabs(q) < TINY && fabs(b) < fabs(a)) {
			//call with r and theta to avoid the problem with thMin
			vstart[1] = theta0;
            trajGenerator->setDerivativeFunction(&derivsUsingRAndTheta);
            trajGenerator->setInitialConditions(vstart);
            kount = trajGenerator->integrateInterval(lambda1, lambda2, dxsav, kmax);
			//odeint(vstart, 3, 0, 1E6, 1E-4, 1E-2, 0., &nok, &nbad, &derivsUsingRAndTheta, &rkqs);
		} else {
			
			//call with r and chi
			findThetaMin();
			vstart[1] = asin(cos(theta0)/cos(thMin));
        	trajGenerator->setDerivativeFunction(&derivsUsingRAndChi);
            trajGenerator->setInitialConditions(vstart);
            kount = trajGenerator->integrateInterval(lambda1, lambda2, dxsav, kmax);
			//odeint(vstart, 3, 0, 1E6, 1E-4, 1E-2, 0., &nok, &nbad, &derivsUsingRAndChi, &rkqs);
			vstart[1] = acos(cos(thMin)*sin(vstart[1]));
			for(int i = 0; i < kount; i++) yp[1][i] = acos(cos(thMin)*sin(yp[1][i]));
		}
		
	} else {
				
		coeff[0] = 1.;
		coeff[1] = 0.;
		coeff[2] = a*a - b*b - q;
		coeff[3] = 2*(pow(a - b, 2) + q);
		coeff[4] = -a*a*q;
		int nsols = rpoly(coeff, 4, zeror, zeroi);
		
		//x_min is the largest positive real root
		double x_min = 0;
		for (int i = 0; i < nsols; i++) {
			if (fabs(zeroi[i]) < TINY) {
				if (zeror[i] > x_min) x_min = zeror[i];
			}
		}
				
		e = (b*b + x_min*x_min)/(b*b - x_min*x_min);
		p = 2.*b*b*x_min/(b*b - x_min*x_min);
		
		vstart[0] = -acos((p/r0 - 1)/e);
		
		if (fabs(q) < TINY && fabs(b) > fabs(a)) {
			//call with psi and theta to avoid the problem with thMin
			vstart[1] = theta0;
            trajGenerator->setDerivativeFunction(&derivsUsingPsiAndTheta);
            trajGenerator->setInitialConditions(vstart);
            kount = trajGenerator->integrateInterval(lambda1, lambda2, dxsav, kmax);
			//odeint(vstart, 3, 0, 1E6, 1E-4, 1E-2, 0., &nok, &nbad, &derivsUsingPsiAndTheta, &rkqs);
		} else {
			//call with psi and chi
			findThetaMin();
			vstart[1] = asin(cos(theta0)/cos(thMin));
            trajGenerator->setDerivativeFunction(&derivsUsingPsiAndChi);
            trajGenerator->setInitialConditions(vstart);
            kount = trajGenerator->integrateInterval(lambda1, lambda2, dxsav, kmax);
			//odeint(vstart, 3, 0, 1E6, 1E-4, 1E-2, 0., &nok, &nbad, &derivsUsingPsiAndChi, &rkqs);
			vstart[1] = acos(cos(thMin)*sin(vstart[1]));
			for(int i = 0; i < kount; i++) yp[1][i] = acos(cos(thMin)*sin(yp[1][i]));
		}
		for(int i = 0; i < kount; i++) yp[0][i] = p/(1 + e*cos(yp[0][i]));
    }
    
    //prune anything going into the horizon or NaN's
    for(int i = 0; i < kount; i++) {
        if (yp[0][i] < rh || yp[0][i] != yp[0][i]) {
            kount = i;
            break;
        }
    }
    
	[self setThetaOut:[NSNumber numberWithDouble:vstart[1]]];
	[self setPhiOut:[NSNumber numberWithDouble:vstart[2]]];
	[self setTrajectoryPointsArray:yp];
	[self setNumberOfPoints:kount];
	
	return self;	
}
#undef TINY


-(void) dealloc {

    free(xp);
    free(yp[0]);
    free(yp);
    delete trajGenerator;
	[super dealloc];
	
}




@end
