//
//  PhotonTrajectory.h
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 8/29/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface PhotonTrajectory : NSObject {
	
	double bhSpin;
	double bhMass;
	double pulsarDistance;
	double massToDistanceRatio;
	double thetaIn;
	double phiIn;
	double carterConst; //Carter parameter in units of M^2
	double impactParam; //Impact parameter in units of M
	double thetaOut;
	double phiOut;
	double theta0;
	double phi0;
	
	int numberOfPoints;
	double **trajectoryPointsArray;
    
}

-(id)calculateTrajectoryPoints;
-(id)calculateTrajectoryParameters;
@property(readwrite, assign) double bhSpin, massToDistanceRatio;
@property(readwrite, assign) double impactParam, carterConst;
@property(readwrite, assign) int numberOfPoints;
@property(readwrite, assign) double **trajectoryPointsArray;

-(double) bhMass;
-(double) pulsarDistance;
-(void) setBhMass:(double)newMass;
-(void) setPulsarDistance:(double)newDistance;
-(NSNumber *) theta0;
-(void) setTheta0:(NSNumber *)value;
-(NSNumber *) phi0;
-(void) setPhi0:(NSNumber *)value;
-(NSNumber *) thetaIn;
-(void) setThetaIn:(NSNumber *)value;
-(NSNumber *) phiIn;
-(void) setPhiIn:(NSNumber *)value;
-(NSNumber *) thetaOut;
-(void) setThetaOut:(NSNumber *)value;
-(NSNumber *) phiOut;
-(void) setPhiOut:(NSNumber *)value;


@end
