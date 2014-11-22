//
//  AppController.m
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 8/29/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import "AppController.h"
#import "DegreesToRadiansTransformer.h"
#import "PhotonTrajectory.h"
#import "PlotterView.h"


@implementation AppController

+(void) initialize {
	[super initialize];
		
	[self initialiseValueTransformers]; 
}

+(void) initialiseValueTransformers {

	DegreesToRadiansTransformer *dToRTransformer = [[[DegreesToRadiansTransformer alloc] init] autorelease];
	[DegreesToRadiansTransformer setValueTransformer:dToRTransformer 
											 forName:@"DegreesToRadiansTransformer"];

}

-(id)init {
	[super init];
	aTraj = [[PhotonTrajectory alloc] init];
	[aTraj setPhiIn:[NSNumber numberWithDouble: 2.0]];
	[aTraj setThetaIn:[NSNumber numberWithDouble: 1.7]];
	[aTraj setBhMass:1.0];
	[aTraj setPulsarDistance:5.0];
	[aTraj setBhSpin:0.0];
		
	[aTraj addObserver:self forKeyPath:@"bhSpin" options:NSKeyValueObservingOptionNew context:NULL];
	[aTraj addObserver:self forKeyPath:@"bhMass" options:NSKeyValueObservingOptionNew context:NULL];
	[aTraj addObserver:self forKeyPath:@"pulsarDistance" options:NSKeyValueObservingOptionNew context:NULL];
	[aTraj addObserver:self forKeyPath:@"thetaIn" options:NSKeyValueObservingOptionNew context:NULL];
	[aTraj addObserver:self forKeyPath:@"phiIn" options:NSKeyValueObservingOptionNew context:NULL];
	[aTraj addObserver:self forKeyPath:@"theta0" options:NSKeyValueObservingOptionNew context:NULL];
	[aTraj addObserver:self forKeyPath:@"phi0" options:NSKeyValueObservingOptionNew context:NULL];
	
	//[self updateView];
	
	return self;
}

-(void)updateView {
	
	[aTraj calculateTrajectoryParameters];
	[aTraj calculateTrajectoryPoints];
	[self setVerticesFromPointer:[aTraj trajectoryPointsArray] andNumberOfVertices:[aTraj numberOfPoints]];
	[theView setNeedsDisplay:YES];
}


-(void)setVerticesFromPointer:(double **)sphCoordArray andNumberOfVertices:(int) n {
	
	int i;
	GLfloat *vertices;
	vertices = (GLfloat *)malloc(3*n*sizeof(*vertices));
	double M = [aTraj bhMass];
	
	for(i = 0; i < n; i++) {
		
		/*Beware: Here I load the	x coord of the traj in the -z of OpenGL,
		 y coord of the traj in the +x of OpenGL,
		 and z coord of the traj in the +y of OpenGL,
		 so that the scene looks like it's in the usual right handed with z upwards system*/
		vertices[i*3] = (GLfloat)M*sphCoordArray[0][i]*sin(sphCoordArray[1][i])*sin(sphCoordArray[2][i]);
		vertices[i*3 + 1] = (GLfloat)M*sphCoordArray[0][i]*cos(sphCoordArray[1][i]);
		vertices[i*3 + 2] = -(GLfloat)M*sphCoordArray[0][i]*sin(sphCoordArray[1][i])*cos(sphCoordArray[2][i]);
	}
	
	[theView setVertexArray:vertices];
	[theView setNumberOfVertices:n];
	
}

- (void)observeValueForKeyPath:(NSString *)keyPath 
					  ofObject:(id)object 
						change:(NSDictionary *)change 
					   context:(void *)context {

	if ([keyPath isEqual:@"bhSpin"] || 
		[keyPath isEqual:@"thetaIn"] || 
		[keyPath isEqual:@"phiIn"]) {
		[self updateView];
	}	
		
	if ([keyPath isEqual:@"bhMass"]) {
		[theView setBhMass:(GLfloat) [[change valueForKey:@"new"] doubleValue]];
		[self updateView];
	}
	
	double newPosition, theta0;
	if ([keyPath isEqual:@"pulsarDistance"]) {
		theta0 = [[aTraj theta0] doubleValue];
		newPosition = [[change valueForKey:@"new"] doubleValue];
		theView.pulsarPositionX = 0.f;
		theView.pulsarPositionY = (GLfloat)newPosition*cos(theta0);
		theView.pulsarPositionZ = -(GLfloat)newPosition*sin(theta0);
		[self updateView];
	}

	double newTheta0, x, y, z, r;
	if ([keyPath isEqual:@"theta0"]) {
		x = [theView pulsarPositionX];
		y = [theView pulsarPositionY];
		z = [theView pulsarPositionZ];
		r = sqrt(x*x + y*y + z*z);
		
		newTheta0 = [[change valueForKey:@"new"] doubleValue];
		theView.pulsarPositionX = 0.f;
		theView.pulsarPositionY = (GLfloat)r*cos(newTheta0);
		theView.pulsarPositionZ = -(GLfloat)r*sin(newTheta0);
		[self updateView];
	}
	


}


@end
