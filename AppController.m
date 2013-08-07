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
#import "PlotterController.h"


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
	[aTraj setMassToDistanceRatio:0.2];
	[aTraj setBhSpin:0.0];
	return self;
}

-(IBAction)getDirection:(id)sender {
	
	[aTraj calculateTrajectoryParameters];
	[aTraj calculateTrajectoryPoints];
	
	if (!plotterController) {
		plotterController = [[PlotterController alloc] init];
	}
	
	[plotterController showWindow:self];
	
	[plotterController setVerticesFromPointer:[aTraj trajectoryPointsArray] andNumberOfVertices:[aTraj numberOfPoints]];
		
	[[plotterController theView] setNeedsDisplay:YES];
	
}


@end
