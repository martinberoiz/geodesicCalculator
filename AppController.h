//
//  AppController.h
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 8/29/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <math.h>

@class PhotonTrajectory;
@class PlotterController;

@interface AppController : NSObject {
	
	PhotonTrajectory *aTraj;
	PlotterController *plotterController;
	
}

-(IBAction)getDirection:(id)sender;
+(void) initialiseValueTransformers;

@end
