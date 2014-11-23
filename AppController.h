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
@class PlotterView;

@interface AppController : NSObject {
	PhotonTrajectory *aTraj;
	IBOutlet PlotterView *theView;
}

-(void)updateView;
+(void) initialiseValueTransformers;
-(void)setVerticesFromPointer:(double **)sphCoordArray andNumberOfVertices:(int) n;


@end
