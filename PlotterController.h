//
//  PlotterController.h
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 9/22/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import <Cocoa/Cocoa.h>
@class PlotterView;

@interface PlotterController : NSWindowController {

	IBOutlet PlotterView *theView;
	
}

-(void) setVerticesFromPointer:(double **)sphCoordArray andNumberOfVertices:(int)n;

@property(readonly) PlotterView *theView;

@end
