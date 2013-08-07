//
//  PlotterController.m
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 9/22/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import "PlotterController.h"
#import "PlotterView.h"

@implementation PlotterController

@synthesize theView;


-(id) init {

	if (![super initWithWindowNibName:@"Plotter"]) {
		return nil;
	}
	return self;
}


-(void)setVerticesFromPointer:(double **)sphCoordArray andNumberOfVertices:(int) n {
	
	int i;
	GLfloat *vertices;
	vertices = (GLfloat *)malloc(3*n*sizeof(*vertices));
	for(i = 0; i < n; i++) {
		
		/*Beware: Here I load the	x coord of the traj in the -z of OpenGL,
									y coord of the traj in the +x of OpenGL,
								and z coord of the traj in the +y of OpenGL,
		 so that the scene looks like it's in the usual right handed with z upwards system*/
		vertices[i*3] = (GLfloat)sphCoordArray[0][i]*sin(sphCoordArray[1][i])*sin(sphCoordArray[2][i]);
		vertices[i*3 + 1] = (GLfloat)sphCoordArray[0][i]*cos(sphCoordArray[1][i]);
		vertices[i*3 + 2] = -(GLfloat)sphCoordArray[0][i]*sin(sphCoordArray[1][i])*cos(sphCoordArray[2][i]);
	}
	
	[theView setVertexArray:vertices];
	[theView setNumberOfVertices:n];
	
}



@end
