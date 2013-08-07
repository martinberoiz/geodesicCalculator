//
//  PlotterView.h
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 9/22/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import <Cocoa/Cocoa.h>

typedef struct {
	GLfloat x;
	GLfloat y;
	GLfloat z;	
} cartesianVector;

@interface PlotterView : NSOpenGLView {

	int numberOfVertices;
	GLfloat *vertexArray;
	GLfloat eyeX, eyeY,eyeZ;
	GLfloat bhMass;
	GLfloat pulsarPositionX, pulsarPositionY, pulsarPositionZ;

}

@property(readwrite, assign) int numberOfVertices;
@property(nonatomic, readwrite) GLfloat* vertexArray;
@property(nonatomic ,readwrite, assign) GLfloat eyeX, eyeY, eyeZ;
@property(readwrite, assign) GLfloat bhMass;
@property(readwrite, assign) GLfloat pulsarPositionX, pulsarPositionY, pulsarPositionZ;

@end
