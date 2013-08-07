//
//  PlotterView.h
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 9/22/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import <Cocoa/Cocoa.h>


@interface PlotterView : NSOpenGLView {

	int numberOfVertices;
	GLfloat *vertexArray;
	GLfloat eyeX, eyeY,eyeZ;

}

@property(readwrite, assign) int numberOfVertices;
@property(readwrite, assign) GLfloat* vertexArray;
@property(readwrite, assign) GLfloat eyeX, eyeY, eyeZ;

@end
