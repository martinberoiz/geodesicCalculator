//
//  PlotterView.m
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 9/22/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import "PlotterView.h"
#import <GLUT/GLUT.h>

void unitSquare(void);

void unitSquare(void) {
	
	glVertex3f(5, 5, 0.0);
	glVertex3f(5, -5, 0.0);
	glVertex3f(-5, -5, 0.0);
	glVertex3f(-5, 5, 0.0);
	
}

@implementation PlotterView

@synthesize numberOfVertices, vertexArray;
@synthesize eyeX, eyeY, eyeZ;


-(void)prepare {
	
	NSOpenGLContext *glContext = [self openGLContext];
	[glContext makeCurrentContext];
	
	glClearColor(0, 0, 0, 0);
	glShadeModel(GL_FLAT);

}

-(id)initWithCoder:(NSCoder *)aDecoder {
	
	self = [super initWithCoder:aDecoder];
	[self setNumberOfVertices:0];
	[self setVertexArray:NULL];
	eyeX = 0;
	eyeY = 0;
	eyeZ = -30;

	[self prepare];
	return self;
}

-(void)reshape {
	NSRect baseRect = [self convertRectToBase:[self bounds]];
	glViewport(0, 0, baseRect.size.width, baseRect.size.height);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-1, 1, -1, 1, 1.5, 100);
	glMatrixMode(GL_MODELVIEW);
	
}

#define TINY 1E-6
-(void)drawRect:(NSRect)dirtyRect {
		
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();
	
	GLfloat upX, upY, upZ;
	if (fabs(eyeY) < TINY) {
		upX = 0;
		upY = 1;
		upZ = 0;
	} else if (eyeY > TINY) {
		upX = -eyeX;
		upZ = -eyeZ;
		upY = (eyeX*eyeX + eyeZ*eyeZ)/eyeY;
	} else if (eyeY < -TINY) {
		upX = eyeX;
		upZ = eyeZ;
		upY = -(eyeX*eyeX + eyeZ*eyeZ)/eyeY;		
	}
	
	//normalize up-vector:
	double upnorm = sqrt(upX*upX + upY*upY + upZ*upZ);
	upX /= upnorm;
	upY /= upnorm;
	upZ /= upnorm;

	gluLookAt(eyeX, eyeY, eyeZ , 0, 0, 0, upX, upY, upZ);
	
	//Draw the axes:
	glColor3f(1.0, 1.0, 1.0);
	glBegin(GL_LINES);
	glVertex3f(-100, 0, 0);
	glVertex3f(100, 0, 0);
	glVertex3f(0, -100, 0);
	glVertex3f(0, 100, 0);
	glVertex3f(0, 0, -100);
	glVertex3f(0, 0, 100);
	glEnd();
	
	
	glColor3f(1.0, 1.0, 0.0);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, vertexArray);
	glDrawArrays(GL_LINE_STRIP, 0, numberOfVertices);
	
	glColor3f(1.0, 0, 1);
	glPushMatrix();
	glRotatef(90, 1, 0, 0);
	glutWireSphere(2.0, 10, 5);
	glPopMatrix();
	
	glDisableClientState(GL_VERTEX_ARRAY);
	
	//Draw the pulsar:
	glColor3f(0.3, 0.3, 0.8);
	glPushMatrix();
	if (numberOfVertices != 0) glTranslatef(vertexArray[0], vertexArray[1], vertexArray[2]);
	else glTranslatef(-5, 0, 0);
	glRotatef(90, 1, 0, 0);
	glutWireSphere(0.5, 10, 5);
	glPopMatrix();

	
	
	glFlush();
	
}
#undef TINY

- (BOOL)acceptsFirstResponder {
	
    return YES;	
}

-(void)mouseDragged:(NSEvent *)theEvent {
	
	NSRect viewBounds = [self bounds];
	float dx = (float)[theEvent deltaX];
	float dy = (float)[theEvent deltaY];
	
	double radiusCameraPosition, thetaCameraPosition, phiCameraPosition;
	radiusCameraPosition = sqrt(eyeX*eyeX + eyeY*eyeY + eyeZ*eyeZ);
	thetaCameraPosition = acos(eyeY/radiusCameraPosition);
	phiCameraPosition = atan2(eyeZ, eyeX);
	
	phiCameraPosition += dx*2*M_PI/viewBounds.size.width;

	thetaCameraPosition += dy*M_PI/viewBounds.size.height;
	if (thetaCameraPosition < 0.01) thetaCameraPosition = 0.01;
	if (thetaCameraPosition > M_PI*0.99) thetaCameraPosition = M_PI*0.99;

	eyeY = (GLfloat)radiusCameraPosition*cos(thetaCameraPosition);
	eyeX = (GLfloat)radiusCameraPosition*sin(thetaCameraPosition)*cos(phiCameraPosition);
	eyeZ = (GLfloat) radiusCameraPosition*sin(thetaCameraPosition)*sin(phiCameraPosition);
	
	[self setNeedsDisplay:YES];
	
}

- (void)magnifyWithEvent:(NSEvent *)theEvent {
	
	double radiusCameraPosition, thetaCameraPosition, phiCameraPosition;
	radiusCameraPosition = sqrt(eyeX*eyeX + eyeY*eyeY + eyeZ*eyeZ);
	thetaCameraPosition = acos(eyeY/radiusCameraPosition);
	phiCameraPosition = atan2(eyeZ, eyeX);

	radiusCameraPosition *= (1 + [theEvent magnification]);
	
	eyeY = (GLfloat)radiusCameraPosition*cos(thetaCameraPosition);
	eyeX = (GLfloat)radiusCameraPosition*sin(thetaCameraPosition)*cos(phiCameraPosition);
	eyeZ = (GLfloat) radiusCameraPosition*sin(thetaCameraPosition)*sin(phiCameraPosition);
	
	[self setNeedsDisplay:YES];
	
}

@end
