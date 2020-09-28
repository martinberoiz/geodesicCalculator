//
//  UnaClase.m
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 9/14/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import "DegreesToRadiansTransformer.h"


@implementation DegreesToRadiansTransformer

+ (Class)transformedValueClass {
	return [NSNumber class]; 
}

+ (BOOL)allowsReverseTransformation { 
	return YES; 
}

- (id)transformedValue:(id)value {
	if (value != nil) {
		return [NSNumber numberWithFloat:[value floatValue]*180./M_PI ];		
	}
	NSBeep();
	return [NSNumber numberWithFloat:0.0];
}

- (id)reverseTransformedValue:(id)value {
	if (value != nil) {
		return [NSNumber numberWithFloat:[value floatValue]*M_PI/180.];		
	}
	NSBeep();
	return [NSNumber numberWithFloat:0.0];	
}

@end
