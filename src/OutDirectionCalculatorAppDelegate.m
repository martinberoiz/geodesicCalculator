//
//  OutDirectionCalculatorAppDelegate.m
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 8/20/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import "OutDirectionCalculatorAppDelegate.h"
#import "AppController.h"


@implementation OutDirectionCalculatorAppDelegate

@synthesize window;

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
	// Insert code here to initialize your application
    [viewController updateView];
}

-(BOOL) applicationShouldTerminateAfterLastWindowClosed:(NSApplication *)theApplication {
	return YES;
}

@end
