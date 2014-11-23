//
//  OutDirectionCalculatorAppDelegate.h
//  OutDirectionCalculator
//
//  Created by Martin Beroiz on 8/20/10.
//  Copyright 2010 UTB. All rights reserved.
//

#import <Cocoa/Cocoa.h>
@class AppController;

@interface OutDirectionCalculatorAppDelegate : NSObject <NSApplicationDelegate> {
    NSWindow *window;
    IBOutlet AppController* viewController;
}

@property (assign) IBOutlet NSWindow *window;

-(BOOL) applicationShouldTerminateAfterLastWindowClosed:(NSApplication *)theApplication;
@end
