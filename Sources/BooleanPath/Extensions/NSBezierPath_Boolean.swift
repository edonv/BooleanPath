//
//  NSBezierPath_Boolean.swift
//  BooleanPath
//
//  Oligin is NSBezierPath+Boolean - Created by Andrew Finnell on 2011/05/31.
//  Copyright 2011 Fortunate Bear, LLC. All rights reserved.
//
//  Based on VectorBoolean - Created by Leslie Titze on 2015/05/19.
//  Copyright (c) 2015 Leslie Titze. All rights reserved.
//
//  Created by Takuto Nakamura on 2019/03/10.
//  Copyright © 2019 Takuto Nakamura. All rights reserved.
//

#if canImport(Cocoa)
import Cocoa

public extension NSBezierPath {
    // 15
    //- (NSBezierPath *) fb_union:(NSBezierPath *)path
    func union(_ path: NSBezierPath) -> NSBezierPath {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = thisGraph.union(with: otherGraph).bezierPath
        result.copyAttributesFrom(self)
        return result
    }
    
    // 24
    //- (NSBezierPath *) fb_intersect:(NSBezierPath *)path
    func intersection(_ path: NSBezierPath) -> NSBezierPath {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = thisGraph.intersect(with: otherGraph).bezierPath
        result.copyAttributesFrom(self)
        return result
    }
    
    // 33
    //- (NSBezierPath *) fb_difference:(NSBezierPath *)path
    func difference(_ path: NSBezierPath) -> NSBezierPath {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = thisGraph.difference(with: otherGraph).bezierPath
        result.copyAttributesFrom(self)
        return result
    }
    
    // 42
    //- (NSBezierPath *) fb_xor:(NSBezierPath *)path
    func xor(_ path: NSBezierPath) -> NSBezierPath {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = thisGraph.xor(with: otherGraph).bezierPath
        result.copyAttributesFrom(self)
        return result
    }
}
#endif
