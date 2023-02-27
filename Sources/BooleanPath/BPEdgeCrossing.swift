//
//  BPEdgeCrossing.swift
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

import CoreGraphics

class BPEdgeCrossing {
    fileprivate var _intersection: BPBezierIntersection
    
    var edge: BPBezierCurve?
    var counterpart: BPEdgeCrossing?
    var fromCrossingOverlap = false
    var entry = false
    var processed = false
    var selfCrossing = false
    var index: Int = 0
    
    init(intersection: BPBezierIntersection) {
        _intersection = intersection
    }
    
    var isProcessed: Bool {
        return processed
    }
    
    var isSelfCrossing: Bool {
        return selfCrossing
    }
    
    var isEntry: Bool {
        return entry
    }
    
    func removeFromEdge() {
        if let edge = edge {
            edge.removeCrossing(self)
        }
    }
    
    var order: Double {
        return parameter
    }
    
    var next: BPEdgeCrossing? {
        edge?.nextCrossing(self)
    }
    
    var previous: BPEdgeCrossing? {
        edge?.previousCrossing(self)
    }
    
    var nextNonself: BPEdgeCrossing? {
        var nextNon = next
        while let n = nextNon, n.isSelfCrossing {
            nextNon = n.next
        }
        return nextNon
    }
    
    var previousNonself: BPEdgeCrossing? {
        var prevNon = previous
        while let p = prevNon, p.isSelfCrossing {
            prevNon = p.previous
        }
        return prevNon
    }
    
    var parameter: Double {
        if edge == _intersection.curve1 {
            return _intersection.parameter1
        } else {
            return _intersection.parameter2
        }
    }
    
    var location: CGPoint {
        return _intersection.location
    }
    
    var curve: BPBezierCurve? {
        return edge
    }
    
    var leftCurve: BPBezierCurve? {
        guard !isAtStart else { return nil }
        
        if edge == _intersection.curve1 {
            return _intersection.curve1LeftBezier
        } else {
            return _intersection.curve2LeftBezier
        }
    }
    
    var rightCurve: BPBezierCurve? {
        guard !isAtEnd else { return nil }
        
        if edge == _intersection.curve1 {
            return _intersection.curve1RightBezier
        } else {
            return _intersection.curve2RightBezier
        }
    }
    
    var isAtStart: Bool {
        if edge == _intersection.curve1 {
            return _intersection.isAtStartOfCurve1
        } else {
            return _intersection.isAtStartOfCurve2
        }
    }
    
    var isAtEnd: Bool {
        if edge == _intersection.curve1 {
            return _intersection.isAtStopOfCurve1
        } else {
            return _intersection.isAtStopOfCurve2
        }
    }
}
