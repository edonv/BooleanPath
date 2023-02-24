//
//  BPEdgeOverlap.swift
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

class BPEdgeOverlap {
    private static let threshold = 1e-2
    
    var edge1: BPBezierCurve
    var edge2: BPBezierCurve
    fileprivate var _range: BPBezierIntersectRange
    
    var range: BPBezierIntersectRange {
        return _range
    }
    
    init(range: BPBezierIntersectRange, edge1: BPBezierCurve, edge2: BPBezierCurve) {
        _range = range
        self.edge1 = edge1
        self.edge2 = edge2
    }
    
    func fitsBefore(_ nextOverlap: BPEdgeOverlap) -> Bool {
        if ProximityMath.areValuesClose(range.parameterRange1.upperBound, value2: 1.0, threshold: BPEdgeOverlap.threshold) {
            let nextEdge = edge1.next
            return nextOverlap.edge1 == nextEdge
                && ProximityMath.areValuesClose(nextOverlap.range.parameterRange1.lowerBound, value2: 0.0, threshold: BPEdgeOverlap.threshold)
        }
        
        return nextOverlap.edge1 == edge1
            && ProximityMath.areValuesClose(nextOverlap.range.parameterRange1.lowerBound, value2: range.parameterRange1.upperBound, threshold: BPEdgeOverlap.threshold)
    }
    
    func fitsAfter(_ previousOverlap: BPEdgeOverlap) -> Bool {
        if ProximityMath.areValuesClose(range.parameterRange1.lowerBound, value2: 0.0, threshold: BPEdgeOverlap.threshold) {
            let previousEdge = edge1.previous
            
            return previousOverlap.edge1 == previousEdge
                && ProximityMath.areValuesClose(previousOverlap.range.parameterRange1.upperBound, value2: 1.0, threshold: BPEdgeOverlap.threshold)
        }
        
        return previousOverlap.edge1 == edge1
            && ProximityMath.areValuesClose(previousOverlap.range.parameterRange1.upperBound, value2: range.parameterRange1.lowerBound, threshold: BPEdgeOverlap.threshold)
    }
    
    func addMiddleCrossing() {
        let intersection = _range.middleIntersection
        
        let ourCrossing = BPEdgeCrossing(intersection: intersection)
        let theirCrossing = BPEdgeCrossing(intersection: intersection)
        
        ourCrossing.counterpart = theirCrossing
        theirCrossing.counterpart = ourCrossing
        
        ourCrossing.fromCrossingOverlap = true
        theirCrossing.fromCrossingOverlap = true
        
        edge1.addCrossing(ourCrossing)
        edge2.addCrossing(theirCrossing)
    }
    
    func doesContainParameter(_ parameter: Double,
                              onEdge edge:BPBezierCurve,
                              startExtends extendsBeforeStart: Bool,
                              endExtends extendsAfterEnd: Bool) -> Bool {
        if extendsBeforeStart && extendsAfterEnd {
            return true
        }
        
        var parameterRange: ClosedRange<Double>
        if edge == edge1 {
            parameterRange = _range.parameterRange1
        } else {
            parameterRange = _range.parameterRange2
        }
        
        let inLeftSide = extendsBeforeStart
            ? parameter >= 0.0
            : parameter > parameterRange.lowerBound
        let inRightSide = extendsAfterEnd
            ? parameter <= 1.0
            : parameter < parameterRange.upperBound
        
        return inLeftSide && inRightSide
    }
}
