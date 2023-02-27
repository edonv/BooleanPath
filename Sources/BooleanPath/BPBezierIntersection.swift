//
//  BPBezierIntersection.swift
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

/// ``BPBezierIntersection`` stores where two bezier curves intersect.
///
/// Initially it just stores the curves and the parameter values where they intersect.
///
/// It can lazily compute the 2D point where they intersect,
/// the left and right parts of the curves relative to
/// the intersection point, and whether the intersection is tangent.
public class BPBezierIntersection {
    static let pointCloseThreshold = 1e-7
    static let parameterCloseThreshold = 1e-4
    
    fileprivate var _location: CGPoint?
    fileprivate var _curve1: BPBezierCurve
    fileprivate var _parameter1: Double
    fileprivate var _curve1LeftBezier: BPBezierCurve?
    fileprivate var _curve1RightBezier: BPBezierCurve?
    fileprivate var _curve2: BPBezierCurve
    fileprivate var _parameter2: Double
    fileprivate var _curve2LeftBezier: BPBezierCurve?
    fileprivate var _curve2RightBezier: BPBezierCurve?
    fileprivate var _tangent: Bool = false
    fileprivate var needToComputeCurve1 = true
    fileprivate var needToComputeCurve2 = true
    
    public var location: CGPoint {
        computeCurve1()
        return _location!
    }
    
    var curve1: BPBezierCurve {
        return _curve1
    }
    
    var parameter1: Double {
        return _parameter1
    }
    
    var curve2: BPBezierCurve {
        return _curve2
    }
    
    var parameter2: Double {
        return _parameter2
    }
    
    //+ (id) intersectionWithCurve1:(FBBezierCurve *)curve1 parameter1:(CGFloat)parameter1 curve2:(FBBezierCurve *)curve2 parameter2:(CGFloat)parameter2;
    //- (id) initWithCurve1:(FBBezierCurve *)curve1 parameter1:(CGFloat)parameter1 curve2:(FBBezierCurve *)curve2 parameter2:(CGFloat)parameter2;
    // let i = FBBezierIntersection(curve1: dvbc1, param1: p1, curve2: dvbc2, param2: p2)
    init(curve1: BPBezierCurve, param1: Double, curve2:BPBezierCurve, param2: Double) {
        _curve1 = curve1
        _parameter1 = param1
        _curve2 = curve2
        _parameter2 = param2
    }
    
    //- (BOOL) isTangent
    /// Returns `true` if intersection is tangent.
    public var isTangent: Bool {
        // If we're at the end of a curve, it's not tangent,
        // so skip all the calculations
        if isAtEndPointOfCurve {
            return false
        }
        computeCurve1()
        computeCurve2()
        
        // Compute the tangents at the intersection.
        let curve1LeftTangent  = PointMath.normalizePoint(PointMath.subtractPoint(_curve1LeftBezier!.controlPoint2, point2: _curve1LeftBezier!.endPoint2))
        let curve1RightTangent = PointMath.normalizePoint(PointMath.subtractPoint(_curve1RightBezier!.controlPoint1, point2: _curve1RightBezier!.endPoint1))
        let curve2LeftTangent  = PointMath.normalizePoint(PointMath.subtractPoint(_curve2LeftBezier!.controlPoint2, point2: _curve2LeftBezier!.endPoint2))
        let curve2RightTangent = PointMath.normalizePoint(PointMath.subtractPoint(_curve2RightBezier!.controlPoint1, point2: _curve2RightBezier!.endPoint1))
        
        // See if the tangents are the same. If so, then we're tangent at the intersection point
        return ProximityMath.arePointsClose(curve1LeftTangent, point2: curve2LeftTangent, threshold: BPBezierIntersection.pointCloseThreshold)
        || ProximityMath.arePointsClose(curve1LeftTangent, point2: curve2RightTangent, threshold: BPBezierIntersection.pointCloseThreshold)
        || ProximityMath.arePointsClose(curve1RightTangent, point2: curve2LeftTangent, threshold: BPBezierIntersection.pointCloseThreshold)
        || ProximityMath.arePointsClose(curve1RightTangent, point2: curve2RightTangent, threshold: BPBezierIntersection.pointCloseThreshold)
    }
    
    //- (FBBezierCurve *) curve1LeftBezier
    var curve1LeftBezier: BPBezierCurve {
        computeCurve1()
        return _curve1LeftBezier!
    }
    
    //- (FBBezierCurve *) curve1RightBezier
    var curve1RightBezier: BPBezierCurve {
        computeCurve1()
        return _curve1RightBezier!
    }
    
    //- (FBBezierCurve *) curve2LeftBezier
    var curve2LeftBezier: BPBezierCurve {
        computeCurve2()
        return _curve2LeftBezier!
    }
    
    //- (FBBezierCurve *) curve2RightBezier
    var curve2RightBezier: BPBezierCurve {
        computeCurve2()
        return _curve2RightBezier!
    }
    
    //- (BOOL) isAtStartOfCurve1
    var isAtStartOfCurve1: Bool {
        return ProximityMath.areValuesClose(_parameter1, value2: 0.0, threshold: BPBezierIntersection.parameterCloseThreshold) || _curve1.isPoint
    }
    
    //- (BOOL) isAtStopOfCurve1
    var isAtStopOfCurve1: Bool {
        return ProximityMath.areValuesClose(_parameter1, value2: 1.0, threshold: BPBezierIntersection.parameterCloseThreshold) || _curve1.isPoint
    }
    
    //- (BOOL) isAtEndPointOfCurve1
    var isAtEndPointOfCurve1: Bool {
        return self.isAtStartOfCurve1 || self.isAtStopOfCurve1
    }
    
    //- (BOOL) isAtStartOfCurve2
    var isAtStartOfCurve2: Bool {
        return ProximityMath.areValuesClose(_parameter2, value2: 0.0, threshold: BPBezierIntersection.parameterCloseThreshold) || _curve2.isPoint
    }
    
    //- (BOOL) isAtStopOfCurve2
    var isAtStopOfCurve2: Bool {
        return ProximityMath.areValuesClose(_parameter2, value2: 1.0, threshold: BPBezierIntersection.parameterCloseThreshold) || _curve2.isPoint
    }
    
    //- (BOOL) isAtEndPointOfCurve2
    var isAtEndPointOfCurve2: Bool {
        return self.isAtStartOfCurve2 || self.isAtStopOfCurve2
    }
    
    //- (BOOL) isAtEndPointOfCurve
    var isAtEndPointOfCurve: Bool {
        return self.isAtEndPointOfCurve1 || self.isAtEndPointOfCurve2
    }
    
    //- (void) computeCurve1
    fileprivate func computeCurve1() {
        if needToComputeCurve1 {
            let pap = _curve1.pointAtParameter(_parameter1)
            _location = pap.point
            _curve1LeftBezier = pap.leftBezierCurve
            _curve1RightBezier = pap.rightBezierCurve
            needToComputeCurve1 = false
        }
    }
    
    //- (void) computeCurve2
    fileprivate func computeCurve2() {
        if needToComputeCurve2 {
            let pap = _curve2.pointAtParameter(_parameter2)
            // not using the point from curve2
            _curve2LeftBezier = pap.leftBezierCurve
            _curve2RightBezier = pap.rightBezierCurve
            needToComputeCurve2 = false
        }
    }
}
