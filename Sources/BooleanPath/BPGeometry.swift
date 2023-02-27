//
//  BPGeometry.swift
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

import SwiftUI
import CoreGraphics

enum ThresholdConstants {
    static var pointCloseness: Double { 1e-10 }
    static var tangentCloseness: Double { 1e-12 }
    static var boundsCloseness: Double { 1e-9 }
}

// MARK: - Point Helpers

enum PointMath {
    static func distanceBetweenPoints(_ point1: CGPoint, point2: CGPoint) -> CGFloat {
        let xDelta = point2.x - point1.x
        let yDelta = point2.y - point1.y
        
        return sqrt(xDelta * xDelta + yDelta * yDelta);
    }
    
    static func distancePointToLine(_ point: CGPoint, lineStartPoint: CGPoint, lineEndPoint: CGPoint) -> CGFloat {
        let lineLength = PointMath.distanceBetweenPoints(lineStartPoint, point2: lineEndPoint)
        if lineLength == 0.0 {
            return 0.0
        }
        
        let xDelta = lineEndPoint.x - lineStartPoint.x
        let yDelta = lineEndPoint.y - lineStartPoint.y
        
        let num = point.x - lineStartPoint.x * xDelta + point.y - lineStartPoint.y * yDelta
        
        let u = num / (lineLength * lineLength)
        
        let intersectionPoint = CGPoint(x: lineStartPoint.x + CGFloat(u * xDelta),
                                        y: lineStartPoint.y + CGFloat(u * yDelta))
        
        return PointMath.distanceBetweenPoints(point, point2: intersectionPoint)
    }
    
    static func addPoint(_ point1: CGPoint, point2: CGPoint) -> CGPoint {
        return CGPoint(x: point1.x + point2.x,
                       y: point1.y + point2.y)
    }
    
    static func unitScalePoint(_ point: CGPoint, scale: CGFloat) -> CGPoint {
        var result = point
        let length = pointLength(point)
        
        if length != 0.0 {
            result.x = CGFloat(result.x * (scale / length))
            result.y = CGFloat(result.y * (scale / length))
        }
        return result
    }
    
    static func scalePoint(_ point: CGPoint, scale: CGFloat) -> CGPoint {
        return CGPoint(x: point.x * scale,
                       y: point.y * scale)
    }
    
    static func dotMultiplyPoint(_ point1: CGPoint, point2: CGPoint) -> CGFloat {
        let dotX = point1.x * point2.x
        let dotY = point1.y * point2.y
        return dotX + dotY
    }
    
    static func subtractPoint(_ point1: CGPoint, point2: CGPoint) -> CGPoint {
        return CGPoint(x: point1.x - point2.x,
                       y: point1.y - point2.y)
    }
    
    static func pointLength(_ point: CGPoint) -> CGFloat {
        let xSq = point.x * point.x
        let ySq = point.y * point.y
        return sqrt(xSq + ySq)
    }
    
    static func pointSquaredLength(_ point: CGPoint) -> CGFloat {
        let xSq = point.x * point.x
        let ySq = point.y * point.y
        return xSq + ySq
    }
    
    static func normalizePoint(_ point: CGPoint) -> CGPoint {
        var result = point
        let length = pointLength(point)
        if length != 0.0 {
            result.x = result.x / length
            result.y = result.y / length
        }
        return result
    }
    
    static func negatePoint(_ point: CGPoint) -> CGPoint {
        return CGPoint(x: -point.x,
                       y: -point.y)
    }
    
    static func roundPoint(_ point: CGPoint) -> CGPoint {
        return CGPoint(x: round(point.x),
                       y: round(point.y))
    }
    
    static func lineNormal(_ lineStart: CGPoint, lineEnd: CGPoint) -> CGPoint {
        return normalizePoint(CGPoint(x: lineStart.y - lineEnd.y, // -(lineEnd.y - lineStart.y)
                                      y: lineEnd.x - lineStart.x))
    }
    
    static func lineMidpoint(_ lineStart: CGPoint, lineEnd: CGPoint) -> CGPoint {
        let distance = distanceBetweenPoints(lineStart, point2: lineEnd)
        let tangent = normalizePoint(subtractPoint(lineEnd, point2: lineStart))
        return addPoint(lineStart, point2: unitScalePoint(tangent, scale: distance / 2.0))
    }
    
    static func rectGetTopLeft(_ rect: CGRect) -> CGPoint {
        return CGPoint(x: rect.minX,
                       y: rect.minY)
    }
    
    static func rectGetTopRight(_ rect: CGRect) -> CGPoint {
        return CGPoint(x: rect.maxX,
                       y: rect.minY)
    }
    
    static func rectGetBottomLeft(_ rect: CGRect) -> CGPoint {
        return CGPoint(x: rect.minX,
                       y: rect.maxY)
    }
    
    static func rectGetBottomRight(_ rect: CGRect) -> CGPoint {
        return CGPoint(x: rect.maxX,
                       y: rect.maxY)
    }
    
    static func expandBoundsByPoint(_ topLeft: inout CGPoint, bottomRight: inout CGPoint, point: CGPoint) {
        if point.x < topLeft.x { topLeft.x = point.x }
        if point.x > bottomRight.x { bottomRight.x = point.x }
        if point.y < topLeft.y { topLeft.y = point.y }
        if point.y > bottomRight.y { bottomRight.y = point.y }
    }
    
    static func unionRect(_ rect1: CGRect, rect2: CGRect) -> CGRect {
        var topLeft = rectGetTopLeft(rect1)
        var bottomRight = rectGetBottomRight(rect1)
        expandBoundsByPoint(&topLeft, bottomRight: &bottomRight, point: rectGetTopLeft(rect2))
        expandBoundsByPoint(&topLeft, bottomRight: &bottomRight, point: rectGetTopRight(rect2))
        expandBoundsByPoint(&topLeft, bottomRight: &bottomRight, point: rectGetBottomRight(rect2))
        expandBoundsByPoint(&topLeft, bottomRight: &bottomRight, point: rectGetBottomLeft(rect2))
        
        return CGRect(x: topLeft.x,
                      y: topLeft.y,
                      width: bottomRight.x - topLeft.x,
                      height: bottomRight.y - topLeft.y)
    }
}

// MARK: - Proximity Helpers

enum ProximityMath {
    static func arePointsClose(_ point1: CGPoint, point2: CGPoint, threshold: CGFloat = ThresholdConstants.pointCloseness) -> Bool {
        return areValuesClose(point1.x, value2: point2.x, threshold: threshold)
            && areValuesClose(point1.y, value2: point2.y, threshold: threshold);
    }
    
    static func areValuesClose<N: BinaryFloatingPoint>(_ value1: N, value2: N, threshold: N = ThresholdConstants.pointCloseness) -> Bool {
        let delta = value1 - value2
        return (delta <= threshold) && (delta >= -threshold)
    }
}

// MARK: - Angle Helpers

/// Helper methods for angles.
enum AngleMath {
    static let twoPi = 2.0 * Double.pi
    static let pi = Double.pi
    static let halfPi = Double.pi / 2
    
    // Normalize the angle between 0 and 2 π
    static func normalizeAngle(_ value: Double) -> Double {
        var value = value
        while value < 0.0 {  value = value + twoPi }
        while value >= twoPi { value = value - twoPi }
        
        return value
    }
    
    // Compute the polar angle from the cartesian point
    static func polarAngle(_ point: CGPoint) -> Double {
        var value = 0.0
        let dpx = Double(point.x)
        let dpy = Double(point.y)
        
        if point.x > 0.0 {
            value = atan(dpy / dpx)
        } else if point.x < 0.0 {
            if point.y >= 0.0 {
                value = atan(dpy / dpx) + pi
            } else {
                value = atan(dpy / dpx) - pi
            }
        } else {
            if point.y > 0.0 {
                value = halfPi
            } else if point.y < 0.0 {
                value = -halfPi
            } else {
                value = 0.0
            }
        }
        
        return normalizeAngle(value)
    }
}

// MARK: - Comparison Range Helpers

enum ComparisonMath {
    static func isValueGreaterThan<N: BinaryFloatingPoint>(_ value: N, minimum: N, threshold: N = ThresholdConstants.tangentCloseness) -> Bool {
        if ProximityMath.areValuesClose(value, value2: minimum, threshold: threshold) {
            return false
        }
        return value > minimum
    }
    
    static func isValueLessThan<N: BinaryFloatingPoint>(_ value: N, maximum: N, threshold: N = ThresholdConstants.tangentCloseness) -> Bool {
        if ProximityMath.areValuesClose(value, value2: maximum, threshold: threshold) {
            return false
        }
        return value < maximum
    }
    
    static func isValueGreaterThanOrEqual<N: BinaryFloatingPoint>(_ value: N, minimum: N, threshold: N = ThresholdConstants.tangentCloseness) -> Bool {
        if ProximityMath.areValuesClose(value, value2: minimum, threshold: threshold) {
            return true
        }
        return value >= minimum
    }
    
    static func isValueLessThanOrEqual<N: BinaryFloatingPoint>(_ value: N, maximum: N, threshold: N = ThresholdConstants.tangentCloseness) -> Bool {
        if ProximityMath.areValuesClose(value, value2: maximum, threshold: threshold) {
            return true
        }
        return value <= maximum
    }
    
    static func angleRangeContainsAngle(_ range: ClosedRange<Double>, angle: Double) -> Bool {
//        return range.contains(angle)
        
        if range.lowerBound <= range.upperBound {
            return ComparisonMath.isValueGreaterThan(angle, minimum: range.lowerBound) && ComparisonMath.isValueLessThan(angle, maximum: range.upperBound)
        }
        
        // The range wraps around 0. See if the angle falls in the first half
        if ComparisonMath.isValueGreaterThan(angle, minimum: range.lowerBound) && angle <= AngleMath.twoPi {
            return true
        }
        
        return angle >= 0.0 && ComparisonMath.isValueLessThan(angle, maximum: range.upperBound)
    }
}

// MARK: Parameter Range Helpers

enum RangeMath {
    static func hasConverged(_ range: ClosedRange<Double>, decimalPlaces: Int) -> Bool {
        let factor = pow(10.0, Double(decimalPlaces))
        let minimum = Int(range.lowerBound * factor)
        let maxiumum = Int(range.upperBound * factor)
        return minimum == maxiumum
    }
    
    static func getSize(_ range: ClosedRange<Double>) -> Double {
        return range.upperBound - range.lowerBound
    }
    
    static func average(_ range: ClosedRange<Double>) -> Double {
        return (range.lowerBound + range.upperBound) / 2.0
    }
    
    static func scaleNormalizedValue(_ range: ClosedRange<Double>, value: Double) -> Double {
        return (range.upperBound - range.lowerBound) * value + range.lowerBound
    }
    
    static func union(_ range1: ClosedRange<Double>, range2: ClosedRange<Double>) -> ClosedRange<Double> {
        return min(range1.lowerBound, range2.lowerBound)...max(range1.upperBound, range2.upperBound)
        //    BPRange(minimum: min(range1.lowerBound, range2.lowerBound),
        //                   maximum: max(range1.upperBound, range2.upperBound))
    }
}

// MARK: - Tangents

struct BPTangentPair {
    var left: CGPoint
    var right: CGPoint
}

struct BPAnglePair {
    var a: Double
    var b: Double
}

// MARK: - Tangent Helpers

enum TangentMath {
    static func areAmbigious(_ edge1Tangents: BPTangentPair, edge2Tangents: BPTangentPair) -> Bool {
        let normalEdge1 = BPTangentPair(left: PointMath.normalizePoint(edge1Tangents.left), right: PointMath.normalizePoint(edge1Tangents.right))
        let normalEdge2 = BPTangentPair(left: PointMath.normalizePoint(edge2Tangents.left), right: PointMath.normalizePoint(edge2Tangents.right))
        
        return ProximityMath.arePointsClose(normalEdge1.left,
                                            point2: normalEdge2.left,
                                            threshold: ThresholdConstants.tangentCloseness)
        
        || ProximityMath.arePointsClose(normalEdge1.left,
                                        point2: normalEdge2.right,
                                        threshold: ThresholdConstants.tangentCloseness)
        
        || ProximityMath.arePointsClose(normalEdge1.right,
                                        point2: normalEdge2.left,
                                        threshold: ThresholdConstants.tangentCloseness)
        
        || ProximityMath.arePointsClose(normalEdge1.right,
                                        point2: normalEdge2.right,
                                        threshold: ThresholdConstants.tangentCloseness)
    }
    
    static func tangentsCross(_ edge1Tangents: BPTangentPair, edge2Tangents: BPTangentPair) -> Bool {
        // Calculate angles for the tangents
        let edge1Angles = BPAnglePair(a: AngleMath.polarAngle(edge1Tangents.left), b: AngleMath.polarAngle(edge1Tangents.right))
        let edge2Angles = BPAnglePair(a: AngleMath.polarAngle(edge2Tangents.left), b: AngleMath.polarAngle(edge2Tangents.right))
        
        // Count how many times edge2 angles appear between the self angles
        let range1 = edge1Angles.a...edge1Angles.b // AngleRange(minimum: edge1Angles.a, maximum: edge1Angles.b)
        var rangeCount1 = 0
        
        if ComparisonMath.angleRangeContainsAngle(range1, angle: edge2Angles.a) {
            rangeCount1 += 1
        }
        if ComparisonMath.angleRangeContainsAngle(range1, angle: edge2Angles.b) {
            rangeCount1 += 1
        }
        
        // Count how many times self angles appear between the edge2 angles
        let range2 = edge1Angles.b...edge1Angles.a // AngleRange(minimum: edge1Angles.b, maximum: edge1Angles.a)
        var rangeCount2 = 0
        
        if ComparisonMath.angleRangeContainsAngle(range2, angle: edge2Angles.a) {
            rangeCount2 += 1
        }
        if ComparisonMath.angleRangeContainsAngle(range2, angle: edge2Angles.b) {
            rangeCount2 += 1
        }
        
        // If each pair of angles split the other two, then the edges cross.
        return rangeCount1 == 1 && rangeCount2 == 1
    }
    
    
    static func lineBoundsMightOverlap(_ bounds1: CGRect, bounds2: CGRect) -> Bool {
        let left = Double(max(bounds1.minX, bounds2.minX))
        let right = Double(min(bounds1.maxX, bounds2.maxX))
        
        if ComparisonMath.isValueGreaterThan(left, minimum: right, threshold: ThresholdConstants.boundsCloseness) {
            return false    // no horizontal overlap
        }
        
        let top = Double(max(bounds1.minY, bounds2.minY))
        let bottom = Double(min(bounds1.maxY, bounds2.maxY))
        return ComparisonMath.isValueLessThanOrEqual(top, maximum: bottom, threshold: ThresholdConstants.boundsCloseness)
    }
}
