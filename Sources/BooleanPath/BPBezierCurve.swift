//
//  BPBezierCurve.swift
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

// FBBezierCurve is one cubic 2D bezier curve.
// It represents one segment of a bezier path, and is where
// the intersection calculation happens

import SwiftUI
import CoreGraphics

struct BPNormalizedLine {
    var a: Double // * x +
    var b: Double // * y +
    var c: Double // constant
    
    init(a: Double, b: Double, c: Double) {
        self.a = a
        self.b = b
        self.c = c
    }
    
    // Create a normalized line such that computing the distance from it is quick.
    //  See:    http://softsurfer.com/Archive/algorithm_0102/algorithm_0102.htm#Distance%20to%20an%20Infinite%20Line
    //          http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/geometry/basic.html
    //
    init(point1: CGPoint, point2: CGPoint) {
        self.a = Double(point1.y - point2.y)
        self.b = Double(point2.x - point1.x)
        self.c = Double(point1.x * point2.y - point2.x * point1.y)
        
        let distance = sqrt(self.b * self.b + self.a * self.a)
        
        // GPC: prevent divide-by-zero from putting NaNs into the values which cause trouble further on. I'm not sure
        // what cases trigger this, but sometimes point1 == point2 so distance is 0.
        if distance != 0.0 {
            self.a /= distance
            self.b /= distance
            self.c /= distance
        } else {
            self.a = 0
            self.b = 0
            self.c = 0
        }
    }
    
    func copyWithOffset(_ offset: Double) -> BPNormalizedLine {
        return BPNormalizedLine(
            a: self.a,
            b: self.b,
            c: self.c + offset)
    }
    
    func distanceFromPoint(_ point: CGPoint) -> Double {
        return a * Double(point.x) + b * Double(point.y) + c;
    }
    
    func intersectionWith(_ other: BPNormalizedLine) -> CGPoint {
        let denominator = (self.a * other.b) - (other.a * self.b)
        
        return CGPoint(
            x: (self.b * other.c - other.b * self.c) / denominator,
            y: (self.a * other.c - other.a * self.c) / denominator)
    }
}

// MARK: - Helper Functions

enum CurveHelpers {
    static func parameterOfPointOnLine(_ lineStart: CGPoint, lineEnd: CGPoint, point: CGPoint) -> Double {
        
    // Note: its assumed you have already checked that point is colinear with the line (lineStart, lineEnd)
        let lineLength: Double = PointMath.distanceBetweenPoints(lineStart, point2: lineEnd)
        let lengthFromStart: Double = PointMath.distanceBetweenPoints(point, point2: lineStart)
        var parameter = lengthFromStart / lineLength
        
    // The only tricky thing here is the sign. Is the point _before_ lineStart, or after lineStart?
        let lengthFromEnd: Double = PointMath.distanceBetweenPoints(point, point2: lineEnd)
        if ProximityMath.areValuesClose(lineLength + lengthFromStart, value2: lengthFromEnd) {
            parameter = -parameter
        }
        return parameter
    }
    
    static func linesIntersect(_ line1Start: CGPoint, line1End: CGPoint, line2Start: CGPoint, line2End: CGPoint, outIntersect: inout CGPoint) -> Bool {
        let line1 = BPNormalizedLine(point1: line1Start, point2: line1End)
        let line2 = BPNormalizedLine(point1: line2Start, point2: line2End)
        outIntersect = line1.intersectionWith(line2)
        if outIntersect.x.isNaN || outIntersect.y.isNaN {
            return false
        }
        outIntersect.y = -outIntersect.y
        return true
    }
/// The three points are a counter-clockwise turn if the return value is greater than 0,
///  clockwise if less than 0, or colinear if 0.    
    static func counterClockwiseTurn(_ point1: CGPoint, point2: CGPoint, point3: CGPoint) -> Double {
        // We're calculating the signed area of the triangle formed by the three points. Well,
        // almost the area of the triangle -- we'd need to divide by 2. But since we only
        // care about the direction (i.e. the sign) dividing by 2 is an unnecessary step.
        // See http://mathworld.wolfram.com/TriangleArea.html for the signed area of a triangle.
        
        let xDeltaA = Double(point2.x - point1.x)
        let yDeltaB = Double(point3.y - point1.y)
        let yDeltaC = Double(point2.y - point1.y)
        let xDeltaD = Double(point3.x - point1.x)
        return xDeltaA * yDeltaB - yDeltaC * xDeltaD
    }
    
    /// Calculate if and where the given line intersects the horizontal line at y.
    static func lineIntersectsHorizontalLine(_ startPoint: CGPoint, endPoint: CGPoint, y: Double, intersectPoint: inout CGPoint) -> Bool {
        // Do a quick test to see if y even falls on the startPoint,endPoint line
        let minY = Double(min(startPoint.y, endPoint.y))
        let maxY = Double(max(startPoint.y, endPoint.y))
        if (y < minY && !ProximityMath.areValuesClose(y, value2: minY)) || (y > maxY && !ProximityMath.areValuesClose(y, value2: maxY)) {
            return false
        }
        
        // There's an intersection here somewhere
        if startPoint.x == endPoint.x {
            intersectPoint = CGPoint(x: startPoint.x, y: CGFloat(y))
        }
        else {
            let slope = Double(endPoint.y - startPoint.y) / Double(endPoint.x - startPoint.x)
            intersectPoint = CGPoint(
                x: CGFloat((y - Double(startPoint.y)) / slope) + startPoint.x,
                y: CGFloat(y))
        }
        return true
    }
    /// Calculate a point on the bezier curve passed in, specifically the point at parameter.
    
    static func bezierWithPoints(_ degree: Int, bezierPoints: [CGPoint], parameter: Double, withCurves: Bool) -> (point: CGPoint, leftCurve: [CGPoint]?, rightCurve: [CGPoint]?) {
        // We're using De Casteljau's algorithm, which not only calculates the point at parameter
        // in a numerically stable way, it also computes the two resulting bezier curves that
        // would be formed if the original were split at the parameter specified.
        //
        // See: http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/de-casteljau.html
        // for an explaination of De Casteljau's algorithm.
    
        // bezierPoints, leftCurve, rightCurve will have a length of degree + 1.
        // degree is the order of the bezier path, which will be cubic (3) most of the time.
    
        // With this algorithm we start out with the points in the bezier path.
        // we assume we'll never get more than a cubic bezier
        var points: [CGPoint] = []
        
        for i in 0 ... degree {
            points.append(bezierPoints[i])
        }
        
        var leftArray = [CGPoint](repeating: CGPoint.zero, count: degree+1)
        var rightArray = [CGPoint](repeating: CGPoint.zero, count: degree+1)
        
        // If the caller is asking for the resulting bezier curves, start filling those in
        if withCurves {
            leftArray[0] = points[0]
            rightArray[degree] = points[degree]
        }
        
        for k in 1 ... degree {
            for i in 0 ... (degree - k) {
                let pxV = (1.0 - parameter) * Double(points[i].x) + parameter * Double(points[i + 1].x)
                points[i].x = CGFloat(pxV)
                let pyV = (1.0 - parameter) * Double(points[i].y) + parameter * Double(points[i + 1].y)
                points[i].y = CGFloat(pyV)
            }
            if withCurves {
                leftArray[k] = points[0]
                rightArray[degree-k] = points[degree-k]
            }
        }
        
    // The point in the curve at parameter ends up in points[0]
        return (point: points[0], leftCurve: leftArray, rightCurve: rightArray)
    }
    
    static func computeCubicFirstDerivativeRoots(_ a: Double, b: Double, c: Double, d: Double) -> [Double] {
        // See http://processingjs.nihongoresources.com/bezierinfo/#bounds for where the formulas come from
        let denominator = -a + 3.0 * b - 3.0 * c + d
        
    // If denominator == 0, fall back to
        if ProximityMath.areValuesClose(denominator, value2: 0.0) {
            let t = (a - b) / (2.0 * (a - 2.0 * b + c))
            return [t];
        }
        
        let numeratorLeft = -a + 2.0 * b - c
        
        let v1 = -a * (c - d)
        let v2 = b * b
        let v3 = b * (c + d)
        let v4 = c * c
        let numeratorRight = -1.0 * sqrt(v1 + v2 - v3 + v4)
        
        let t1 = (numeratorLeft + numeratorRight) / denominator
        let t2 = (numeratorLeft - numeratorRight) / denominator
        // NOTE: This should be changed to return a tuple
        return [t1, t2]
    }
    
    static func gaussQuadratureBaseForCubic(_ t: Double, p1: Double, p2: Double, p3: Double, p4: Double) -> Double {
        let t1 = (-3.0 * p1) + (9.0 * p2) - (9.0 * p3) + (3.0 * p4)
        let t2 = t * t1 + 6.0 * p1 - 12.0 * p2 + 6.0 * p3
        
        return t * t2 - 3.0 * p1 + 3.0 * p2
        //return t * (t * (-3 * p1 + 9 * p2 - 9 * p3 + 3 * p4) + 6 * p1 + 12 * p2 + 3 * p3) - 3 * p1 + 3 * p2
    }
    
    static func gaussQuadratureFOfTForCubic(_ t: Double, p1: CGPoint, p2: CGPoint, p3: CGPoint, p4: CGPoint) -> Double {
        let baseX = CurveHelpers.gaussQuadratureBaseForCubic(t,
                                                             p1: Double(p1.x),
                                                             p2: Double(p2.x),
                                                             p3: Double(p3.x),
                                                             p4: Double(p4.x))
        let baseY = CurveHelpers.gaussQuadratureBaseForCubic(t,
                                                             p1: Double(p1.y),
                                                             p2: Double(p2.y),
                                                             p3: Double(p3.y),
                                                             p4: Double(p4.y))
        
        return sqrt(baseX * baseX + baseY * baseY)
    }
    
    static func gaussQuadratureComputeCurveLengthForCubic(_ z: Double, steps: Int, p1: CGPoint, p2: CGPoint, p3: CGPoint, p4: CGPoint) -> Double {
        let z2 = Double(z / 2.0)
        var sum: Double = 0.0
        for i in 0 ..< steps {
            let correctedT: Double = z2 * BPLegendreGaussAbscissaeValues[steps][i] + z2
            sum += BPLegendreGaussWeightValues[steps][i] * CurveHelpers.gaussQuadratureFOfTForCubic(correctedT, p1: p1, p2: p2, p3: p3, p4: p4)
        }
        return z2 * sum
    }
    
    static func sign(_ value: CGFloat) -> Int {
        if value < 0.0 {
            return -1
        } else {
            return 1
        }
    }
    
    static func countBezierCrossings(_ bezierPoints: [CGPoint], degree: Int) -> Int {
        var count: Int = 0
        var sign = CurveHelpers.sign(bezierPoints[0].y)
        
        var previousSign = sign
        for i in 1...degree {
            sign = CurveHelpers.sign(bezierPoints[i].y)
            if sign != previousSign {
                count += 1
            }
            previousSign = sign;
        }
        return count
    }
    
    static let findBezierRootsMaximumDepth = 64
    
    /// FBIsControlPolygonFlatEnough Usage:
    /// var intersectAt = CGPoint(x:0,y:0)
    /// var boolGotit = FBIsControlPolygonFlatEnough(points,degree,&intersectAt)
    static func isControlPolygonFlatEnough(_ bezierPoints: [CGPoint], degree: Int, intersectionPoint: inout CGPoint) -> Bool {
        
        // 2^-63
        let BPFindBezierRootsErrorThreshold = CGFloat(ldexpf(Float(1.0), Int32(-1 * (CurveHelpers.findBezierRootsMaximumDepth - 1))))
        let line = BPNormalizedLine(point1: bezierPoints[0], point2: bezierPoints[degree])
        
        // Find the bounds around the line
        var belowDistance = 0.0
        var aboveDistance = 0.0
        for i in 1 ..< degree {
            let distance = line.distanceFromPoint(bezierPoints[i])
            if distance > aboveDistance {
                aboveDistance = distance
            }
            if distance < belowDistance {
                belowDistance = distance
            }
        }
        
        let zeroLine = BPNormalizedLine(a: 0.0, b: 1.0, c: 0.0)
        let aboveLine = line.copyWithOffset(-aboveDistance)
        let intersect1 = zeroLine.intersectionWith(aboveLine)
        
        let belowLine = line.copyWithOffset(-belowDistance)
        let intersect2 = zeroLine.intersectionWith(belowLine)
        
        let error = max(intersect1.x, intersect2.x) - min(intersect1.x, intersect2.x)
        if error < BPFindBezierRootsErrorThreshold {
            intersectionPoint = zeroLine.intersectionWith(line)
            return true
        }
        return false
    }
    
    static func findBezierRootsWithDepth(_ bezierPoints: [CGPoint], degree: Int, depth: Int, perform: (_ root: Double) -> Void) {
        let crossingCount = CurveHelpers.countBezierCrossings(bezierPoints, degree: degree)
        if crossingCount == 0 {
            return
        }
        else if crossingCount == 1 {
            if depth >= CurveHelpers.findBezierRootsMaximumDepth {
                let root = Double(bezierPoints[0].x + bezierPoints[degree].x) / 2.0
                perform(root)
                return
            }
            var intersectionPoint = CGPoint.zero
            if CurveHelpers.isControlPolygonFlatEnough(bezierPoints, degree: degree, intersectionPoint: &intersectionPoint) {
                perform(Double(intersectionPoint.x))
                return
            }
        }
        
        // Subdivide and try again
        let bwp = CurveHelpers.bezierWithPoints(degree, bezierPoints: bezierPoints, parameter: 0.5, withCurves: true)
        CurveHelpers.findBezierRootsWithDepth(bwp.leftCurve!, degree: degree, depth: depth + 1, perform: perform)
        CurveHelpers.findBezierRootsWithDepth(bwp.rightCurve!, degree: degree, depth: depth + 1, perform: perform)
    }
    
    static func findBezierRoots(_ bezierPoints: [CGPoint], degree: Int, perform: (_ root: Double) -> Void) {
        CurveHelpers.findBezierRootsWithDepth(bezierPoints, degree: degree, depth: 0, perform: perform)
    }
}

// MARK: - Convex Hull

enum ConvexHullHelpers {
    static func doPointsTurnWrongDirection(_ point1: CGPoint, point2: CGPoint, point3: CGPoint) -> Bool {
        let area = CurveHelpers.counterClockwiseTurn(point1, point2: point2, point3: point3)
        return ProximityMath.areValuesClose(area, value2: 0.0) || area < 0.0
    }
    
    static func buildFromPoints(_ inPoints: [CGPoint]) -> (hull: [CGPoint], hullLength: Int) {
        // Compute the convex hull for this bezier curve. The convex hull is made up of the end and control points.
        //  The hard part is determine the order they go in, and if any are inside or colinear with the convex hull.
    
        // Uses the Monotone chain algorithm:
        //  http://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
    
        // Start with all the end and control points in any order.
        let numberOfPoints = 4
        var points = inPoints
        
        // Sort points ascending x, if equal compare y
        //  Bubble sort, which should be ok with a max of 4 elements, and the fact that our one current use case
        //  already has them in ascending X order (i.e. should be just comparisons to verify)
        var sortLength = numberOfPoints
        repeat {
            var newSortLength = 0
            for i in 1 ..< sortLength {
                if ( points[i - 1].x > points[i].x || (ProximityMath.areValuesClose(points[i - 1].x, value2: points[i].x) && points[i - 1].y > points[i].y) ) {
                    let tempPoint = points[i]
                    points[i] = points[i - 1]
                    points[i - 1] = tempPoint
                    newSortLength = i
                }
            }
            sortLength = newSortLength
        } while sortLength > 0
        
        // Create the results
        var filledInIx = 0
        var results = [CGPoint](repeating: CGPoint.zero, count: 8)
        
        // Build lower hull
        for i in 0..<numberOfPoints {
            while filledInIx >= 2 && ConvexHullHelpers.doPointsTurnWrongDirection(results[filledInIx - 2], point2: results[filledInIx - 1], point3: points[i]) {
                filledInIx -= 1
            }
            results[filledInIx] = points[i];
            filledInIx += 1;
        }
        
        // Build upper hull
        let thresholdIndex = filledInIx + 1
        for i in (0 ... numberOfPoints - 2).reversed() {
            while filledInIx >= thresholdIndex && ConvexHullHelpers.doPointsTurnWrongDirection(results[filledInIx - 2], point2: results[filledInIx - 1], point3: points[i]) {
                filledInIx -= 1
            }
            
            results[filledInIx] = points[i];
            filledInIx += 1;
        }
        
        return (hull: results, hullLength: filledInIx - 1)
    }
    
    // MARK: ILLUSTRATE_CALL_TO_BP
    static func buildFromPoints() {
        // In C this was:
        //  NSUInteger convexHullLength = 0;
        //  NSPoint convexHull[8] = {};
        //  FBConvexHullBuildFromPoints(distanceBezierPoints, convexHull, &convexHullLength);
        let distanceBezierPoints = [CGPoint](repeating: CGPoint.zero, count: 4)
        
        let (_, convexHullLength) = ConvexHullHelpers.buildFromPoints(distanceBezierPoints)
        
        if convexHullLength < 0 {
            print("YIKES")
        }
    }
}

// MARK: - BPBezierCurveData

struct BPBezierCurveLocation {
    var parameter: Double
    var distance: Double
    
    init(parameter: Double, distance: Double) {
        self.parameter = parameter
        self.distance = distance
    }
}

class BPBezierCurveData {
    var endPoint1: CGPoint
    var controlPoint1: CGPoint
    var controlPoint2: CGPoint
    var endPoint2: CGPoint
    
    var isStraightLine : Bool		// GPC: flag when curve came from a straight line segment
    
    var length : Double?         // cached value
    fileprivate var _bounds : CGRect? // cached value
    var _isPoint : Bool?          // cached value
    var _boundingRect : CGRect?   // cached value
    
    // FBBezierCurveDataMake
    // let some_data = FBBezierCurveData(endPoint1, controlPoint1, controlPoint2, endPoint2, isStraightLine)
    //
    init(endPoint1 : CGPoint, controlPoint1 : CGPoint, controlPoint2 : CGPoint, endPoint2 : CGPoint, isStraightLine: Bool) {
        self.endPoint1 = endPoint1
        self.controlPoint1 = controlPoint1
        self.controlPoint2 = controlPoint2
        self.endPoint2 = endPoint2
        self.isStraightLine = isStraightLine
        
        self.length = BPBezierCurveData.invalidLength
    }
    
    init(cloning clone: BPBezierCurveData) {
        self.endPoint1 = clone.endPoint1
        self.controlPoint1 = clone.controlPoint1
        self.controlPoint2 = clone.controlPoint2
        self.endPoint2 = clone.endPoint2
        self.isStraightLine = clone.isStraightLine
        
        self.length = clone.length
    }
    
    /// Instead of using the C language trick of setting a BOOL to a negative number
    /// to see whether a test has ever been performed, we are using a Swift optional
    ///
    /// static const BOOL FBBezierCurveDataInvalidIsPoint = -1;
    /// if _data.isPoint != nil, let isPoint = _data.isPoint {
    ///     if isPoint {
    ///         println("Got it true")
    ///     } else {
    ///         println("Got it false")
    ///     }
    /// } else {
    ///     println("Have to make sure we set it NOW")
    /// }
    var isPoint: Bool {
        // If the two end points are close together, then we're a point.
        // Ignore the control points.
        let BPClosenessThreshold = 1e-5
        
        // check cached value
        if let _isPoint {
            return _isPoint
        }
        
        _isPoint = ProximityMath.arePointsClose(endPoint1, point2: endPoint2, threshold: BPClosenessThreshold)
            && ProximityMath.arePointsClose(endPoint1, point2: controlPoint1, threshold: BPClosenessThreshold)
            && ProximityMath.arePointsClose(endPoint1, point2: controlPoint2, threshold: BPClosenessThreshold)
        
        return _isPoint!;
    }
    
    // static NSRect FBBezierCurveDataBoundingRect(FBBezierCurveData *me)
    // Becomes:
    // let bnds = curve_data.boundingRect()
    var boundingRect: CGRect {
        // Use the cache if we have one
        if let _boundingRect {
            return _boundingRect
        }
        
        let left = min(endPoint1.x, min(controlPoint1.x, min(controlPoint2.x, endPoint2.x)))
        let top = min(endPoint1.y, min(controlPoint1.y, min(controlPoint2.y, endPoint2.y)))
        let right = max(endPoint1.x, max(controlPoint1.x, max(controlPoint2.x, endPoint2.x)))
        let bottom = max(endPoint1.y, max(controlPoint1.y, max(controlPoint2.y, endPoint2.y)))
        
        _boundingRect = CGRect(x: left, y: top, width: right - left, height: bottom - top)
        
        return _boundingRect!
    }
    
    // static NSRect FBBezierCurveDataBounds(FBBezierCurveData* me)
    // Becomes:
    // let bnds = curve_data.bounds()
    var bounds: CGRect {
        // Use the cache if we have one
        if let _bounds {
            return _bounds
        }
        
        var bounds = CGRect.zero
        
        if isStraightLine {
            var topLeft = endPoint1
            var bottomRight = topLeft
            PointMath.expandBoundsByPoint(&topLeft, bottomRight: &bottomRight, point: endPoint2)
            
            bounds = CGRect(x: topLeft.x, y: topLeft.y, width: bottomRight.x - topLeft.x, height: bottomRight.y - topLeft.y)
        } else {
            // Start with the end points
            var (topLeft, _, _) = pointAtParameter(0)
            var bottomRight = topLeft
            let (lastPoint, _, _) = pointAtParameter(1)
            
            PointMath.expandBoundsByPoint(&topLeft, bottomRight: &bottomRight, point: lastPoint);
            
            // Find the roots, which should be the extremities
            let xRoots: [Double] = CurveHelpers.computeCubicFirstDerivativeRoots(Double(endPoint1.x), b: Double(controlPoint1.x), c: Double(controlPoint2.x), d: Double(endPoint2.x))
            
            for i in 0 ..< xRoots.count {
                let t = xRoots[i]
                if t < 0 || t > 1 {
                    continue
                }
                let (location, _, _) = pointAtParameter(t)
                PointMath.expandBoundsByPoint(&topLeft, bottomRight: &bottomRight, point: location)
            }
            
            let yRoots: [Double] = CurveHelpers.computeCubicFirstDerivativeRoots(Double(endPoint1.y), b: Double(controlPoint1.y), c: Double(controlPoint2.y), d: Double(endPoint2.y))
            for i in 0 ..< yRoots.count {
                let t = yRoots[i]
                if t < 0 || t > 1 {
                    continue
                }
                let (location, _, _) = pointAtParameter(t)
                PointMath.expandBoundsByPoint(&topLeft, bottomRight: &bottomRight, point: location)
            }
            
            bounds = CGRect(x: topLeft.x, y: topLeft.y, width: bottomRight.x - topLeft.x, height: bottomRight.y - topLeft.y)
        }
        
        _bounds = bounds
        return bounds
    }
    
    // static CGFloat FBBezierCurveDataGetLengthAtParameter(FBBezierCurveData* me, CGFloat parameter)
    // Becomes:
    // let someLength = curve_data.getLengthAtParameter(CGFloat parameter)
    //
    func getLengthAtParameter(_ parameter: Double) -> Double {
        
        // Use the cached value if at all possible
        if parameter == 1.0 && length != BPBezierCurveData.invalidLength,
            let length {
            return length
        }
        
        // If it's a line, use that equation instead
        var calculatedLength = BPBezierCurveData.invalidLength
        if isStraightLine {
            calculatedLength = PointMath.distanceBetweenPoints(endPoint1, point2: endPoint2) * parameter
        } else {
            calculatedLength = CurveHelpers.gaussQuadratureComputeCurveLengthForCubic(Double(parameter), steps: 12, p1: endPoint1, p2: controlPoint1, p3: controlPoint2, p4: endPoint2)
        }
        
        // If possible, update our cache
        if parameter == 1.0 {
            length = calculatedLength
        }
        
        return calculatedLength
    }
    
    func reversed() -> BPBezierCurveData {
        return BPBezierCurveData(endPoint1: endPoint2, controlPoint1: controlPoint2, controlPoint2: controlPoint1, endPoint2: endPoint1, isStraightLine: isStraightLine)
    }
    
    // static CGFloat FBBezierCurveDataGetLength(FBBezierCurveData* me)
    // Becomes:
    // let someLength = curve_data.getLength()
    func getLength() -> Double {
        return getLengthAtParameter(1.0)
    }
    
    // NSPoint FBBezierCurveDataPointAtParameter(FBBezierCurveData me, CGFloat param, FBB-C-Data *leftBezCurve, FBB-C-Data *rightBezCurve)
    // Becomes:
    // let (point,left,right) = curve_data.pointAtParameter(param_float)
    func pointAtParameter(_ parameter: Double) -> (point: CGPoint, leftCurve: BPBezierCurveData?, rightCurve: BPBezierCurveData?) {
        // This method is a simple wrapper around the BezierWithPoints() helper function. It computes the 2D point at the given parameter,
        //  and (optionally) the resulting curves that splitting at the parameter would create.
        let points = [endPoint1, controlPoint1, controlPoint2, endPoint2]
        
        let bwp = CurveHelpers.bezierWithPoints(3, bezierPoints: points, parameter: parameter, withCurves: true);
        
        var leftBCData: BPBezierCurveData? = nil
        if let leftA = bwp.leftCurve {
            leftBCData = BPBezierCurveData(endPoint1: leftA[0], controlPoint1: leftA[1], controlPoint2: leftA[2], endPoint2: leftA[3], isStraightLine: isStraightLine)
        }
        
        var rightBCData: BPBezierCurveData? = nil
        if let rightA = bwp.rightCurve {
            rightBCData = BPBezierCurveData(endPoint1: rightA[0], controlPoint1: rightA[1], controlPoint2: rightA[2], endPoint2: rightA[3], isStraightLine: isStraightLine)
        }
        
        return (point: bwp.point, leftCurve: leftBCData, rightCurve: rightBCData)
    }
    
    // static FBBezierCurveData FBBezierCurveDataSubcurveWithRange(FBBezierCurveData me, FBRange range)
    // Becomes:
    // let new_curve_data = curve_data.subcurveWithRange(range)
    func subcurveWithRange(_ range: ClosedRange<Double>) -> BPBezierCurveData {
        let upperCurve = self.pointAtParameter(range.lowerBound).rightCurve!
        
        if range.lowerBound == 1.0 {
            return upperCurve           // avoid the divide by zero below
        }
        
        // We need to adjust the maximum parameter to fit on the new curve before we split again
        let adjustedMaximum = (range.upperBound - range.lowerBound) / (1.0 - range.lowerBound)
        
        let lowerCurve = upperCurve.pointAtParameter(adjustedMaximum).leftCurve!
        
        return lowerCurve
    }
    
    // static FBNormalizedLine FBBezierCurveDataRegularFatLineBounds(FBBezierCurveData me, FBRange *range)
    // Becomes:
    // let normalized_line = curve_data.regularFatLineBounds(&range)
    func regularFatLineBounds() -> (line: BPNormalizedLine, range: ClosedRange<Double>) {
        // Create the fat line based on the end points
        let line = BPNormalizedLine(point1: endPoint1, point2: endPoint2)
        
        // Compute the bounds of the fat line. The fat line bounds should entirely encompass the
        //  bezier curve. Since we know the convex hull entirely compasses the curve, just take
        //  all four points that define this cubic bezier curve. Compute the signed distances of
        //  each of the end and control points from the fat line, and that will give us the bounds.
        
        // In this case we know that the end points are on the line, thus their distances will be 0.
        //  So we can skip computing those and just use 0.
        let controlPoint1Distance = line.distanceFromPoint(controlPoint1)
        let controlPoint2Distance = line.distanceFromPoint(controlPoint2)
        
        let minim = min(controlPoint1Distance, min(controlPoint2Distance, 0.0))
        let maxim = max(controlPoint1Distance, max(controlPoint2Distance, 0.0))
        
        return (line, minim...maxim)
    }
    
    // static FBNormalizedLine FBBezierCurveDataPerpendicularFatLineBounds(FBBezierCurveData me, FBRange *range)
    // Becomes:
    // let normalized_line = curve_data.perpendicularFatLineBounds(&range)
    
    func perpendicularFatLineBounds() -> (line: BPNormalizedLine, range: ClosedRange<Double>) {
        
        // Create a fat line that's perpendicular to the line created by the two end points.
        let normal = PointMath.lineNormal(endPoint1, lineEnd: endPoint2)
        let startPoint = PointMath.lineMidpoint(endPoint1, lineEnd: endPoint2)
        let endPoint = PointMath.addPoint(startPoint, point2: normal)
        let line = BPNormalizedLine(point1: startPoint, point2: endPoint)
        
        // Compute the bounds of the fat line. The fat line bounds should entirely encompass the
        //  bezier curve. Since we know the convex hull entirely compasses the curve, just take
        //  all four points that define this cubic bezier curve. Compute the signed distances of
        //  each of the end and control points from the fat line, and that will give us the bounds.
        let controlPoint1Distance = line.distanceFromPoint(controlPoint1)
        let controlPoint2Distance = line.distanceFromPoint(controlPoint2)
        let point1Distance = line.distanceFromPoint(endPoint1)
        let point2Distance = line.distanceFromPoint(endPoint2)
        
        let minim = min(controlPoint1Distance, min(controlPoint2Distance, min(point1Distance, point2Distance)))
        let maxim = max(controlPoint1Distance, max(controlPoint2Distance, max(point1Distance, point2Distance)))
        
        return (line, minim...maxim)
    }
    
    // static FBRange FBBezierCurveDataClipWithFatLine(FBBezierCurveData me, FBNormalizedLine fatLine, FBRange bounds)
    // Becomes:
    // let new_range = curve_data.clipWithFatLine(fatline, bounds)
    func clipWithFatLine(_ fatLine: BPNormalizedLine, bounds: ClosedRange<Double>) -> ClosedRange<Double> {
        // This method computes the range of self that could possibly intersect with
        // the fat line passed in (and thus with the curve enclosed by the fat line).
        //
        // To do that, we first compute the signed distance of all our points (end and control)
        // from the fat line, and map those onto a bezier curve at evenly spaced intervals
        // from [0..1]. The parts of the distance bezier that fall inside of the fat line bounds
        // correspond to the parts of ourself that could potentially intersect with the other curve.
        // Ideally, we'd calculate where the distance bezier intersected the horizontal lines
        // representing the fat line bounds. However, computing those intersections is hard and
        // costly. So instead we'll compute the convex hull, and intersect those lines with the
        // fat line bounds. The intersection with the lowest x coordinate will be the minimum,
        // and the intersection with the highest x coordinate will be the maximum.
        
        // The convex hull (for cubic beziers) is the four points that define the curve.
        // A useful property of the convex hull is that the entire curve lies inside of it.
        
        // First calculate bezier curve points distance from the fat line that's clipping us
        let distanceBezierPoints: [CGPoint] = [
            CGPoint(x: 0.0, y: fatLine.distanceFromPoint(endPoint1)),
            CGPoint(x: 1.0 / 3.0, y: fatLine.distanceFromPoint(controlPoint1)),
            CGPoint(x: 2.0 / 3.0, y: fatLine.distanceFromPoint(controlPoint2)),
            CGPoint(x: 1.0, y: fatLine.distanceFromPoint(endPoint2))
        ]
        
        let (convexHull, convexHullLength) = ConvexHullHelpers.buildFromPoints(distanceBezierPoints)
        
        // Find intersections of convex hull with the fat line bounds
        var range = 1.0...0.0
        
        for i in 0 ..< convexHullLength {
            // Pull out the current line on the convex hull
            let indexOfNext = i < (convexHullLength - 1) ? i + 1 : 0
            
            let startPoint = convexHull[i]
            let endPoint = convexHull[indexOfNext]
            var intersectionPoint = CGPoint.zero
            
            // See if the segment of the convex hull intersects with the minimum fat line bounds
            if CurveHelpers.lineIntersectsHorizontalLine(startPoint, endPoint: endPoint, y: bounds.lowerBound, intersectPoint: &intersectionPoint) {
                if Double(intersectionPoint.x) < range.lowerBound {
                    range = Double(intersectionPoint.x)...range.upperBound
//                    range.lowerBound = Double(intersectionPoint.x)
                }
                if Double(intersectionPoint.x) > range.upperBound {
                    range = range.lowerBound...Double(intersectionPoint.x)
//                    range.upperBound = Double(intersectionPoint.x)
                }
            }
            
            // See if this segment of the convex hull intersects with the maximum fat line bounds
            if CurveHelpers.lineIntersectsHorizontalLine(startPoint, endPoint: endPoint, y: bounds.upperBound, intersectPoint: &intersectionPoint) {
                if Double(intersectionPoint.x) < range.lowerBound {
                    range = Double(intersectionPoint.x)...range.upperBound
//                    range.lowerBound = Double(intersectionPoint.x)
                }
                if Double(intersectionPoint.x) > range.upperBound {
                    range = range.lowerBound...Double(intersectionPoint.x)
//                    range.upperBound = Double(intersectionPoint.x)
                }
            }
            
            // We want to be able to refine t even if the convex hull lies completely inside the bounds.
            // This also allows us to be able to use range of [1..0] as a sentinel value meaning the
            // convex hull lies entirely outside of bounds, and the curves don't intersect.
            if Double(startPoint.y) < bounds.upperBound && Double(startPoint.y) > bounds.lowerBound {
                if Double(startPoint.x) < range.lowerBound {
                    range = Double(startPoint.x)...range.upperBound
//                    range.lowerBound = Double(startPoint.x)
                }
                if Double(startPoint.x) > range.upperBound {
                    range = range.lowerBound...Double(startPoint.x)
//                    range.upperBound = Double(startPoint.x)
                }
            }
        }
        
        // Check for bad values
        if range.lowerBound.isInfinite || range.lowerBound.isNaN || range.upperBound.isInfinite || range.upperBound.isNaN {
            // equivalent to: something went wrong, so I don't know
            range = 0...1
        }
        
        return range
    }
    
    func convertSelfAndPoint(_ point: CGPoint) -> [CGPoint] {
        let selfPoints: [CGPoint] = [endPoint1, controlPoint1, controlPoint2, endPoint2]
        
        // c[i] in the paper
        let distanceFromPoint = [
            PointMath.subtractPoint(selfPoints[0], point2: point),
            PointMath.subtractPoint(selfPoints[1], point2: point),
            PointMath.subtractPoint(selfPoints[2], point2: point),
            PointMath.subtractPoint(selfPoints[3], point2: point)
        ]
        
        // d[i] in the paper
        let weightedDelta = [
            PointMath.scalePoint(PointMath.subtractPoint(selfPoints[1], point2: selfPoints[0]), scale: 3),
            PointMath.scalePoint(PointMath.subtractPoint(selfPoints[2], point2: selfPoints[1]), scale: 3),
            PointMath.scalePoint(PointMath.subtractPoint(selfPoints[3], point2: selfPoints[2]), scale: 3)
        ]
        
        // Precompute the dot product of distanceFromPoint and weightedDelta in order to speed things up
        var precomputedTable: [[Double]] = [
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ]
        for row in 0 ..< 3 {
            for column in 0 ..< 4 {
                precomputedTable[row][column] = PointMath.dotMultiplyPoint(weightedDelta[row], point2: distanceFromPoint[column])
            }
        }
        
        // Precompute some of the values to speed things up
        let BPZ: [[Double]] = [
            [1.0, 0.6, 0.3, 0.1],
            [0.4, 0.6, 0.6, 0.4],
            [0.1, 0.3, 0.6, 1.0]
        ]
        
        // create our output array
        var bezierPoints = [CGPoint](repeating: CGPoint.zero, count: 6)
        
        // Set the x values of the bezier points
        for i in 0 ..< 6 {
            bezierPoints[i] = CGPoint(x: CGFloat(i) / 5.0, y: 0)
        }
        
        // Finally set the y values of the bezier points
        let n = 3
        let m = n - 1
        for k in 0 ... (n + m) {
            let lowerBound = max(0, k - m)
            let upperBound = min(k, n)
            for i in lowerBound ... upperBound {
                let j = k - i
                bezierPoints[i + j].y += CGFloat(precomputedTable[j][i] * BPZ[j][i])
            }
        }
        
        return bezierPoints
    }
    
    func closestLocationToPoint(_ point: CGPoint) -> BPBezierCurveLocation {
        let bezierPoints = convertSelfAndPoint(point)
        
        var distance = PointMath.distanceBetweenPoints(endPoint1, point2: point)
        var parameter = 0.0
        
        CurveHelpers.findBezierRoots(bezierPoints, degree: 5) { root in
            let location = self.pointAtParameter(root).point
            let theDistance = PointMath.distanceBetweenPoints(location, point2: point)
            if theDistance < distance {
                distance = theDistance
                parameter = root
            }
        }
        
        let lastDistance = PointMath.distanceBetweenPoints(endPoint2, point2: point)
        if lastDistance < distance {
            distance = lastDistance
            parameter = 1.0
        }
        
        let location = BPBezierCurveLocation(parameter: parameter, distance: distance)
        
        return location
    }
    
    // 893
    //static BOOL FBBezierCurveDataIsEqualWithOptions(FBBezierCurveData me, FBBezierCurveData other, CGFloat threshold)
    func isEqual(to other: BPBezierCurveData, threshold: Double) -> Bool {
        guard !(isPoint || other.isPoint)
                || isStraightLine == other.isStraightLine else { return false }
        
        if isStraightLine {
            return ProximityMath.arePointsClose(endPoint1, point2: other.endPoint1, threshold: threshold)
                && ProximityMath.arePointsClose(endPoint2, point2: other.endPoint2, threshold: threshold)
        } else {
            return ProximityMath.arePointsClose(endPoint1, point2: other.endPoint1, threshold: threshold)
                && ProximityMath.arePointsClose(controlPoint1, point2: other.controlPoint1, threshold: threshold)
                && ProximityMath.arePointsClose(controlPoint2, point2: other.controlPoint2, threshold: threshold)
                && ProximityMath.arePointsClose(endPoint2, point2: other.endPoint2, threshold: threshold)
        }
    }
}

// MARK: Private BPBezierCurveData functions

// TODO: make functions in extension private
extension BPBezierCurveData {
    private static var invalidLength: Double { -1.0 }
    
    // static FBBezierCurveData FBBezierCurveDataBezierClipWithBezierCurve(FBBezierCurveData me, FBBezierCurveData curve, FBBezierCurveData originalCurve, FBRange *originalRange, BOOL *intersects)
    // Becomes:
    // let (new_curve, success) = curve_data.bezierClipWithBezierCurve(curve, &originalCurve, &originalRange)
    static func bezierClipWithBezierCurve(_ me: BPBezierCurveData,
                                          curve: BPBezierCurveData,
                                          originalCurve: BPBezierCurveData,
                                          originalRange: inout ClosedRange<Double>) -> (clipped: BPBezierCurveData, intersects: Bool) {
        // This method does the clipping of self. It removes the parts of self that we can determine don't intersect
        //  with curve. It'll return the clipped version of self, update originalRange which corresponds to the range
        //  on the original curve that the return value represents. Finally, it'll set the intersects out parameter
        //  to yes or no depending on if the curves intersect or not.
        
        // Clipping works as follows:
        //  Draw a line through the two endpoints of the other curve, which we'll call the fat line. Measure the
        //  signed distance between the control points on the other curve and the fat line. The distance from the line
        //  will give us the fat line bounds. Any part of our curve that lies further away from the fat line than the
        //  fat line bounds we know can't intersect with the other curve, and thus can be removed.
        
        // We actually use two different fat lines. The first one uses the end points of the other curve, and the second
        //  one is perpendicular to the first. Most of the time, the first fat line will clip off more, but sometimes the
        //  second proves to be a better fat line in that it clips off more. We use both in order to converge more quickly.
        
        // Compute the regular fat line using the end points, then compute the range that could still possibly intersect
        //  with the other curve
        let (fatLine, fatLineBounds) = curve.regularFatLineBounds()
        let regularClippedRange = me.clipWithFatLine(fatLine, bounds: fatLineBounds)
        
        // A range of [1, 0] is a special sentinel value meaning "they don't intersect".
        // If they don't, bail early to save time
        if regularClippedRange.lowerBound == 1.0 && regularClippedRange.upperBound == 0.0 {
            return (clipped: me, intersects: false)
        }
        
        // Just in case the regular fat line isn't good enough, try the perpendicular one
        let (perpendicularLine, perpendicularLineBounds) = curve.perpendicularFatLineBounds()
        let perpendicularClippedRange = me.clipWithFatLine(perpendicularLine, bounds: perpendicularLineBounds)
        
        if perpendicularClippedRange.lowerBound == 1.0 && perpendicularClippedRange.upperBound == 0.0 {
            return (clipped: me, intersects: false)
        }
        
        // Combine to form Voltron.
        // Take the intersection of the regular fat line range and the perpendicular one.
        let clippedRange = max(regularClippedRange.lowerBound, perpendicularClippedRange.lowerBound)...min(regularClippedRange.upperBound, perpendicularClippedRange.upperBound)
        
        // Right now the clipped range is relative to ourself, not the original curve.
        // So map the newly clipped range onto the original range
        let newRange = RangeMath.scaleNormalizedValue(originalRange, value: clippedRange.lowerBound)...RangeMath.scaleNormalizedValue(originalRange, value: clippedRange.upperBound)
        
        originalRange = newRange
        
        // Actually divide the curve, but be sure to use the original curve. This helps with errors building up.
        return (clipped: originalCurve.subcurveWithRange(originalRange), intersects: true)
    }
    
    // static BOOL FBBezierCurveDataIsPoint(FBBezierCurveData *me)
    // Becomes:
    // if curve_data.isPoint()
    // SEE NEW ACCESSOR ABOVE: func isPoint() -> Bool
    
    
    // static NSRect FBBezierCurveDataBoundingRect(FBBezierCurveData *me)
    // Becomes:
    // let bnds = curve_data.boundingRect()
    // SEE NEW ACCESSOR ABOVE: func boundingRect() -> CGRect
    
    
    // static NSRect FBBezierCurveDataBounds(FBBezierCurveData* me)
    // Becomes:
    // let bnds = curve_data.bounds()
    // SEE NEW ACCESSOR ABOVE: func bounds() -> CGRect
    
    
    // 759
    //static void FBBezierCurveDataRefineIntersectionsOverIterations(NSUInteger iterations, FBRange *usRange, FBRange *themRange, FBBezierCurveData originalUs, FBBezierCurveData originalThem, FBBezierCurveData us, FBBezierCurveData them, FBBezierCurveData nonpointUs, FBBezierCurveData nonpointThem)
    // Becomes:
    // refineIntersectionsOverIterations(blah, blah, blah.....)
    // TODO: Need to decide whether this needs to be static
    
    static func refineIntersectionsOverIterations(_ iterations: Int,
                                                  usRange: inout ClosedRange<Double>,
                                                  themRange: inout ClosedRange<Double>,
                                                  originalUs: BPBezierCurveData,
                                                  originalThem: BPBezierCurveData,
                                                  us: inout BPBezierCurveData,
                                                  them: inout BPBezierCurveData,
                                                  nonpointUs: inout BPBezierCurveData,
                                                  nonpointThem: inout BPBezierCurveData) {
        for _ in 0..<iterations {
            var intersects = false
            
            (us, intersects) = bezierClipWithBezierCurve(us, curve: them, originalCurve: originalUs, originalRange: &usRange)
            if !intersects {
                (us, intersects) = bezierClipWithBezierCurve(nonpointUs, curve: nonpointThem, originalCurve: originalUs, originalRange: &usRange)
            }
            
            (them, intersects) = bezierClipWithBezierCurve(them, curve: us, originalCurve: originalThem, originalRange: &themRange)
            if !intersects {
                (them, intersects) = bezierClipWithBezierCurve(nonpointThem, curve: nonpointUs, originalCurve: originalThem, originalRange: &themRange)
            }
            
            if !them.isPoint {
                nonpointThem = them
            }
            
            if !us.isPoint {
                nonpointUs = us
            }
        }
    }
    
    // 777
    // static FBBezierCurveData FBBezierCurveDataClipLineOriginalCurve(FBBezierCurveData me, FBBezierCurveData originalCurve, FBBezierCurveData curve, FBRange *originalRange, FBBezierCurveData otherCurve, BOOL *intersects)
    // Becomes:
    // let (clippedCurve, intersects) = clipLineOriginalCurve(blah, blah, blah.....)
    
    /**
     **Objective-C:**
     
     FBBezierCurveData clippedCurve = FBBezierCurveDataClipLineOriginalCurve(FBBezierCurveData me, FBBezierCurveData originalCurve, FBBezierCurveData curve, FBRange *originalRange, FBBezierCurveData otherCurve, BOOL *intersects)
     
     **Becomes:**
     
     let (clippedCurve, intersects) = clipLineOriginalCurve(blah, blah, blah.....)
     
     - returns: **clippedCurve** - a curve clipped by another
     - returns: **intersects** - a success flag
     */
    static func clipLineOriginalCurve(_ originalCurve: BPBezierCurveData,
                                      curve: BPBezierCurveData,
                                      originalRange: inout ClosedRange<Double>,
                                      otherCurve: BPBezierCurveData) -> (clippedCurve: BPBezierCurveData, intersects: Bool) {
        let themOnUs1 = CurveHelpers.parameterOfPointOnLine(curve.endPoint1, lineEnd: curve.endPoint2, point: otherCurve.endPoint1)
        let themOnUs2 = CurveHelpers.parameterOfPointOnLine(curve.endPoint1, lineEnd: curve.endPoint2, point: otherCurve.endPoint2)
        let clippedRange = max(0.0, min(themOnUs1, themOnUs2))...min(1.0, max(themOnUs1, themOnUs2))
        
        if clippedRange.lowerBound > clippedRange.upperBound {
            // No intersection
            return (curve, false)
        }
        
        // Right now the clipped range is relative to ourself, not the original curve,
        // so map the newly clipped range onto the original range.
        originalRange = RangeMath.scaleNormalizedValue(
            originalRange,
            value: clippedRange.lowerBound
        )...RangeMath.scaleNormalizedValue(
            originalRange,
            value: clippedRange.upperBound
        )
        
        return (originalCurve.subcurveWithRange(originalRange) , true)
    }
    
    /*
     func test(c1: FBBezierCurveData, c2: FBBezierCurveData, inout r: FBRange, c3: FBBezierCurveData) {
     clipLineOriginalCurve(c1, c2, &r, c3)
     }
     */
    
    // 796
    // static BOOL FBBezierCurveDataCheckLinesForOverlap(FBBezierCurveData me, FBRange *usRange, FBRange *themRange, FBBezierCurveData originalUs, FBBezierCurveData originalThem, FBBezierCurveData *us, FBBezierCurveData *them)
    
    /**
     **Objective-C:**
     
     BOOL overlaps = FBBezierCurveDataCheckLinesForOverlap(FBBezierCurveData me, FBRange *usRange, FBRange *themRange, FBBezierCurveData originalUs, FBBezierCurveData originalThem, FBBezierCurveData *us, FBBezierCurveData *them)
     
     **Becomes:**
     
     let overlaps = checkLinesForOverlap(blah, blah, blah.....)
     
     */
    @discardableResult
    static func checkLinesForOverlap(_ me: BPBezierCurveData,
                                     usRange: inout ClosedRange<Double>,
                                     themRange: inout ClosedRange<Double>,
                                     originalUs: BPBezierCurveData,
                                     originalThem: BPBezierCurveData,
                                     us: inout BPBezierCurveData,
                                     them: inout BPBezierCurveData) -> Bool {
        // First see if its possible for them to overlap at all
        if !TangentMath.lineBoundsMightOverlap(us.bounds, bounds2: them.bounds) {
            return false
        }
        
        // Are all 4 points in a single line?
        let errorThreshold = 1e-7
        
        let isColinear = ProximityMath.areValuesClose(CurveHelpers.counterClockwiseTurn(us.endPoint1,
                                                                                        point2: us.endPoint2,
                                                                                        point3: them.endPoint1),
                                                      value2: 0.0,
                                                      threshold: errorThreshold)
        && ProximityMath.areValuesClose(CurveHelpers.counterClockwiseTurn(us.endPoint1,
                                                                          point2: us.endPoint2,
                                                                          point3: them.endPoint2),
                                        value2: 0.0,
                                        threshold: errorThreshold)
        
        if !isColinear {
            return false
        }
        
        var intersects = false
        (us, intersects) = clipLineOriginalCurve(originalUs, curve: us, originalRange: &usRange, otherCurve: them)
        
        if !intersects {
            return false
        }
        (them, intersects) = clipLineOriginalCurve(originalThem, curve: them, originalRange: &themRange, otherCurve: us)
        return intersects
    }
    
    // 906
    //static BOOL FBBezierCurveDataAreCurvesEqual(FBBezierCurveData me, FBBezierCurveData other)
    static func curvesAreEqual(_ me: BPBezierCurveData, other: BPBezierCurveData) -> Bool {
        guard !(me.isPoint || other.isPoint)
                || me.isStraightLine == other.isStraightLine else { return false }
        
        let endPointThreshold = 1e-4
        let controlPointThreshold = 1e-1
        if me.isStraightLine {
            return ProximityMath.arePointsClose(me.endPoint1, point2: other.endPoint1, threshold: endPointThreshold)
            && ProximityMath.arePointsClose(me.endPoint2, point2: other.endPoint2, threshold: endPointThreshold)
        } else {
            return ProximityMath.arePointsClose(me.endPoint1, point2: other.endPoint1, threshold: endPointThreshold)
            && ProximityMath.arePointsClose(me.controlPoint1, point2: other.controlPoint1, threshold: controlPointThreshold)
            && ProximityMath.arePointsClose(me.controlPoint2, point2: other.controlPoint2, threshold: controlPointThreshold)
            && ProximityMath.arePointsClose(me.endPoint2, point2: other.endPoint2, threshold: endPointThreshold)
        }
    }
    
    // 931
    //static FBBezierCurveData FBBezierCurveDataReversed(FBBezierCurveData me)
    static func reversed(_ me: BPBezierCurveData) -> BPBezierCurveData {
        return BPBezierCurveData(endPoint1: me.endPoint2, controlPoint1: me.controlPoint2, controlPoint2: me.controlPoint1, endPoint2: me.endPoint1, isStraightLine: me.isStraightLine)
    }
    
    // 936
    // static BOOL FBBezierCurveDataCheckForOverlapRange(FBBezierCurveData me, FBBezierIntersectRange **intersectRange, FBRange *usRange, FBRange *themRange, FBBezierCurve* originalUs, FBBezierCurve* originalThem, FBBezierCurveData us, FBBezierCurveData them)
    @discardableResult
    static func checkForOverlapRange(_ me: BPBezierCurveData,
                                     intersectRange: inout BPBezierIntersectRange?,
                                     usRange: ClosedRange<Double>,
                                     themRange: ClosedRange<Double>,
                                     originalUs: BPBezierCurve,
                                     originalThem: BPBezierCurve,
                                     us: BPBezierCurveData,
                                     them: BPBezierCurveData) -> Bool {
        if curvesAreEqual(us, other: them) {
            intersectRange = BPBezierIntersectRange(curve1: originalUs, parameterRange1: usRange, curve2: originalThem, parameterRange2: themRange, reversed: false)
            return true
        } else if curvesAreEqual(us, other: them.reversed()) {
            intersectRange = BPBezierIntersectRange(curve1: originalUs, parameterRange1: usRange, curve2: originalThem, parameterRange2: themRange, reversed: true)
            return true
        }
        
        return false
    }
    
    // 952
    //static FBBezierCurveData FBBezierCurveDataFindPossibleOverlap(FBBezierCurveData me, FBBezierCurveData originalUs, FBBezierCurveData them, FBRange *possibleRange)
    static func findPossibleOverlap(_ me: BPBezierCurveData,
                                    originalUs: BPBezierCurveData,
                                    them: BPBezierCurveData,
                                    possibleRange: inout ClosedRange<Double>) -> BPBezierCurveData {
        let themOnUs1 = originalUs.closestLocationToPoint(them.endPoint1)
        let themOnUs2 = originalUs.closestLocationToPoint(them.endPoint2)
        let range = min(themOnUs1.parameter, themOnUs2.parameter)...max(themOnUs1.parameter, themOnUs2.parameter)
        
        possibleRange = range;
        
        return originalUs.subcurveWithRange(range)
    }
    
    // 961
    //static BOOL FBBezierCurveDataCheckCurvesForOverlapRange(FBBezierCurveData me, FBBezierIntersectRange **intersectRange, FBRange *usRange, FBRange *themRange, FBBezierCurve* originalUs, FBBezierCurve* originalThem, FBBezierCurveData us, FBBezierCurveData them)
    static func checkCurvesForOverlapRange(_ me: BPBezierCurveData,
                                           intersectRange: inout BPBezierIntersectRange?,
                                           usRange: inout ClosedRange<Double>,
                                           themRange: inout ClosedRange<Double>,
                                           originalUs: BPBezierCurve,
                                           originalThem: BPBezierCurve,
                                           us: BPBezierCurveData,
                                           them: BPBezierCurveData) -> Bool {
        
        if checkForOverlapRange(me, intersectRange: &intersectRange, usRange: usRange, themRange: themRange, originalUs: originalUs, originalThem: originalThem, us: us, them: them) {
            return true
        }
        var usSubcurveRange = 0.0...0.0
        let usSubcurve = findPossibleOverlap(me, originalUs: originalUs.data, them: them, possibleRange: &usSubcurveRange)
        
        var themSubcurveRange = 0.0...0.0
        let themSubcurve = findPossibleOverlap(me, originalUs: originalThem.data, them: us, possibleRange: &themSubcurveRange)
        
        let threshold = 1e-4
        if usSubcurve.isEqual(to: themSubcurve, threshold: threshold)
            || usSubcurve.isEqual(to: reversed(themSubcurve), threshold: threshold) {
            usRange = usSubcurveRange
            themRange = themSubcurveRange
            
            return checkForOverlapRange(me, intersectRange: &intersectRange, usRange: usRange, themRange: themRange, originalUs: originalUs, originalThem: originalThem, us: usSubcurve, them: themSubcurve);
        }
        
        return false
    }
    
    // 982
    //static void FBBezierCurveDataCheckNoIntersectionsForOverlapRange(FBBezierCurveData me, FBBezierIntersectRange **intersectRange, FBRange *usRange, FBRange *themRange, FBBezierCurve* originalUs, FBBezierCurve* originalThem, FBBezierCurveData us, FBBezierCurveData them, FBBezierCurveData nonpointUs, FBBezierCurveData nonpointThem)
    static func checkNoIntersectionsForOverlapRange(_ me: BPBezierCurveData,
                                                    intersectRange: inout BPBezierIntersectRange?,
                                                    usRange: inout ClosedRange<Double>,
                                                    themRange: inout ClosedRange<Double>,
                                                    originalUs: BPBezierCurve,
                                                    originalThem: BPBezierCurve,
                                                    us: inout BPBezierCurveData,
                                                    them: inout BPBezierCurveData,
                                                    nonpointUs: BPBezierCurveData,
                                                    nonpointThem: BPBezierCurveData) {
        
        if us.isStraightLine && them.isStraightLine {
            checkLinesForOverlap(me, usRange: &usRange, themRange: &themRange, originalUs: originalUs.data, originalThem: originalThem.data, us: &us, them: &them)
        }
        
        checkForOverlapRange(me, intersectRange: &intersectRange, usRange: usRange, themRange: themRange, originalUs: originalUs, originalThem: originalThem, us: us, them: them)
    }
    
    // 990
    //static BOOL FBBezierCurveDataCheckForStraightLineOverlap(FBBezierCurveData me, FBBezierIntersectRange **intersectRange, FBRange *usRange, FBRange *themRange, FBBezierCurve* originalUs, FBBezierCurve* originalThem, FBBezierCurveData us, FBBezierCurveData them, FBBezierCurveData nonpointUs, FBBezierCurveData nonpointThem)
    static func straightLineOverlap(_ me: BPBezierCurveData,
                                    intersectRange: inout BPBezierIntersectRange?,
                                    usRange: inout ClosedRange<Double>,
                                    themRange: inout ClosedRange<Double>,
                                    originalUs: BPBezierCurve,
                                    originalThem: BPBezierCurve,
                                    us: inout BPBezierCurveData,
                                    them: inout BPBezierCurveData,
                                    nonpointUs: BPBezierCurveData,
                                    nonpointThem: BPBezierCurveData) -> Bool {
        var hasOverlap = false
        
        if us.isStraightLine && them.isStraightLine {
            hasOverlap = checkLinesForOverlap(me, usRange: &usRange, themRange: &themRange, originalUs: originalUs.data, originalThem: originalThem.data, us: &us, them: &them)
        }
        
        if hasOverlap {
            hasOverlap = checkForOverlapRange(me, intersectRange: &intersectRange, usRange: usRange, themRange: themRange, originalUs: originalUs, originalThem: originalThem, us: us, them: them)
        }
        
        return hasOverlap
    }
    
    static func pfRefineParameter(_ me: BPBezierCurveData,
                                  parameter: CGFloat,
                                  point: CGPoint) -> CGFloat {
        // Use Newton's Method to refine our parameter. In general, that formula is:
        //
        //  parameter = parameter - f(parameter) / f'(parameter)
        //
        // In our case:
        //
        //  f(parameter) = (Q(parameter) - point) * Q'(parameter) = 0
        //
        // Where Q'(parameter) is tangent to the curve at Q(parameter) and orthogonal to [Q(parameter) - P]
        //
        // Taking the derivative gives us:
        //
        //  f'(parameter) = (Q(parameter) - point) * Q''(parameter) + Q'(parameter) * Q'(parameter)
        //
        
        let bezierPoints: [CGPoint] = [me.endPoint1, me.controlPoint1, me.controlPoint2, me.endPoint2]
        
        // Compute Q(parameter)
        let qAtParameter = CurveHelpers.bezierWithPoints(3, bezierPoints: bezierPoints, parameter: parameter, withCurves: false).point
        
        // Compute Q'(parameter)
        let qPrimePoints: [CGPoint] = [
            CGPoint(x:(bezierPoints[1].x - bezierPoints[0].x) * 3.0,
                    y:(bezierPoints[1].y - bezierPoints[0].y) * 3.0),
            CGPoint(x:(bezierPoints[2].x - bezierPoints[1].x) * 3.0,
                    y:(bezierPoints[2].y - bezierPoints[1].y) * 3.0),
            CGPoint(x:(bezierPoints[3].x - bezierPoints[2].x) * 3.0,
                    y:(bezierPoints[3].y - bezierPoints[2].y) * 3.0)
        ]
        let qPrimeAtParameter = CurveHelpers.bezierWithPoints(2, bezierPoints: qPrimePoints, parameter: parameter, withCurves: false).point
        
        // Compute Q''(parameter)
        let qPrimePrimePoints: [CGPoint] = [
            CGPoint(
                x: (qPrimePoints[1].x - qPrimePoints[0].x) * 2.0,
                y: (qPrimePoints[1].y - qPrimePoints[0].y) * 2.0
            ),
            CGPoint(
                x: (qPrimePoints[2].x - qPrimePoints[1].x) * 2.0,
                y: (qPrimePoints[2].y - qPrimePoints[1].y) * 2.0
            )
        ]
        let qPrimePrimeAtParameter = CurveHelpers.bezierWithPoints(1, bezierPoints: qPrimePrimePoints, parameter: parameter, withCurves: false).point
        
        // Compute f(parameter) and f'(parameter)
        let qMinusPoint = PointMath.subtractPoint(qAtParameter, point2: point)
        let fAtParameter = PointMath.dotMultiplyPoint(qMinusPoint, point2: qPrimeAtParameter)
        let fPrimeAtParameter = PointMath.dotMultiplyPoint(qMinusPoint, point2: qPrimePrimeAtParameter)
            + PointMath.dotMultiplyPoint(qPrimeAtParameter, point2: qPrimeAtParameter)
        
        // Newton's method!
        return parameter - (fAtParameter / fPrimeAtParameter)
    }
    
    // 1050
    //static FBBezierIntersectRange *FBBezierCurveDataMergeIntersectRange(FBBezierIntersectRange *intersectRange, FBBezierIntersectRange *otherIntersectRange)
    static func mergeIntersectRange(_ intersectRange: BPBezierIntersectRange?,
                                    otherIntersectRange: BPBezierIntersectRange?) -> BPBezierIntersectRange? {
        guard let otherIntersectRange else { return intersectRange }
        guard let intersectRange else { return otherIntersectRange }
        
        intersectRange.merge(otherIntersectRange)
        
        return intersectRange
    }
    
    @discardableResult
    static func intersectionsWithStraightLines(_ me: BPBezierCurveData,
                                               curve: BPBezierCurveData,
                                               usRange: inout ClosedRange<Double>,
                                               themRange: inout ClosedRange<Double>,
                                               originalUs: BPBezierCurve,
                                               originalThem: BPBezierCurve,
                                               stop: inout Bool,
                                               outputBlock: (_ intersect: BPBezierIntersection) -> (setStop: Bool, stopValue:Bool)) -> Bool {
        if !me.isStraightLine || !curve.isStraightLine {
            return false
        }
        
        var intersectionPoint = CGPoint.zero
        let intersects = CurveHelpers.linesIntersect(me.endPoint1, line1End: me.endPoint2, line2Start: curve.endPoint1, line2End: curve.endPoint2, outIntersect: &intersectionPoint)
        
        if !intersects {
            return false
        }
        
        let meParam = CurveHelpers.parameterOfPointOnLine(me.endPoint1, lineEnd: me.endPoint2, point: intersectionPoint)
        if ComparisonMath.isValueLessThan(meParam, maximum: 0.0) || ComparisonMath.isValueGreaterThan(meParam, minimum: 1.0) {
            return false
        }
        
        let curveParam = CurveHelpers.parameterOfPointOnLine(curve.endPoint1, lineEnd: curve.endPoint2, point: intersectionPoint)
        if ComparisonMath.isValueLessThan(curveParam, maximum: 0.0) || ComparisonMath.isValueGreaterThan(curveParam, minimum: 1.0) {
            return false
        }
        
        let intersect = BPBezierIntersection(curve1: originalUs, param1: meParam, curve2: originalThem, param2: curveParam)
        
        let stopResults = outputBlock(intersect)
        if stopResults.setStop {
            stop = stopResults.stopValue
        }
        
        return true
    }
}

// ========================================================
// MARK: ---- MAIN SPLIT FUNCTION ----
// ========================================================

internal func pfIntersectionsWithBezierCurve(_ me: BPBezierCurveData,
                                             curve: BPBezierCurveData,
                                             usRange: inout ClosedRange<Double>,
                                             themRange: inout ClosedRange<Double>,
                                             originalUs: BPBezierCurve,
                                             originalThem: BPBezierCurve,
                                             intersectRange: inout BPBezierIntersectRange?,
                                             depth: Int,
                                             stop: inout Bool,
                                             outputBlock: (_ intersect: BPBezierIntersection) -> (setStop: Bool, stopValue:Bool)) {
    // This is the main work loop.
    // At a high level this method sits in a loop and removes sections (ranges)
    // of the two bezier curves that it knows don't intersect.
    // (how it knows that is covered in the appropriate method).
    // The idea is to whittle the curves down to the point where they do intersect.
    //
    // When the range where they intersect converges (i.e. matches to 6 decimal places)
    // or there are more than 500 attempts, the loop stops.
    //
    // A special case is when we're not able to remove at least 20% of
    // the curves on a given interation.
    // In that case we assume there are likely multiple intersections, so we
    // divide one of curves in half, and recurse on the two halves.
    
    let places = 6 // How many decimals place to calculate the solution out to
    let maxIterations = 500 // how many iterations to allow before we just give up
    let maxDepth = 10 // how many recursive calls to allow before we just give up
    let minimumChangeNeeded = 0.20 // how much to clip off for a given iteration minimum before we subdivide the curve
    
    // us is self, but will become clipped down to where the intersection is (perhaps)
    var us = BPBezierCurveData(cloning: me)
    
    // them is the other curve we're intersecting with, but clipped down to where the intersection is
    var them = BPBezierCurveData(cloning: curve)
    
    var nonpointUs = BPBezierCurveData(cloning: us)
    var nonpointThem = BPBezierCurveData(cloning: them)
    
    // Horizontal and vertical lines are somewhat special cases,
    // and the math doesn't always work out that great.
    // For example, two vertical lines that overlap will kick out
    // as intersecting at the endpoints.
    //
    // Try to detect that kind of overlap at the start.
    if BPBezierCurveData.straightLineOverlap(me, intersectRange: &intersectRange, usRange: &usRange, themRange: &themRange, originalUs: originalUs, originalThem: originalThem, us: &us, them: &them, nonpointUs: nonpointUs, nonpointThem: nonpointThem) {
        return
    }
    
    if us.isStraightLine && them.isStraightLine {
        BPBezierCurveData.intersectionsWithStraightLines(me, curve: curve, usRange: &usRange, themRange: &themRange, originalUs: originalUs, originalThem: originalThem, stop: &stop, outputBlock: outputBlock)
        return
    }
    
    let originalUsData = BPBezierCurveData(cloning: originalUs.data)
    let originalThemData = BPBezierCurveData(cloning: originalThem.data)
    
    // Don't check for convergence until we actually see if we intersect or not.
    // i.e. Make sure we go through at least once, otherwise the results don't mean anything.
    //
    // Be sure to stop as soon as either range converges, otherwise calculations
    // for the other range goes funky because one curve is essentially a point.
    var iterations = 0
    var hadConverged = true
    
    while (iterations < maxIterations
           && (iterations == 0 || (!RangeMath.hasConverged(usRange, decimalPlaces: places) || !RangeMath.hasConverged(themRange, decimalPlaces: places))) ) {
        // Remember what the current range is so we can calculate how much it changed later
        let previousUsRange = usRange
        let previousThemRange = themRange
        
        // Remove the range from ourselves that doesn't intersect with them.
        // If the other curve is already a point, use the previous iteration's
        //  copy of them so calculations still work.
        var intersects = false
        if !them.isPoint {
            nonpointThem = them
        }
        (us, intersects) = BPBezierCurveData.bezierClipWithBezierCurve(nonpointUs, curve: nonpointThem, originalCurve: originalUsData, originalRange: &usRange);
        if !intersects {
            BPBezierCurveData.checkNoIntersectionsForOverlapRange(me, intersectRange: &intersectRange, usRange: &usRange, themRange: &themRange, originalUs: originalUs, originalThem: originalThem, us: &us, them: &them, nonpointUs: nonpointUs, nonpointThem: nonpointThem)
            return  // If they don't intersect at all stop now
        }
        if iterations > 0 && (us.isPoint || them.isPoint) {
            break
        }
        
        // Remove the range of them that doesn't intersect with us
        if !us.isPoint {
            nonpointUs = us
        } else if iterations == 0 {
            // If the first time through, "us" was reduced to a point, then we're never
            // going to know if the curves actually intersect, even if both ranges converge.
            //
            // The ranges can converge on the parameters on each respective curve that is
            // closest to the other. But without being clipped to a smaller range,
            // the algorithm won't necessarily detect that they don't actually intersect
            hadConverged = false
        }
        
        (them, intersects) = BPBezierCurveData.bezierClipWithBezierCurve(nonpointThem, curve: nonpointUs, originalCurve: originalThemData, originalRange: &themRange)
        if !intersects {
            BPBezierCurveData.checkNoIntersectionsForOverlapRange(me, intersectRange: &intersectRange, usRange: &usRange, themRange: &themRange, originalUs: originalUs, originalThem: originalThem, us: &us, them: &them, nonpointUs: nonpointUs, nonpointThem: nonpointThem)
            return  // If they don't intersect at all stop now
        }
        if iterations > 0 && (us.isPoint || them.isPoint) {
            break
        }
        
        // See if either of curves ranges is reduced by less than 20%.
        let percentChangeInUs = (RangeMath.getSize(previousUsRange) - RangeMath.getSize(usRange)) / RangeMath.getSize(previousUsRange)
        let percentChangeInThem = (RangeMath.getSize(previousThemRange) - RangeMath.getSize(themRange)) / RangeMath.getSize(previousThemRange)
        
        var didNotSplit = false
        
        if percentChangeInUs < minimumChangeNeeded && percentChangeInThem < minimumChangeNeeded {
            // We're not converging fast enough, likely because there are
            // multiple intersections here.
            
            // Or the curves are the same, check for that first.
            if BPBezierCurveData.checkCurvesForOverlapRange(me, intersectRange: &intersectRange, usRange: &usRange, themRange: &themRange, originalUs: originalUs, originalThem: originalThem, us: us, them: them) {
                return
            }
            
            // Divide and conquer. Divide the longer curve in half, and recurse
            if RangeMath.getSize(usRange) > RangeMath.getSize(themRange) {
                // Since our remaining range is longer, split the remains of us in half at the midway point
                var usRange1 = usRange.lowerBound...((usRange.lowerBound + usRange.upperBound) / 2.0)
                let us1 = originalUsData.subcurveWithRange(usRange1)
                // make a local copy because it'll get modified when we recurse
                var themRangeCopy1 = themRange
                
                var usRange2 = ((usRange.lowerBound + usRange.upperBound) / 2.0)...usRange.upperBound
                let us2 = originalUsData.subcurveWithRange(usRange2)
                // make a local copy because it'll get modified when we recurse
                var themRangeCopy2 = themRange
                
                let range1ConvergedAlready = RangeMath.hasConverged(usRange1, decimalPlaces: places) && RangeMath.hasConverged(themRange, decimalPlaces: places)
                let range2ConvergedAlready = RangeMath.hasConverged(usRange2, decimalPlaces: places) && RangeMath.hasConverged(themRange, decimalPlaces: places);
                
                if !range1ConvergedAlready && !range2ConvergedAlready && depth < maxDepth {
                    // Compute the intersections between the two halves of us and them
                    var leftIntersectRange: BPBezierIntersectRange? // = nil
                    pfIntersectionsWithBezierCurve(us1, curve: them, usRange: &usRange1, themRange: &themRangeCopy1, originalUs: originalUs, originalThem: originalThem, intersectRange: &leftIntersectRange, depth: depth + 1, stop: &stop, outputBlock: outputBlock)
                    
                    intersectRange = BPBezierCurveData.mergeIntersectRange(intersectRange, otherIntersectRange: leftIntersectRange)
                    
                    if stop {
                        return
                    }
                    
                    var rightIntersectRange: BPBezierIntersectRange? // = nil
                    pfIntersectionsWithBezierCurve(us2, curve: them, usRange: &usRange2, themRange: &themRangeCopy2, originalUs: originalUs, originalThem: originalThem, intersectRange: &rightIntersectRange, depth: depth + 1, stop: &stop, outputBlock: outputBlock);
                    
                    intersectRange = BPBezierCurveData.mergeIntersectRange(intersectRange, otherIntersectRange: rightIntersectRange)
                    return
                } else {
                    didNotSplit = true
                }
            } else {
                // Since their remaining range is longer, split the
                // remains of them in half at the midway point
                var themRange1 = themRange.lowerBound...((themRange.lowerBound + themRange.upperBound) / 2.0)
                let them1 = originalThemData.subcurveWithRange(themRange1)
                // make a local copy because it'll get modified when we recurse
                var usRangeCopy1 = usRange 
                
                var themRange2 = ((themRange.lowerBound + themRange.upperBound) / 2.0)...themRange.upperBound
                let them2 = originalThemData.subcurveWithRange(themRange2)
                // make a local copy because it'll get modified when we recurse
                var usRangeCopy2 = usRange 
                
                let range1ConvergedAlready = RangeMath.hasConverged(themRange1, decimalPlaces: places) && RangeMath.hasConverged(usRange, decimalPlaces: places)
                let range2ConvergedAlready = RangeMath.hasConverged(themRange2, decimalPlaces: places) && RangeMath.hasConverged(usRange, decimalPlaces: places)
                
                if !range1ConvergedAlready && !range2ConvergedAlready && depth < maxDepth {
                    // Compute the intersections between the two halves of them and us
                    var leftIntersectRange: BPBezierIntersectRange? // = nil
                    pfIntersectionsWithBezierCurve(us, curve: them1, usRange: &usRangeCopy1, themRange: &themRange1, originalUs: originalUs, originalThem: originalThem, intersectRange: &leftIntersectRange, depth: depth + 1, stop: &stop, outputBlock: outputBlock)
                    
                    intersectRange = BPBezierCurveData.mergeIntersectRange(intersectRange, otherIntersectRange: leftIntersectRange)
                    
                    if stop {
                        return
                    }
                    
                    var rightIntersectRange: BPBezierIntersectRange? // = nil
                    pfIntersectionsWithBezierCurve(us, curve: them2, usRange: &usRangeCopy2, themRange: &themRange2, originalUs: originalUs, originalThem: originalThem, intersectRange: &rightIntersectRange, depth: depth + 1, stop: &stop, outputBlock: outputBlock)
                    
                    intersectRange = BPBezierCurveData.mergeIntersectRange(intersectRange, otherIntersectRange: rightIntersectRange)
                    
                    return
                } else {
                    didNotSplit = true
                }
            }
            
            if didNotSplit && (RangeMath.getSize(previousUsRange) - RangeMath.getSize(usRange)) == 0 && (RangeMath.getSize(previousThemRange) - RangeMath.getSize(themRange)) == 0 {
                // We're not converging at _all_ and we can't split, so we need to bail out.
                return // no intersections
            }
        }
        
        iterations += 1
    }
    
    // It's possible that one of the curves has converged, but the other hasn't.
    // Since the math becomes wonky once a curve becomes a point, the loop stops
    // as soon as either curve converges.
    //
    // However for our purposes we need _both_ curves to converge; that is to say
    // we need the parameter for each curve where they intersect.
    // Fortunately, since one curve did converge we know the 2D point where they
    // converge, plus we have a reasonable approximation for the parameter for
    // the curve that didn't.
    // That means we can use Newton's method to refine the parameter of the curve
    // that didn't converge.
    if !RangeMath.hasConverged(usRange, decimalPlaces: places) || !RangeMath.hasConverged(themRange, decimalPlaces: places) {
        if BPBezierCurveData.checkCurvesForOverlapRange(me,
                                                        intersectRange: &intersectRange,
                                                        usRange: &usRange,
                                                        themRange: &themRange,
                                                        originalUs: originalUs,
                                                        originalThem: originalThem,
                                                        us: originalUsData,
                                                        them: originalThemData) {
            return
        }
        
        // We bail out of the main loop as soon as we know things intersect, but before
        // the math falls apart. Unfortunately sometimes this means we don't always get
        // the best estimate of the parameters.
        //
        // Below we fall back to Netwon's method, but it's accuracy is dependant on our
        // previous calculations.
        //
        // So here assume things intersect and just try to tighten up the parameters.
        // If the math falls apart because everything's a point, that's OK since we
        // already have a "reasonable" estimation of the parameters.
        
        BPBezierCurveData.refineIntersectionsOverIterations(3,
                                                            usRange: &usRange,
                                                            themRange: &themRange,
                                                            originalUs: originalUsData,
                                                            originalThem: originalThemData,
                                                            us: &us,
                                                            them: &them,
                                                            nonpointUs: &nonpointUs,
                                                            nonpointThem: &nonpointThem)
        
        // Sometimes we need a little more precision. Be careful though, because
        // in some cases trying for more makes the math fall apart.
        
        if !RangeMath.hasConverged(usRange, decimalPlaces: places) || !RangeMath.hasConverged(themRange, decimalPlaces: places) {
            BPBezierCurveData.refineIntersectionsOverIterations(4,
                                                                usRange: &usRange,
                                                                themRange: &themRange,
                                                                originalUs: originalUsData,
                                                                originalThem: originalThemData,
                                                                us: &us,
                                                                them: &them,
                                                                nonpointUs: &nonpointUs,
                                                                nonpointThem: &nonpointThem)
        }
    }
    
    if RangeMath.hasConverged(usRange, decimalPlaces: places) && !RangeMath.hasConverged(themRange, decimalPlaces: places) {
        // Refine the them range since it didn't converge
        let intersectionPoint = originalUsData.pointAtParameter(RangeMath.average(usRange)).point
        
        // Although the range didn't converge, it should be a reasonable
        // approximation which is all Newton needs
        var refinedParameter = RangeMath.average(themRange)
        for _ in 0..<3 {
            refinedParameter = BPBezierCurveData.pfRefineParameter(originalThemData, parameter: refinedParameter, point: intersectionPoint)
            refinedParameter = min(themRange.upperBound, max(themRange.lowerBound, refinedParameter))
        }
        themRange = refinedParameter...refinedParameter
//        themRange.lowerBound = refinedParameter
//        themRange.upperBound = refinedParameter
        hadConverged = false
    } else if !RangeMath.hasConverged(usRange, decimalPlaces: places) && RangeMath.hasConverged(themRange, decimalPlaces: places) {
        // Refine the us range since it didn't converge
        let intersectionPoint = originalThemData.pointAtParameter(RangeMath.average(themRange)).point
        
        // Although the range didn't converge, it should be a reasonable
        // approximation which is all Newton needs
        var refinedParameter = RangeMath.average(usRange)
        for _ in 0 ..< 3 {
            refinedParameter = BPBezierCurveData.pfRefineParameter(originalUsData, parameter: refinedParameter, point: intersectionPoint)
            refinedParameter = min(usRange.upperBound, max(usRange.lowerBound, refinedParameter))
        }
        themRange = refinedParameter...refinedParameter
//        usRange.lowerBound = refinedParameter
//        usRange.upperBound = refinedParameter
        hadConverged = false
    }
    
    // If it never converged and we stopped because of our loop max,
    // assume overlap or something else. Bail.
    if (!RangeMath.hasConverged(usRange, decimalPlaces: places) || !RangeMath.hasConverged(themRange, decimalPlaces: places)) && iterations >= maxIterations {
        BPBezierCurveData.checkForOverlapRange(me, intersectRange: &intersectRange, usRange: usRange, themRange: themRange, originalUs: originalUs, originalThem: originalThem, us: us, them: them)
        return
    }
    
    if !hadConverged {
        // Since one of them didn't converge, we need to make sure they actually intersect.
        // Compute the point from both and compare
        let intersectionPoint = originalUsData.pointAtParameter(RangeMath.average(usRange)).point
        let checkPoint = originalThemData.pointAtParameter(RangeMath.average(themRange)).point
        let threshold = 1e-3
        
        if !ProximityMath.arePointsClose(intersectionPoint, point2: checkPoint, threshold: threshold) {
            return
        }
    }
    
    // Return the final intersection, which we represent by the original curves
    // and the parameters where they intersect.
    //
    // The parameter values are useful later in the boolean operations,
    // plus it allows us to do lazy calculations.
    let intersection = BPBezierIntersection(curve1: originalUs, param1: RangeMath.average(usRange), curve2:originalThem, param2: RangeMath.average(themRange))
    
    let stopResults = outputBlock(intersection)
    if stopResults.setStop {
        stop = stopResults.stopValue
    }
}

// MARK: - BPBezierCurve

/// The main purpose of this class is to compute the intersections of two bezier
/// curves. It does this using the bezier clipping algorithm, described in
/// "Curve intersection using Bezier clipping" by TW Sederberg and T Nishita.
/// https://cagd.cs.byu.edu/~tom/papers/bezclip.pdf
public class BPBezierCurve {
    // MARK: edge extensions
    fileprivate var _startShared = false
    fileprivate var _contour: BPBezierContour?
    fileprivate var _index: Int = 0
    var crossings: [BPEdgeCrossing] = []
    
    // 89 of Edge extension
    var index: Int {
        get {
            return _index
        }
        set {
            _index = newValue
        }
    }
    
    // 99 of Edge extension
    var isStartShared: Bool {
        return _startShared
    }
    
    var startShared: Bool {
        get {
            return _startShared
        }
        set {
            _startShared = newValue
        }
    }
    
    // 109 of Edge extension
    var contour: BPBezierContour? {
        get {
            return _contour
        }
        set {
            _contour = newValue
        }
    }
    
    private(set) var data: BPBezierCurveData
    
    init(endPoint1: CGPoint, controlPoint1: CGPoint, controlPoint2: CGPoint, endPoint2: CGPoint) {
        data = BPBezierCurveData(
            endPoint1: endPoint1,
            controlPoint1: controlPoint1,
            controlPoint2: controlPoint2,
            endPoint2: endPoint2,
            isStraightLine: false)
    }
    
    init(startPoint: CGPoint, endPoint: CGPoint) {
        // Convert the line into a bezier curve to keep our intersection algorithm general (i.e. only
        //  has to deal with curves, not lines). As long as the control points are colinear with the
        //  end points, it'll be a line. But for consistency sake, we put the control points inside
        //  the end points, 1/3 of the total distance away from their respective end point.
        let distance = PointMath.distanceBetweenPoints(startPoint, point2: endPoint)
        let leftTangent = PointMath.normalizePoint(PointMath.subtractPoint(endPoint, point2: startPoint))
        
        data = BPBezierCurveData(
            endPoint1: startPoint,
            controlPoint1: PointMath.addPoint(startPoint, point2: PointMath.unitScalePoint(leftTangent, scale: distance / 3.0)),
            controlPoint2: PointMath.addPoint(startPoint, point2: PointMath.unitScalePoint(leftTangent, scale: 2.0 * distance / 3.0)),
            endPoint2: endPoint,
            isStraightLine: true)
    }
    
    init(curveData: BPBezierCurveData) {
        data = curveData;
    }
    
    var endPoint1: CGPoint {
        get {
            return data.endPoint1
        }
    }
    
    var endPoint2: CGPoint {
        get {
            return data.endPoint2
        }
    }
    
    var controlPoint1: CGPoint {
        get {
            return data.controlPoint1
        }
    }
    
    var controlPoint2: CGPoint {
        get {
            return data.controlPoint2
        }
    }
    
    var isStraightLine: Bool {
        get {
            return data.isStraightLine
        }
    }
    
    // MARK: - bezierCurvesFromPathElements
    // 1336
    
    public class func bezierCurvesFromPathElements(_ elements: [Path.Element]) -> [BPBezierCurve] {
        // Helper method to easily convert a bezier path into an array of BPBezierCurves.
        // Very straight-forward, only lines are a special case.
        var startPoint: CGPoint?
        
//        let bezier = LRTPathWrapper(path)
        var bezierCurves: [BPBezierCurve] = []
        
        var previousPoint = CGPoint.zero
        
        for item in elements {
            switch item {
            case let .move(v):
                previousPoint = v
                startPoint = v
            case let .line(v):
                // Convert lines to bezier curves as well.
                // Just set control point to be in the line formed by the end points
                bezierCurves.append(BPBezierCurve(startPoint: previousPoint, endPoint: v))
                previousPoint = v
            case .quadCurve(let to, let control):
                let twoThirds: CGFloat = 2.0 / 3.0
                // lastPoint + twoThirds * (via - lastPoint)
                let control1 = PointMath.addPoint(previousPoint, point2: PointMath.scalePoint(PointMath.subtractPoint(control, point2: previousPoint), scale: twoThirds))
                // toPt + twoThirds * (via - toPt)
                let control2 = PointMath.addPoint(to, point2: PointMath.scalePoint(PointMath.subtractPoint(control, point2: to), scale: twoThirds))
                bezierCurves.append(BPBezierCurve(endPoint1: previousPoint, controlPoint1: control1, controlPoint2: control2, endPoint2: to))
                previousPoint = to
            case .curve(let to, let control1, let control2):
                bezierCurves.append(BPBezierCurve(endPoint1: previousPoint, controlPoint1: control1, controlPoint2: control2, endPoint2: to))
                previousPoint = to
            case .closeSubpath:
                // Create a line back to the start if required
                if let startPoint = startPoint {
                    if !previousPoint.equalTo(startPoint) {
                        bezierCurves.append(BPBezierCurve(startPoint: previousPoint, endPoint: startPoint))
                    }
                }
                startPoint = nil
                previousPoint = CGPoint.zero
            }
        }
        
        return bezierCurves
    }
    
    // 1376
    class func bezierCurveWithLineStartPoint(_ startPoint: CGPoint, endPoint: CGPoint) -> BPBezierCurve {
        return BPBezierCurve(startPoint: startPoint, endPoint: endPoint)
    }
    
    // 1381
    class func bezierCurveWithEndPoint1(_ endPoint1: CGPoint, controlPoint1: CGPoint, controlPoint2: CGPoint, endPoint2: CGPoint) -> BPBezierCurve {
        return BPBezierCurve(endPoint1: endPoint1, controlPoint1: controlPoint1, controlPoint2: controlPoint2, endPoint2: endPoint2)
    }
    
    // 1386
    class func bezierCurveWithBezierCurveData(_ data: BPBezierCurveData) -> BPBezierCurve {
        return BPBezierCurve(curveData: data)
    }
    
//    static func dataIsEqual(_ me: BPBezierCurveData, other: BPBezierCurveData) -> Bool {
//        return me.isEqualWithOptions(other, threshold: 1e-10)
//    }
    
    func isEqual(to other: BPBezierCurve) -> Bool {
        data.isEqual(to: other.data, threshold: 1e-10)
    }
    
    func doesHaveIntersectionsWithBezierCurve(_ curve: BPBezierCurve) -> Bool {
        // (intersect: FBBezierIntersection) -> (setStop: Bool, stopValue:Bool)
        var count = 0
        var unusedRange: BPBezierIntersectRange?
        
        intersectionsWithBezierCurve(curve, overlapRange: &unusedRange) {
            (intersect: BPBezierIntersection) -> (setStop: Bool, stopValue:Bool) in
            count += 1
            return (setStop:true, stopValue:true) // Only need the one
        }
        return count > 0;
    }
    
    // 1457
    //- (void) intersectionsWithBezierCurve:(FBBezierCurve *)curve overlapRange:(FBBezierIntersectRange **)intersectRange withBlock:(FBCurveIntersectionBlock)block
    public func intersectionsWithBezierCurve(_ curve: BPBezierCurve, overlapRange: inout BPBezierIntersectRange?,
                                             withBlock block: (_ intersect: BPBezierIntersection) -> (setStop: Bool, stopValue:Bool)) {
        // For performance reasons, do a quick bounds check to see if these even might intersect
        guard TangentMath.lineBoundsMightOverlap(data.boundingRect, bounds2: curve.data.boundingRect),
              TangentMath.lineBoundsMightOverlap(data.bounds,       bounds2: curve.data.bounds) else { return }
        
        var usRange: ClosedRange<Double> = 0...1
        var themRange: ClosedRange<Double> = 0...1
        var stop = false
        pfIntersectionsWithBezierCurve(data, curve: curve.data, usRange: &usRange, themRange: &themRange, originalUs: self, originalThem: curve, intersectRange: &overlapRange, depth: 0, stop: &stop, outputBlock: block)
    }
    
    // 1473
    //- (FBBezierCurve *) subcurveWithRange:(FBRange)range
    func subcurveWithRange(_ range: ClosedRange<Double>) -> BPBezierCurve {
        return BPBezierCurve(curveData: data.subcurveWithRange(range))
    }
    
    // 1478
    //- (void) splitSubcurvesWithRange:(FBRange)range left:(FBBezierCurve **)leftCurve middle:(FBBezierCurve **)middleCurve right:(FBBezierCurve **)rightCurve
    func splitSubcurvesWithRange(_ range: ClosedRange<Double>, left: Bool, middle: Bool, right: Bool) -> (left: BPBezierCurve?, mid: BPBezierCurve?, right: BPBezierCurve?) {
        // Return a bezier curve representing the parameter range specified.
        // We do this by splitting twice:
        //   once on the minimum, then splitting the result of that on the maximum.
        var leftResult: BPBezierCurve?
        var midResult: BPBezierCurve?
        var rightResult: BPBezierCurve?
        
        // Start with the left side curve
        var remainingCurve: BPBezierCurveData?
        if range.lowerBound == 0.0 {
            remainingCurve = data
        } else {
            var leftCurveData: BPBezierCurveData?
            let pap = data.pointAtParameter(range.lowerBound)
            leftCurveData = pap.leftCurve
            remainingCurve = pap.rightCurve
            if left {
                if let leftCurveData = leftCurveData {
                    leftResult = BPBezierCurve(curveData: leftCurveData)
                }
            }
        }
        // Special case where we start at the end
        if range.lowerBound == 1.0 {
            if middle {
                if let remainingCurve = remainingCurve {
                    midResult = BPBezierCurve(curveData: remainingCurve)
                }
            }
        } else if let remainingCurve = remainingCurve {
            // We need to adjust the maximum parameter to
            // fit on the new curve before we split again
            let adjustedMaximum = (range.upperBound - range.lowerBound) / (1.0 - range.lowerBound)
            let pap = remainingCurve.pointAtParameter(adjustedMaximum)
            if middle {
                if let curveData = pap.leftCurve {
                    midResult = BPBezierCurve(curveData: curveData)
                }
            }
            if right {
                if let curveData = pap.rightCurve {
                    rightResult = BPBezierCurve(curveData: curveData)
                }
            }
        }
        return (left: leftResult, mid: midResult, right: rightResult)
    }
    
    // 1516
    //- (FBBezierCurve *) reversedCurve
    func reversedCurve() -> BPBezierCurve {
        return BPBezierCurve(curveData: BPBezierCurveData.reversed(data))
    }
    
    // 1521
    //- (NSPoint) pointAtParameter:(CGFloat)parameter leftBezierCurve:(FBBezierCurve **)leftBezierCurve rightBezierCurve:(FBBezierCurve **)rightBezierCurve
    func pointAtParameter(_ parameter: Double) -> (point: CGPoint, leftBezierCurve: BPBezierCurve?, rightBezierCurve: BPBezierCurve?) {
        var leftBezierCurve: BPBezierCurve?
        var rightBezierCurve: BPBezierCurve?
        
        let pap = data.pointAtParameter(parameter)
        if let leftData = pap.leftCurve {
            leftBezierCurve = BPBezierCurve(curveData: leftData)
        }
        if let rightData = pap.rightCurve {
            rightBezierCurve = BPBezierCurve(curveData: rightData)
        }
        return (point: pap.point, leftBezierCurve: leftBezierCurve, rightBezierCurve: rightBezierCurve)
    }
    
    // 1535
    //- (CGFloat) refineParameter:(CGFloat)parameter forPoint:(NSPoint)point
    fileprivate func refineParameter(_ parameter: Double, forPoint point: CGPoint) -> Double {
        return BPBezierCurveData.pfRefineParameter(self.data, parameter: parameter, point: point)
    }
    
    // 1540
    //- (CGFloat) length
    func length() -> Double {
        return data.getLength()
    }
    
    // 1545
    //- (CGFloat) lengthAtParameter:(CGFloat)parameter
    func lengthAtParameter(_ parameter: Double) -> Double {
        return data.getLengthAtParameter(parameter)
    }
    
    // 1550
    //- (BOOL) isPoint
    var isPoint: Bool {
        return data.isPoint
    }
    // 1555
    //- (FBBezierCurveLocation) closestLocationToPoint:(NSPoint)point
    func closestLocationToPoint(_ point: CGPoint) -> BPBezierCurveLocation {
        return data.closestLocationToPoint(point)
    }
    
    // 1560
    //- (NSRect) bounds
    var bounds: CGRect {
        return data.bounds
    }
    
    // 1565
    //- (NSRect) boundingRect
    var boundingRect: CGRect {
        return data.boundingRect
    }
    
    func pointFromRightOffset(_ offset: Double) -> CGPoint {
        var offset = offset
        let len = length()
        offset = min(offset, len)
        let time = 1.0 - (offset / len)
        return data.pointAtParameter(time).point
    }
    
    // 1578
    //- (NSPoint) pointFromLeftOffset:(CGFloat)offset
    func pointFromLeftOffset(_ offset: Double) -> CGPoint {
        var offset = offset
        let len = length()
        offset = min(offset, len)
        let time = offset / len
        return data.pointAtParameter(time).point
    }
    
    // 1586
    //- (NSPoint) tangentFromRightOffset:(CGFloat)offset
    func tangentFromRightOffset(_ offset: Double) -> CGPoint {
        var offset = offset
        if data.isStraightLine && !data.isPoint {
            return PointMath.subtractPoint(data.endPoint1, point2: data.endPoint2)
        }
        
        if offset == 0.0 && !data.controlPoint2.equalTo(data.endPoint2) {
            return PointMath.subtractPoint(data.controlPoint2, point2: data.endPoint2)
        } else {
            let len = length()
            if offset == 0.0 {
                offset = min(1.0, len)
            }
            let time = 1.0 - (offset / len)
            let pap = data.pointAtParameter(time)
            if let curve = pap.leftCurve {
                return PointMath.subtractPoint(curve.controlPoint2, point2: curve.endPoint2)
            }
        }
        
        // nothing else worked!
        return CGPoint.zero
    }
    
    // 1607
    //- (NSPoint) tangentFromLeftOffset:(CGFloat)offset
    func tangentFromLeftOffset(_ offset: Double) -> CGPoint {
        var offset = offset
        if data.isStraightLine && !data.isPoint {
            return PointMath.subtractPoint(data.endPoint2, point2: data.endPoint1)
        }
        
        if offset == 0.0 && !data.controlPoint1.equalTo(data.endPoint1) {
            return PointMath.subtractPoint(data.controlPoint1, point2: data.endPoint1)
        } else {
            let len = length()
            if offset == 0.0 {
                offset = min(1.0, len)
            }
            let time = offset / len
            let pap = data.pointAtParameter(time)
            if let curve = pap.rightCurve {
                return PointMath.subtractPoint(curve.controlPoint2, point2: curve.endPoint2)
            }
        }
        
        // nothing else worked!
        return CGPoint.zero
    }
    
    #if canImport(Cocoa)
    // 1628
    //- (NSBezierPath *) bezierPath
    var bezierPath: NSBezierPath {
        let path = NSBezierPath()
        path.move(to: endPoint1)
        path.curve(to: endPoint2, controlPoint1: controlPoint1, controlPoint2: controlPoint2)
        return path
    }
    #endif
    
    #if canImport(UIKit)
    var bezierPath: UIBezierPath {
        let path = UIBezierPath()
        path.move(to: endPoint1)
        path.addCurve(to: endPoint2, controlPoint1: controlPoint1, controlPoint2: controlPoint2)
        return path
    }
    #endif
    
    var path: Path {
        var path = Path()
        path.move(to: endPoint1)
        path.addCurve(to: endPoint2, control1: controlPoint1, control2: controlPoint2)
        return path
    }
    
    // 1636
    //- (FBBezierCurve *) clone
    func clone() -> BPBezierCurve {
        return BPBezierCurve(curveData: data)
    }
}

extension BPBezierCurve: Equatable {
    public static func ==(lhs: BPBezierCurve, rhs: BPBezierCurve) -> Bool {
        return lhs.data.isEqual(to: rhs.data, threshold: ThresholdConstants.pointCloseness)
    }
}

extension BPBezierCurve: CustomDebugStringConvertible {
    public var debugDescription: String {
        return String(format: "<BPBezierCurve (%.18f, %.18f)-[%.18f, %.18f] curve to [%.18f, %.18f]-(%.18f, %.18f)>",
                      data.endPoint1.x, data.endPoint1.y, data.controlPoint1.x, data.controlPoint1.y,
                      data.controlPoint2.x, data.controlPoint2.y, data.endPoint2.x, data.endPoint2.y)
    }
}

extension BPBezierCurve: CustomStringConvertible {
    public var description: String {
        return "<\(data.endPoint1.x), \(data.endPoint1.y), \(data.controlPoint1.x), \(data.controlPoint1.y), \(data.controlPoint2.x), \(data.controlPoint2.y), \(data.endPoint2.x), \(data.endPoint2.y)>"
    }
}
