//
//  BPBezierCurve_Edge.swift
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

private enum EdgeMath {
    typealias LRCurvePair = (leftCurve: BPBezierCurve, rightCurve: BPBezierCurve)
    
    static func findEdge1TangentCurves(_ edge: BPBezierCurve, intersection: BPBezierIntersection) -> LRCurvePair {
        var leftCurve: BPBezierCurve, rightCurve: BPBezierCurve
        if intersection.isAtStartOfCurve1 {
            leftCurve = edge.previousNonpoint
            rightCurve = edge
        } else if intersection.isAtStopOfCurve1 {
            leftCurve = edge
            rightCurve = edge.nextNonpoint
        } else {
            leftCurve = intersection.curve1LeftBezier
            rightCurve = intersection.curve1RightBezier
        }
        return (leftCurve: leftCurve, rightCurve: rightCurve)
    }
    
    static func findEdge2TangentCurves(_ edge: BPBezierCurve, intersection: BPBezierIntersection) -> LRCurvePair {
        var leftCurve: BPBezierCurve, rightCurve: BPBezierCurve
        
        if intersection.isAtStartOfCurve2 {
            leftCurve = edge.previousNonpoint
            rightCurve = edge
        } else if intersection.isAtStopOfCurve2 {
            leftCurve = edge
            rightCurve = edge.nextNonpoint
        } else {
            leftCurve = intersection.curve2LeftBezier
            rightCurve = intersection.curve2RightBezier
        }
        
        return (leftCurve: leftCurve, rightCurve: rightCurve)
    }
    
    static func computeEdgeTangents(_ leftCurve: BPBezierCurve, rightCurve: BPBezierCurve, offset: Double, edgeTangents: inout BPTangentPair) {
        edgeTangents.left = leftCurve.tangentFromRightOffset(offset)
        edgeTangents.right = rightCurve.tangentFromLeftOffset(offset)
    }
    
    static func computeEdge1RangeTangentCurves(_ edge: BPBezierCurve, intersectRange: BPBezierIntersectRange) -> LRCurvePair {
        var leftCurve: BPBezierCurve, rightCurve: BPBezierCurve
        if intersectRange.isAtStartOfCurve1 {
            leftCurve = edge.previousNonpoint
        } else {
            leftCurve = intersectRange.curve1LeftBezier
        }
        if intersectRange.isAtStopOfCurve1 {
            rightCurve = edge.nextNonpoint
        } else {
            rightCurve = intersectRange.curve1RightBezier
        }
        return (leftCurve: leftCurve, rightCurve: rightCurve)
    }
    
    static func computeEdge2RangeTangentCurves(_ edge: BPBezierCurve, intersectRange: BPBezierIntersectRange) -> LRCurvePair {
        var leftCurve: BPBezierCurve, rightCurve: BPBezierCurve
        
        if intersectRange.isAtStartOfCurve2 {
            leftCurve = edge.previousNonpoint
        } else {
            leftCurve = intersectRange.curve2LeftBezier
        }
        if intersectRange.isAtStopOfCurve2 {
            rightCurve = edge.nextNonpoint
        } else {
            rightCurve = intersectRange.curve2RightBezier
        }
        return (leftCurve: leftCurve, rightCurve: rightCurve)
    }
}

extension BPBezierCurve {
    func addCrossing(_ crossing: BPEdgeCrossing) {
        crossing.edge = self
        crossings.append(crossing)
        sortCrossings()
    }
    
    func removeCrossing(_ crossing: BPEdgeCrossing) {
        for (index, element) in crossings.enumerated() {
            if element === crossing {
                crossings.remove(at: index)
                break
            }
        }
        sortCrossings()
    }
    
    func removeAllCrossings() {
        crossings.removeAll()
    }
    
    var next: BPBezierCurve {
        var nxt: BPBezierCurve = self
        
        if let contour = contour {
            if contour.edges.count > 0 {
                let nextIndex = index + 1
                if nextIndex >= contour.edges.count {
                    nxt = contour.edges.first!
                } else {
                    nxt = contour.edges[nextIndex]
                }
            }
        }
        return nxt
    }
    
    var previous: BPBezierCurve {
        var prev: BPBezierCurve = self
        
        if let contour = contour {
            if contour.edges.count > 0 {
                if index == 0 {
                    prev = contour.edges.last!
                } else {
                    prev = contour.edges[index - 1]
                }
            }
        }
        return prev
    }
    
    var nextNonpoint: BPBezierCurve {
        var edge = self.next
        while edge.isPoint {
            edge = edge.next
        }
        return edge
    }
    
    var previousNonpoint: BPBezierCurve {
        var edge = self.previous
        while edge.isPoint {
            edge = edge.previous
        }
        return edge
    }
    
    var hasCrossings: Bool {
        return !crossings.isEmpty
    }
    
    func crossingsWithBlock(_ block: (_ crossing: BPEdgeCrossing) -> (setStop: Bool, stopValue:Bool)) {
        for crossing in crossings {
            let (set, val) = block(crossing)
            if set && val {
                break
            }
        }
    }
    
    func crossingsCopyWithBlock(_ block: (_ crossing: BPEdgeCrossing) -> (setStop: Bool, stopValue:Bool)) {
        let crossingsCopy = crossings
        for crossing in crossingsCopy {
            let (set, val) = block(crossing)
            if set && val {
                break
            }
        }
    }
    
    func nextCrossing(_ crossing: BPEdgeCrossing) -> BPEdgeCrossing? {
        if crossing.index < crossings.count - 1 {
            return crossings[crossing.index + 1]
        } else {
            return nil
        }
    }
    
    func previousCrossing(_ crossing: BPEdgeCrossing) -> BPEdgeCrossing? {
        if crossing.index > 0 {
            return crossings[crossing.index - 1]
        } else {
            return nil
        }
    }
    
    func intersectingEdgesWithBlock(_ block: (_ intersectingEdge: BPBezierCurve) -> Void) {
        
        crossingsWithBlock() {
            (crossing: BPEdgeCrossing) -> (setStop: Bool, stopValue:Bool) in
            if !crossing.isSelfCrossing {
                if let crossingCounterpartEdge = crossing.counterpart?.edge {
                    block(crossingCounterpartEdge)
                }
            }
            return (false, false)
        }
    }
    
    func selfIntersectingEdgesWithBlock(_ block: (_ intersectingEdge: BPBezierCurve) -> Void) {
        crossingsWithBlock() {
            (crossing: BPEdgeCrossing) -> (setStop: Bool, stopValue:Bool) in
            
            if crossing.isSelfCrossing {
                if let crossingCounterpartEdge = crossing.counterpart?.edge {
                    block(crossingCounterpartEdge)
                }
            }
            return (false, false)
        }
    }
    
    var firstCrossing: BPEdgeCrossing? {
        return crossings.first
    }
    
    var lastCrossing: BPEdgeCrossing? {
        return crossings.last
    }
    
    var firstNonselfCrossing: BPEdgeCrossing? {
        var first = firstCrossing
        while first != nil && first!.isSelfCrossing {
            first = first?.next
        }
        return first
    }
    
    var lastNonselfCrossing: BPEdgeCrossing? {
        var last = lastCrossing
        while last != nil && last!.isSelfCrossing {
            last = last?.previous
        }
        return last
    }
    
    var hasNonselfCrossings: Bool {
        for crossing in crossings {
            if !crossing.isSelfCrossing {
                return true
            }
        }
        return false
    }
    
    func crossesEdge(_ edge2: BPBezierCurve, atIntersection intersection: BPBezierIntersection) -> Bool {
        // Continue if intersection isn't tangent. If it's tangent, then it can't be true anyway.
        guard !intersection.isTangent else { return false }
        
        // Continue if intersection is at the end point of the curve. If it's not, then it can't be false anyway.
        guard intersection.isAtEndPointOfCurve else { return true }
        
        var edge1Tangents = BPTangentPair(left: .zero, right: .zero)
        var edge2Tangents = BPTangentPair(left: .zero, right: .zero)
        var offset = 0.0
        
        let (edge1LeftCurve, edge1RightCurve) = EdgeMath.findEdge1TangentCurves(self, intersection: intersection)
        let edge1Length = min(edge1LeftCurve.length(), edge1RightCurve.length())
        
        let (edge2LeftCurve, edge2RightCurve) = EdgeMath.findEdge2TangentCurves(edge2, intersection: intersection)
        let edge2Length = min(edge2LeftCurve.length(), edge2RightCurve.length())
        
        let maxOffset = min(edge1Length, edge2Length)
        
        repeat {
            EdgeMath.computeEdgeTangents(edge1LeftCurve, rightCurve: edge1RightCurve, offset: offset, edgeTangents: &edge1Tangents)
            EdgeMath.computeEdgeTangents(edge2LeftCurve, rightCurve: edge2RightCurve, offset: offset, edgeTangents: &edge2Tangents)
            
            offset += 1.0
        } while TangentMath.areAmbigious(edge1Tangents, edge2Tangents: edge2Tangents) && offset < maxOffset
        
        return TangentMath.tangentsCross(edge1Tangents, edge2Tangents: edge2Tangents)
    }
    
    func crossesEdge(_ edge2: BPBezierCurve, atIntersectRange intersectRange: BPBezierIntersectRange) -> Bool {
        var edge1Tangents = BPTangentPair(left: CGPoint.zero, right: CGPoint.zero)
        var edge2Tangents = BPTangentPair(left: CGPoint.zero, right: CGPoint.zero)
        var offset = 0.0
        
        let (edge1LeftCurve, edge1RightCurve) = EdgeMath.computeEdge1RangeTangentCurves(self, intersectRange: intersectRange)
        
        let edge1Length = min(edge1LeftCurve.length(), edge1RightCurve.length())
        
        let (edge2LeftCurve, edge2RightCurve) = EdgeMath.computeEdge2RangeTangentCurves(edge2, intersectRange: intersectRange)
        let edge2Length = min(edge2LeftCurve.length(), edge2RightCurve.length())
        
        let maxOffset = min(edge1Length, edge2Length);
        
        repeat {
            EdgeMath.computeEdgeTangents(edge1LeftCurve, rightCurve: edge1RightCurve, offset: offset, edgeTangents: &edge1Tangents)
            EdgeMath.computeEdgeTangents(edge2LeftCurve, rightCurve: edge2RightCurve, offset: offset, edgeTangents: &edge2Tangents)
            
            offset += 1.0
        } while TangentMath.areAmbigious(edge1Tangents, edge2Tangents: edge2Tangents) && offset < maxOffset
        
        return TangentMath.tangentsCross(edge1Tangents, edge2Tangents: edge2Tangents);
    }
    
    // MARK: Private funcs
    
    fileprivate func sortCrossings() {
        crossings.sort(by: { $0.order < $1.order })
        for (index, crossing) in crossings.enumerated() {
            crossing.index = index
        }
    }
}
