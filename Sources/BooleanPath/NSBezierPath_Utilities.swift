//
//  NSBezierPath_Utilities.swift
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


public struct NSBezierElement {
    var kind: CGPathElementType
    var point: CGPoint
    var controlPoints: [CGPoint]
}
import SwiftUI

let BPDebugPointSize: CGFloat = 10.0
let BPDebugSmallPointSize: CGFloat = 3.0

#if canImport(Cocoa)
import Cocoa

public extension NSBezierPath {
    var elements: [Path.Element] {
        var points = [CGPoint](repeating: .zero, count: 3)
        return (0..<self.elementCount).reduce(into: []) { partialResult, i in
            let element: Path.Element
            switch self.element(at: i, associatedPoints: &points) {
            case .moveTo:
                element = .move(to: CGPoint(x: points[0].x, y: points[0].y))
            case .lineTo:
                element = .line(to: CGPoint(x: points[0].x, y: points[0].y))
            case .curveTo:
                element = .curve(to: CGPoint(x: points[2].x, y: points[2].y),
                                 control1: CGPoint(x: points[0].x, y: points[0].y),
                                 control2: CGPoint(x: points[1].x, y: points[1].y))
            case .closePath:
                element = .closeSubpath
            @unknown default:
                element = .closeSubpath
            }
            
            partialResult.append(element)
        }
    }
    
    var cgPath: CGPath {
        let path: CGMutablePath = CGMutablePath()
        var points = [NSPoint](repeating: NSPoint.zero, count: 3)
        
        for i in (0 ..< self.elementCount) {
            switch self.element(at: i, associatedPoints: &points) {
            case .moveTo:
                path.move(to: CGPoint(x: points[0].x, y: points[0].y))
            case .lineTo:
                path.addLine(to: CGPoint(x: points[0].x, y: points[0].y))
            case .curveTo:
                path.addCurve(to: CGPoint(x: points[2].x, y: points[2].y),
                              control1: CGPoint(x: points[0].x, y: points[0].y),
                              control2: CGPoint(x: points[1].x, y: points[1].y))
            case .closePath:
                path.closeSubpath()
            @unknown default:
                path.closeSubpath()
            }
        }
        
        return path
    }
    
    var dividedPaths: [NSBezierPath] {
        var path: NSBezierPath?
        var paths = [NSBezierPath]()
        var points = [CGPoint](repeating: .zero, count: 3)
        for i in (0 ..< self.elementCount) {
            switch self.element(at: i, associatedPoints: &points) {
            case .moveTo:
                path = NSBezierPath()
                path?.move(to: CGPoint(x: points[0].x, y: points[0].y))
            case .lineTo:
                path?.line(to: CGPoint(x: points[0].x, y: points[0].y))
            case .curveTo:
                path?.curve(to: NSPoint(x: points[2].x, y: points[2].y),
                            controlPoint1: NSPoint(x: points[0].x, y: points[0].y),
                            controlPoint2: NSPoint(x: points[1].x, y: points[1].y))
            case .closePath:
                if let path = path {
                    paths.append(path)
                }
                path = nil
            @unknown default:
                if let path = path {
                    paths.append(path)
                }
                path = nil
            }
        }
        return paths
    }
    
    func copyAttributesFrom(_ path: NSBezierPath) {
        self.lineWidth = path.lineWidth
        self.lineCapStyle = path.lineCapStyle
        self.lineJoinStyle = path.lineJoinStyle
        self.miterLimit = path.miterLimit
        self.flatness = path.flatness
    }
    
    func callPath() {
        var points = [NSPoint](repeating: NSPoint.zero, count: 3)
        for i in (0 ..< self.elementCount) {
            switch self.element(at: i, associatedPoints: &points) {
            case .moveTo:
                Swift.print("moveTo: (\(points[0].x), \(points[0].y))")
            case .lineTo:
                Swift.print("lineTo: (\(points[0].x), \(points[0].y))")
            case .curveTo:
                Swift.print("curveTo: (\(points[2].x), \(points[2].y)), (\(points[0].x), \(points[0].y)), (\(points[1].x), \(points[1].y))")
            case .closePath:
                Swift.print("close")
            @unknown default:
                Swift.print("close")
            }
        }
    }
}

// MARK: - Static Func Inits for Debugging

public extension NSBezierPath {
    static func circleAtPoint(_ point: CGPoint) -> NSBezierPath {
        let rect = CGRect(
            x: point.x - BPDebugPointSize * 0.5,
            y: point.y - BPDebugPointSize * 0.5,
            width: BPDebugPointSize,
            height: BPDebugPointSize);
        return NSBezierPath(ovalIn: rect)
    }
    
    static func rectAtPoint(_ point: CGPoint) -> NSBezierPath {
        let rect = CGRect(
            x: point.x - BPDebugPointSize * 0.5,
            y: point.y - BPDebugPointSize * 0.5,
            width: BPDebugPointSize,
            height: BPDebugPointSize);
        return NSBezierPath(rect: rect)
    }
    
    static func smallCircleAtPoint(_ point: CGPoint) -> NSBezierPath {
        let rect = CGRect(
            x: point.x - BPDebugSmallPointSize * 0.5,
            y: point.y - BPDebugSmallPointSize * 0.5,
            width: BPDebugSmallPointSize,
            height: BPDebugSmallPointSize);
        return NSBezierPath(ovalIn: rect)
    }
    
    static func smallRectAtPoint(_ point: CGPoint) -> NSBezierPath {
        let rect = CGRect(
            x: point.x - BPDebugSmallPointSize * 0.5,
            y: point.y - BPDebugSmallPointSize * 0.5,
            width: BPDebugSmallPointSize,
            height: BPDebugSmallPointSize);
        return NSBezierPath(rect: rect)
    }
    
    static func triangleAtPoint(_ point: CGPoint, direction tangent: CGPoint) -> NSBezierPath {
        let endPoint = BPAddPoint(point, point2: BPScalePoint(tangent, scale: BPDebugPointSize * 1.5))
        let normal1 = BPLineNormal(point, lineEnd: endPoint)
        let normal2 = CGPoint(x: -normal1.x, y: -normal1.y)
        let basePoint1 = BPAddPoint(point, point2: BPScalePoint(normal1, scale: BPDebugPointSize * 0.5))
        let basePoint2 = BPAddPoint(point, point2: BPScalePoint(normal2, scale: BPDebugPointSize * 0.5))
        let path = NSBezierPath()
        path.move(to: basePoint1)
        path.line(to: endPoint)
        path.line(to: basePoint2)
        path.line(to: basePoint1)
        path.close()
        return path
    }
}
#endif


