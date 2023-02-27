//
//  Path_Utilities.swift
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

public extension Path {
    var elementCount: Int {
        var count = 0
        Path(cgPath).forEach { _ in count += 1 }
        return count
    }
    
    var elements: [Path.Element] {
        var arr = [Path.Element]()
        self.forEach { el in
            arr.append(el)
        }
        return arr
    }
    
    var dividedPaths: [Path] {
        var path: Path?
        return elements.reduce(into: []) { partialResult, element in
            switch element {
            case .move(let to):
                path = Path()
                path?.move(to: to)
            case .line(let to):
                path?.addLine(to: to)
            case .curve(let to, let control1, let control2):
                path?.addCurve(to: to,
                               control1: control1,
                               control2: control2)
            case .quadCurve(let to, let control):
                path?.addQuadCurve(to: to,
                                   control: control)
            case .closeSubpath:
                if let path = path {
                    partialResult.append(path)
                }
                path = nil
            }
        }
    }
    
//    func copyAttributesFrom(_ path: Path) {
//        path.strokedPath(<#T##style: StrokeStyle##StrokeStyle#>)
//        
//        self.lineWidth = path.lineWidth
//        self.lineCapStyle = path.lineCapStyle
//        self.lineJoinStyle = path.lineJoinStyle
//        self.miterLimit = path.miterLimit
//        self.flatness = path.flatness
//    }
    
    func callPath() {
        for element in elements {
            switch element {
            case .move(let to):
                print("moveTo: (\(to.x), \(to.y))")
            case .line(let to):
                print("lineTo: (\(to.x), \(to.y))")
            case .curve(let to, let control1, let control2):
                print("curveTo: (\(to.x), \(to.y)), (\(control1.x), \(control1.y)), (\(control2.x), \(control2.y))")
            case .quadCurve(let to, let control):
                print("curveTo: (\(to.x), \(to.y)), (\(control.x), \(control.y))")
            case .closeSubpath:
                print("close")
            }
        }
    }
}

// MARK: - Static Func Inits for Debugging

public extension Path {
    static func circleAtPoint(_ point: CGPoint) -> Path {
        let rect = CGRect(
            x: point.x - BPDebugPointSize * 0.5,
            y: point.y - BPDebugPointSize * 0.5,
            width: BPDebugPointSize,
            height: BPDebugPointSize);
        return Path(ellipseIn: rect)
    }
    
    static func rectAtPoint(_ point: CGPoint) -> Path {
        let rect = CGRect(
            x: point.x - BPDebugPointSize * 0.5,
            y: point.y - BPDebugPointSize * 0.5,
            width: BPDebugPointSize,
            height: BPDebugPointSize);
        return Path(rect)
    }
    
    static func smallCircleAtPoint(_ point: CGPoint) -> Path {
        let rect = CGRect(
            x: point.x - BPDebugSmallPointSize * 0.5,
            y: point.y - BPDebugSmallPointSize * 0.5,
            width: BPDebugSmallPointSize,
            height: BPDebugSmallPointSize);
        return Path(ellipseIn: rect)
    }
    
    static func smallRectAtPoint(_ point: CGPoint) -> Path {
        let rect = CGRect(
            x: point.x - BPDebugSmallPointSize * 0.5,
            y: point.y - BPDebugSmallPointSize * 0.5,
            width: BPDebugSmallPointSize,
            height: BPDebugSmallPointSize);
        return Path(rect)
    }
    
    static func triangleAtPoint(_ point: CGPoint, direction tangent: CGPoint) -> Path {
        let endPoint = PointMath.addPoint(point, point2: PointMath.scalePoint(tangent, scale: BPDebugPointSize * 1.5))
        let normal1 = PointMath.lineNormal(point, lineEnd: endPoint)
        let normal2 = CGPoint(x: -normal1.x, y: -normal1.y)
        let basePoint1 = PointMath.addPoint(point, point2: PointMath.scalePoint(normal1, scale: BPDebugPointSize * 0.5))
        let basePoint2 = PointMath.addPoint(point, point2: PointMath.scalePoint(normal2, scale: BPDebugPointSize * 0.5))
        var path = Path()
        path.move(to: basePoint1)
        path.addLine(to: endPoint)
        path.addLine(to: basePoint2)
        path.addLine(to: basePoint1)
        path.closeSubpath()
        return path
    }

}
