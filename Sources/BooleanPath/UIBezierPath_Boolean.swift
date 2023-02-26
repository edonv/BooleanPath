//
//  UIBezierPath_Boolean.swift
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

#if canImport(UIKit)
import UIKit

public extension UIBezierPath {
    func union(_ path: UIBezierPath) -> UIBezierPath {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = thisGraph.union(with: otherGraph).bezierPath
        result.copyAttributesFrom(self)
        return result
    }
    
    func intersection(_ path: UIBezierPath) -> UIBezierPath {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = thisGraph.intersect(with: otherGraph).bezierPath
        result.copyAttributesFrom(self)
        return result
    }
    
    func subtraction(_ path: UIBezierPath) -> UIBezierPath {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = thisGraph.difference(with: otherGraph).bezierPath
        result.copyAttributesFrom(self)
        return result
    }
    
    func difference(_ path: UIBezierPath) -> UIBezierPath {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = UIBezierPath()
        result.append(thisGraph.difference(with: otherGraph).bezierPath)
        result.append(otherGraph.difference(with: thisGraph).bezierPath)
        result.copyAttributesFrom(self)
        return result
    }
}
#endif
