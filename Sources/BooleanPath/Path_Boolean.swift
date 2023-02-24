//
//  Path_Boolean.swift
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
    func union(_ path: Path) -> Path {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = thisGraph.union(with: otherGraph).path
//        result.copyAttributesFrom(self)
        return result
    }
    
    func intersection(_ path: Path) -> Path {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = thisGraph.intersect(with: otherGraph).path
//        result.copyAttributesFrom(self)
        return result
    }
    
    func subtraction(_ path: Path) -> Path {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        let result = thisGraph.subtract(with: otherGraph).path
//        result.copyAttributesFrom(self)
        return result
    }
    
    func difference(_ path: Path) -> Path {
        let thisGraph = BPBezierGraph(path: self)
        let otherGraph = BPBezierGraph(path: path)
        var result = Path()
        result.addPath(thisGraph.subtract(with: otherGraph).path)
        result.addPath(otherGraph.subtract(with: thisGraph).path)
//        result.copyAttributesFrom(self)
        return result
    }
}
