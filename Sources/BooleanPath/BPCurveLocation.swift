//
//  BPCurveLocation.swift
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

import Foundation

class BPCurveLocation {
    var graph: BPBezierGraph?
    var contour: BPBezierContour?
    fileprivate var _edge: BPBezierCurve
    fileprivate var _parameter: Double
    fileprivate var _distance: Double
    
    init(edge: BPBezierCurve, parameter: Double, distance: Double) {
        _edge = edge
        _parameter = parameter
        _distance = distance
    }
    
    var edge: BPBezierCurve {
        return _edge
    }
    
    var parameter: Double {
        return _parameter
    }
    
    var distance: Double {
        return _distance
    }
}

