//
//  Int+Extensions.swift
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

extension Int {
    /// Our utility extension to ease identification
    /// of places where we test for even numbers
    public var isEven: Bool {
        return self.isMultiple(of: 2)
    }
    
    /// Our utility extension to ease identification
    /// of places where we test for odd numbers
    public var isOdd: Bool {
        return !self.isMultiple(of: 2)
    }
}
