# BooleanPath for macOS and iOS
Add boolean operations to `Cocoa`'s `NSBezierPath`, `UIKit`'s `UIBezierPath`, and `SwiftUI`'s `Path`.

## About BooleanPath
This is a rework and update of [BooleanPath](https://github.com/Kyome22/BooleanPath) written by Kyome22, which is a rewrite of [VectorBoolean](https://github.com/lrtitze/Swift-VectorBoolean) written by Leslie Titze.  
BooleanPath is written by Swift for macOS and iOS.

## Installation
### Swift Package Manager
*(need to add instructions for Swift Package Manager)*

## Demo

The sample code is in the project.

![sample](https://github.com/Kyome22/BooleanPath/blob/master/images/sample.png)

## Usage (Example)
*(needs to be updated from original example)*

```swift
import Cocoa
import BooleanPath

let rectPath = NSBezierPath(rect: NSRect(x: 10, y: 30, width: 60, height: 60))
let circlePath = NSBezierPath(ovalIn: NSRect(x: 25, y: 15, width: 50, height: 50))
  
// Union        
let unionPath: NSBezierPath = rectPath.union(circlePath)
unionPath.fill()

// Intersection
let intersectionPath: NSBezierPath = rectPath.intersection(circlePath)
intersectionPath.fill()
        
// Difference
let differencePath: NSBezierPath = rectPath.difference(circlePath)
differencePath.fill()
        
// X_Or
let xorPath: NSBezierPath = rectPath.xor(circlePath)
xorPath.fill()
```
## Todo's

- [ ] Add documentation from original Swift-VectorBoolean.
