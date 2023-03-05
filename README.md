# *Disclaimer!!*
This was meant as a way to update the original package to a Swift Package and add compatability with `UIBezierPath` and `SwiftUI`'s `Path`. After my initial reorganizing of the codebase, it isn't working. I've decided not to spend more time on this, as Swift has added boolean functions to [CGPath](https://developer.apple.com/documentation/coregraphics/cgpath) [here](https://developer.apple.com/documentation/coregraphics/cgpath#2874290) for iOS 16+/macOS 13+. Meanwhile, the *original* original library works with minimal work [found here](https://github.com/edonv/VectorBoolean/blob/master/README.md) (it works as a Swift Package). Feel free to use or fork either library as you'd like.


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

- [ ] ~~Add documentation from original Swift-VectorBoolean.~~
