/*
 See notes in Spaces file
 
 TODO: Optimise
 */

import Foundation

public struct Vector: VectorType {
   
    

    public var space: VectorSpace
    public var elements: [Complex]
    
    public init(elements: [Complex], in space: VectorSpace) {
        assert(elements.count == space.dimension, "Array is not the same dimension as the space")
        self.elements = elements
        self.space = space
    }
    
    public init(elements: [Double], in space: VectorSpace) {
        assert(elements.count == space.dimension, "Array is not the same dimension as the space")
        self.elements = elements.map{Complex($0)}
        self.space = space
    }

    
    
    public init(in space: VectorSpace) {
        self.elements = Array(repeating: Complex(real: 0.0), count: space.dimension)
        self.space = space
    }
    
    
    public subscript(index: Int) -> Complex {
        get {
            assert(index < elements.count, "Index out of range getting vector value")
            return elements[index]
        }
        set {
            assert(index < elements.count, "Index out of range setting vector value")
            elements[index] = newValue
        }
    }
    // Closed under scalar mutliplication w.r.t the scalar field over which it is defined
    public static func * (left: Self, right: Complex) -> Self {
        return Self(elements: scalarBinaryOperation(array: left.elements, by: right,operation: *), in: left.space)
    }
    
    public static func * (left: Self, right: Double) -> Self {
        return Self(elements: scalarBinaryOperation(array: left.elements, by: Complex(real: right), operation: *), in: left.space)
    }
    public static func * (left: Complex, right: Self) -> Self {
        return Self(elements: scalarBinaryOperation(array: right.elements, by: left,operation: *), in: right.space)
    }
    
    public static func * (left: Double, right: Self) -> Self {
        return Self(elements: scalarBinaryOperation(array: right.elements, by: Complex(real: left),operation: *), in: right.space)
    }
    
    public static func / (left: Self, right: Complex) -> Self {
        return Self(elements: scalarBinaryOperation(array: left.elements, by: right,operation: /), in: left.space)
    }
    
    public static func / (left: Self, right: Double) -> Self {
        return Self(elements: scalarBinaryOperation(array: left.elements, by: Complex(real: right),operation: /), in: left.space)
    }
}

extension Vector: CustomStringConvertible {
    public var description: String {
        var output = "Vector in space \(space.description) with identifier: \(space.identifier)\n"
        for value in elements {
            output.append("\(value)\n")
        }
        return output
    }
}

extension Vector  {
    public func innerProduct(dualVector: Self) -> Complex {
        var sum = Complex(real: 0.0)
        checkInSameSpace(self, dualVector)
        for i in 0 ..< elements.count {
            sum =  sum + elements[i] *  dualVector[i].conjugate
        }
        return sum
    }
}

extension Vector {
    public func outerProduct(with v: Self) -> Matrix {
        assert (self.space == v.space)
        var output = Matrix(in: self.space)
        for i in 0 ..< elements.count {
            for j in 0 ..< elements.count {
                output[i,j] = self[i] * v[j].conjugate
            }
        }
        return output
    }
}

extension Vector: Equatable {
    public static func == (lhs: Vector, rhs: Vector) -> Bool {
        assert(lhs.space == rhs.space, "Cannot equate two vectors in different spaces")
        
        for i in 0..<lhs.space.dimension {
            if lhs.elements[i] == rhs.elements[i] { continue }
            return false
        }
        
        return true
    }
}

extension Vector: AdaptiveSteppable {
    public func checkEquivalenceUnderRelativeTolerance(_ other: Vector, relativeTol: Double) -> Bool {
        assert(space == other.space, "Cannot check equivalence of vectors in different spaces")
        
        for i in 0..<space.dimension {
            let lhs = elements[i]
            let rhs = other.elements[i]
            
            if (lhs-rhs).modulus > relativeTol { return false}
        }
        
        return true
    }
    
    
}

//  Created by M J Everitt on 18/01/2022.
