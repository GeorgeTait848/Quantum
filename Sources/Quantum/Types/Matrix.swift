/*
 See notes in Spaces file
 
 Note that a matrix must live in a vector space.
 
 /*
  In the accelerate branch, we remove generics arising from Scalar protocol.
  The motivation behind this protocol was, scalar is a concept and not subject to precision (i.e one should have
  the ability to use Float and Double interchangably. However, using the vDSP api of Accelerate requires the specific
  precision that comes with using concrete data types. The increase in performance from using this library supercedes the
  motivation for scalar protocol. If the Accelerate framework can be shown to have significant performance increases in this code,
  these changes will be used in the main branch.
  
  To maintain quantum functionality, all vectors and operators will be assumed to be complex, and of double precision. 
  */
 TODO: Optimise
 KEY:
 DCMM - divide and conquer matrix multiplication
 BFMM - brute force matrix multiplication
 SMM - strassen's matrix multiplication
 */
import Foundation
import Accelerate

// Importantly (Liskov) Matrix can be a VectorType but it should not extend Vector.
public struct Matrix: VectorType, OperatorType {

    public static func == (lhs: Matrix, rhs: Matrix) -> Bool {
        
        assert(lhs.space == rhs.space, "Two matrices in different spaces cannot be equal")
        for row in 0..<lhs.space.dimension {
            for col in 0..<lhs.space.dimension {
                if lhs[row,col] == rhs[row,col] { continue }
                return false
            }
        }
        return true
    }
    
    // to satisfy collection for use with integrators.
    public subscript(position: Int) -> Complex {
        return elements[position]
    }
    
    public var space: VectorSpace
    public var elements: [Complex]
    
    public init(elements: [Complex], in space: VectorSpace) {
        assert(elements.count == (space.dimension * space.dimension), "Matrix is not the same dimension as the space")
        self.elements = elements
        self.space = space
    }
    
    public init(elements: [Double], in space: VectorSpace) {
        assert(elements.count == (space.dimension * space.dimension), "Matrix is not the same dimension as the space")
        self.elements = elements.map{Complex($0)}
        self.space = space
    }
    public init(in space: VectorSpace) {
        self.elements = Array(repeating: Complex(real: 0.0), count: space.dimension * space.dimension)
        self.space = space
    }


    public subscript(row: Int, col: Int) -> Complex {
        
        get {
            let index = atIndex(row: row, column: col, nColumns: space.dimension )
            assert(index < elements.count, "Index out of range getting vector value")
            return elements[index]
        }
        set {
            let index = atIndex(row: row, column: col, nColumns: space.dimension )
            assert(index < elements.count, "Index out of range setting vector value")
            elements[index] = newValue
        }
    }
    
    func indexIsValid(row: Int, column: Int) -> Bool {
        return row >= 0 && row < space.dimension && column >= 0 && column < space.dimension
    }
    
    public static func * (lhs: Self, rhs: Self) -> Self {

        return lhs.bruteForceMatrixMultiplication(rhs)
    }
    
    
    public func bruteForceMatrixMultiplication(_ rhs: Matrix) -> Matrix {
        
        assert(rhs.space == self.space, "Cannot multiply two matrices of different spaces")
        
        var output = Self(in: self.space)
        let dim = self.space.dimension

        for i in 0 ..< dim {
            for j in 0 ..< dim {
                for k in 0 ..< dim {
                    output[i, j] = output[i, j] + self[i, k] * rhs[k, j]
                }
            }
        }
        return output
    }
    
    private func accelerateMatrixMultiplication(_ rhs: Matrix) -> Matrix {
        var output = Self(in: self.space)
        
        
        return output
    }

    
    
    public static func * (lhs: Matrix, rhs: Vector) -> Vector {
        checkInSameSpace(lhs,rhs)

        var output = Vector(in: lhs.space)
        // there are much more efficent ways to do this - coding for calrity
        for i in 0 ..< lhs.space.dimension {
            for j in 0 ..< lhs.space.dimension {
                output[i] = output[i] + lhs[i, j] * rhs[j]
            }
        }
        return output
    }
    
    //The below functions need to be replaced using accelerate
    //different functions for double and complex
    // Closed under scalar mutliplication w.r.t the scalar field over which it is defined
    public static func * (left: Self, right: Complex) -> Self {
        return Self(elements: scalarBinaryOperation(array: left.elements, by: right,operation: *), in: left.space)
    }
    
    public static func * (left: Matrix, right: Double) -> Matrix {
        return Self(elements: scalarBinaryOperation(array: left.elements, by: Complex(real: right),operation: *), in: left.space)
    }
    
    
    public static func * (left: Complex, right: Self) -> Self {
        return Self(elements: scalarBinaryOperation(array: right.elements, by: left,operation: *), in: right.space)
    }
    
    public static func * (left: Double, right: Matrix) -> Matrix {
        return Self(elements: scalarBinaryOperation(array: right.elements, by: Complex(real: left),operation: *), in: right.space)
    }
    
    
    public static func / (left: Self, right: Complex) -> Self {
        return Self(elements: scalarBinaryOperation(array: left.elements, by: right,operation: /), in: left.space)
    }
    
    public static func / (left: Matrix, right: Double) -> Matrix {
        return Self(elements: scalarBinaryOperation(array: left.elements, by: Complex(real: right),operation: /), in: left.space)
    }
    
    
}

extension Matrix: CustomStringConvertible {
    public var description: String {
        var output = "Operator in space \(space.description) with identifier: \(space.identifier))\n"
        for i in 0 ..< self.space.dimension {
            for j in 0 ..< self.space.dimension {
                output.append("\(self[i, j])")
                if (j < self.space.dimension - 1) {
                    output.append(" , ")
                }
            }
            output.append("\n")
        }
        return output
    }
}
extension Matrix  {
    public func transpose () -> Matrix {
        let dim = space.dimension
        var output = Self(in: self.space)
        for i in 0 ..< dim {
            for j in 0 ..< dim {
                output[j, i] = self[i, j]
            }
        }
        return output
    }
}


infix operator =~= : ComparisonPrecedence
extension Matrix {
    
    public static func =~= (lhs: Matrix, rhs: Matrix) -> Bool {
        return lhs.approximateEquals(rhs, testPrecision: 1.0e-6 )
    }
    public func approximateEquals(_ other: Matrix, testPrecision: Double) -> Bool {
        let spacetest = self.space == other.space
        var valuestest = true
        for i in 0 ..< self.elements.count {
            let difference = self.elements[i] - other.elements[i]
            let temptest = difference.modulus < testPrecision
            valuestest = valuestest && temptest
        }
        
        return spacetest && valuestest
    }
    
}



extension Matrix {
    public func hermitianAdjoint () -> Matrix {
        let dim = space.dimension
        var output = Self(in: self.space)
        for i in 0 ..< dim {
            for j in 0 ..< dim {
                output[j, i] = self[i, j].conjugate
            }
        }
        return output
    }
}
//  Created by M J Everitt on 18/01/2022.

