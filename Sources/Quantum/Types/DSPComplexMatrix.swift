//
//  DSPDoubleComplexMatrix.swift
//  
//
//  Created by George Tait on 30/04/2024.
//

import Foundation
import Accelerate

public struct DSPComplexMatrix {
    
    var elements: DSPDoubleSplitComplex
    var space: VectorSpace
    
    public init (elements: DSPDoubleSplitComplex, in space: VectorSpace) {
        self.elements = elements
        self.space = space
    }
    
    public mutating func accelerateVectorMult(_ rhs: Vector) -> Vector {
        
        assert(space == rhs.space, "Cannot multiply matrix and vector of different spaces.")
        
        var output = Vector(in: space)
        
        let dim = vDSP_Length(space.dimension)
        let rhsNumCols = vDSP_Length(1)
        
        let rhsReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: Int(dim))
        _ = rhsReal.initialize(from: rhs.elements.map{$0.real})
        
        let rhsImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: Int(dim))
        _ = rhsImag.initialize(from: rhs.elements.map{$0.imag})
        
        var rhs_calc = DSPDoubleSplitComplex(realp: rhsReal.baseAddress!, imagp: rhsImag.baseAddress!)
        
        let outReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: Int(dim))
        let outImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: Int(dim))
        
        var outElem = DSPDoubleSplitComplex(realp: outReal.baseAddress!, imagp: outImag.baseAddress!)
        
        vDSP_zmmulD(&elements, 1, &rhs_calc, 1, &outElem, 1, dim, rhsNumCols, dim)
        
        for i in 0..<space.dimension {
            output.elements[i] = Complex(real: outElem.realp[i], imag: outElem.imagp[i])
        }
        
        return output
    }
    
    public static func + (lhs: DSPComplexMatrix, rhs: DSPComplexMatrix) -> DSPComplexMatrix {
        
        assert(lhs.space == rhs.space, "Cannot add two matrices in different spaces")
        
        let dim = lhs.space.dimension
        let outReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim*dim)
        let outImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim*dim)
        
        var outElem = DSPDoubleSplitComplex(realp: outReal.baseAddress!, imagp: outImag.baseAddress!)
        vDSP.add(lhs.elements, to: rhs.elements, count: dim*dim, result: &outElem)
        
        return DSPComplexMatrix(elements: outElem, in: lhs.space)
    }
    
    public static func - (lhs: DSPComplexMatrix, rhs: DSPComplexMatrix) -> DSPComplexMatrix {
        
        assert(lhs.space == rhs.space, "Cannot add two matrices in different spaces")
        
        let dim = lhs.space.dimension
        let outReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim*dim)
        let outImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim*dim)
        
        var outElem = DSPDoubleSplitComplex(realp: outReal.baseAddress!, imagp: outImag.baseAddress!)
        vDSP.subtract(rhs.elements, from: rhs.elements, count: dim*dim, result: &outElem)
        
        return DSPComplexMatrix(elements: outElem, in: lhs.space)
    }
    
}


