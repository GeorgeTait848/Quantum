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
    
    public init (fromComplexArray elem: [Complex], in space: VectorSpace) {
    
        let dim = space.dimension
        let real = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim*dim)
        _ = real.initialize(from: elem.map{$0.real})
        
        let imag = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim*dim)
        _ = imag.initialize(from: elem.map{$0.imag})
        
        let dsp_elem = DSPDoubleSplitComplex(realp: real.baseAddress!, imagp: imag.baseAddress!)
        self.elements = dsp_elem
        self.space = space
        
    }
    
    public static func * (lhs: inout DSPComplexMatrix, rhs: inout DSPComplexVector) -> DSPComplexVector {
        
        assert(lhs.space == rhs.space, "Cannot multiply matrix and vector of different spaces.")
        
        let dim = vDSP_Length(lhs.space.dimension)
        let rhsNumCols = vDSP_Length(1)
        
        let outReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: Int(dim))
        let outImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: Int(dim))
        
        var outElem = DSPDoubleSplitComplex(realp: outReal.baseAddress!, imagp: outImag.baseAddress!)
        
        vDSP_zmmulD(&lhs.elements, 1, &rhs.elements, 1, &outElem, 1, dim, rhsNumCols, dim)
        
        return DSPComplexVector(elements: outElem, in: lhs.space)
        
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
    
    public static func * (lhs: inout DSPComplexMatrix, rhs: Complex) -> DSPComplexMatrix {
        
        let dim = vDSP_Length(lhs.space.dimension)
        
        let right_real = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        right_real.initialize(to: rhs.real)
        
        let right_imag = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        right_imag.initialize(to: rhs.imag)
        var right = DSPDoubleSplitComplex(realp: right_real, imagp: right_imag)
        
        
        let out_real = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim*dim))
        let out_imag = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim*dim))
        var out_dsp = DSPDoubleSplitComplex(realp: out_real, imagp: out_imag)
        
        vDSP_zvzsmlD(&lhs.elements, 1, &right, &out_dsp, 1, dim*dim)
        
        
        return DSPComplexMatrix(elements: out_dsp, in: lhs.space)
    }
    
    public static func * (lhs: Complex, rhs: inout DSPComplexMatrix) -> DSPComplexMatrix {
        
        return rhs * lhs
    }
    
    public static func * (lhs: inout DSPComplexMatrix, rhs: Double) -> DSPComplexMatrix {
        
        return lhs * Complex(real: rhs)
    }
    
    public static func * (lhs: Double, rhs: inout DSPComplexMatrix) -> DSPComplexMatrix {
        
        return Complex(real: lhs) * rhs
    }
    
    public static func / (lhs: inout DSPComplexMatrix, rhs: Complex) -> DSPComplexMatrix {
        
        let norm = rhs.norm
        
        let div_real = rhs.real/norm
        let div_imag = -rhs.imag/norm
        
        let dim = vDSP_Length(lhs.space.dimension)
        
        let div_real_dsp = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        div_real_dsp.initialize(to: div_real)
        
        let div_imag_dsp = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        div_imag_dsp.initialize(to: div_imag)
        
        var div = DSPDoubleSplitComplex(realp: div_real_dsp, imagp: div_imag_dsp)
        
        let out_real = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim*dim))
        let out_imag = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim*dim))
        var out_dsp = DSPDoubleSplitComplex(realp: out_real, imagp: out_imag)
        
        vDSP_zvzsmlD(&lhs.elements, 1, &div, &out_dsp, 1, dim*dim)
        
        return DSPComplexMatrix(elements: out_dsp, in: lhs.space)
    }
    
    public static func / (lhs: inout DSPComplexMatrix, rhs: Double) -> DSPComplexMatrix {
        
        let dim = vDSP_Length(lhs.space.dimension)
        
        let div_real_dsp = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        div_real_dsp.initialize(to: 1.0/rhs)
        
        let div_imag_dsp = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        div_imag_dsp.initialize(to: 0.0)
        var div = DSPDoubleSplitComplex(realp: div_real_dsp, imagp: div_imag_dsp)
        
        
        let out_real = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim*dim))
        let out_imag = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim*dim))
        var out_dsp = DSPDoubleSplitComplex(realp: out_real, imagp: out_imag)
        
        vDSP_zvzsmlD(&lhs.elements, 1, &div, &out_dsp, 1, dim*dim)
        
        return DSPComplexMatrix(elements: out_dsp, in: lhs.space)
        
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


