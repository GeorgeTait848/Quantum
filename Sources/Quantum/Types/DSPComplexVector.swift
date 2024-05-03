//
//  DSPComplexVector.swift
//  BenchmarkQuantum
//
//  Created by George Tait on 03/05/2024.
//

import Foundation


import Foundation
import Accelerate

public struct DSPComplexVector {
    
    var elements: DSPDoubleSplitComplex
    var space: VectorSpace
    
    public init (elements: DSPDoubleSplitComplex, in space: VectorSpace) {
        self.elements = elements
        self.space = space
    }
    
    public init (fromComplexArray elem: [Complex], in space: VectorSpace) {
        let dim = space.dimension
        let real = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim)
        _ = real.initialize(from: elem.map{$0.real})
        
        let imag = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim)
        _ = imag.initialize(from: elem.map{$0.imag})
        
        let dsp_elem = DSPDoubleSplitComplex(realp: real.baseAddress!, imagp: imag.baseAddress!)
        self.elements = dsp_elem
        self.space = space
    }
    
    public subscript (row: Int) -> Complex {
        
        get {
            return Complex(real: elements.realp[row], imag: elements.imagp[row])
        }
    }
   
    
    public static func * (lhs: inout DSPComplexVector, rhs: Complex) -> DSPComplexVector {
        
        let dim = vDSP_Length(lhs.space.dimension)
        
        let right_real = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        right_real.initialize(to: rhs.real)
        
        let right_imag = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        right_imag.initialize(to: rhs.imag)
        var right = DSPDoubleSplitComplex(realp: right_real, imagp: right_imag)
        
        
        let out_real = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim))
        let out_imag = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim))
        var out_dsp = DSPDoubleSplitComplex(realp: out_real, imagp: out_imag)
        
        vDSP_zvzsmlD(&lhs.elements, 1, &right, &out_dsp, 1, dim)
        
        
        return DSPComplexVector(elements: out_dsp, in: lhs.space)
    }
    
    public static func * (lhs: Complex, rhs: inout DSPComplexVector) -> DSPComplexVector {
        
        let dim = vDSP_Length(rhs.space.dimension)
        
        let left_real = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        left_real.initialize(to: lhs.real)
        
        let left_imag = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        left_imag.initialize(to: lhs.imag)
        var left = DSPDoubleSplitComplex(realp: left_real, imagp: left_imag)
        
        
        let out_real = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim))
        let out_imag = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim))
        var out_dsp = DSPDoubleSplitComplex(realp: out_real, imagp: out_imag)
        
        vDSP_zvzsmlD(&rhs.elements, 1, &left, &out_dsp, 1, dim)
        
        return DSPComplexVector(elements: out_dsp, in: rhs.space)
    }
    
    public static func * (lhs: inout DSPComplexVector, rhs: Double) -> DSPComplexVector {
        
        return lhs * Complex(real: rhs)
    }
    
    public static func * (lhs: Double, rhs: inout DSPComplexVector) -> DSPComplexVector {
        
        return Complex(real: lhs) * rhs
    }
    
    
    
    public static func / (lhs: inout DSPComplexVector, rhs: Complex) -> DSPComplexVector {
        
        let norm = rhs.norm
        
        let div_real = rhs.real/norm
        let div_imag = -rhs.imag/norm
        
        let dim = vDSP_Length(lhs.space.dimension)
        
        let div_real_dsp = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        div_real_dsp.initialize(to: div_real)
        
        let div_imag_dsp = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        div_imag_dsp.initialize(to: div_imag)
        var div = DSPDoubleSplitComplex(realp: div_real_dsp, imagp: div_imag_dsp)
        
        
        let out_real = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim))
        let out_imag = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim))
        var out_dsp = DSPDoubleSplitComplex(realp: out_real, imagp: out_imag)
        
        vDSP_zvzsmlD(&lhs.elements, 1, &div, &out_dsp, 1, dim)
        
        return DSPComplexVector(elements: out_dsp, in: lhs.space)
    }
    
    public static func / (lhs: inout DSPComplexVector, rhs: Double) -> DSPComplexVector {
        
        let dim = vDSP_Length(lhs.space.dimension)
        
        let div_real_dsp = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        div_real_dsp.initialize(to: 1.0/rhs)
        
        let div_imag_dsp = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        div_imag_dsp.initialize(to: 0.0)
        var div = DSPDoubleSplitComplex(realp: div_real_dsp, imagp: div_imag_dsp)
        
        
        let out_real = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim))
        let out_imag = UnsafeMutablePointer<Double>.allocate(capacity: Int(dim))
        var out_dsp = DSPDoubleSplitComplex(realp: out_real, imagp: out_imag)
        
        vDSP_zvzsmlD(&lhs.elements, 1, &div, &out_dsp, 1, dim)
        
        return DSPComplexVector(elements: out_dsp, in: lhs.space)
        
    }
    

    
    public mutating func innerProduct(_ rhs: inout DSPComplexVector) -> Complex {
        
        assert(space == rhs.space, "Cannot take inner product of vectors in different spaces.")
        
        let out_real = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        let out_imag = UnsafeMutablePointer<Double>.allocate(capacity: 1)
        var out_dsp = DSPDoubleSplitComplex(realp: out_real, imagp: out_imag)
        
        var output = Complex(real: 0.0)
        
        let dim = vDSP_Length(space.dimension)
        
   
        vDSP_zidotprD(&elements, 1, &rhs.elements, 1, &out_dsp, dim)
        output.real = out_real.pointee
        output.imag = out_imag.pointee
        
        return output
    }
    
    public static func + (lhs: DSPComplexVector, rhs: DSPComplexVector) -> DSPComplexVector {
        
        assert(lhs.space == rhs.space, "Cannot add two matrices in different spaces")
        
        let dim = lhs.space.dimension
        let outReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim)
        let outImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim)
        
        var outElem = DSPDoubleSplitComplex(realp: outReal.baseAddress!, imagp: outImag.baseAddress!)
        vDSP.add(lhs.elements, to: rhs.elements, count: dim, result: &outElem)
        
        return DSPComplexVector(elements: outElem, in: lhs.space)
    }
    
    public static func - (lhs: DSPComplexVector, rhs: DSPComplexVector) -> DSPComplexVector {
        
        assert(lhs.space == rhs.space, "Cannot add two matrices in different spaces")
        
        let dim = lhs.space.dimension
        let outReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim)
        let outImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: dim)
        
        var outElem = DSPDoubleSplitComplex(realp: outReal.baseAddress!, imagp: outImag.baseAddress!)
        
        /*
         By the argument names, one would assume subtracting lhs from rhs is rhs-lhs,
         but is in fact other way around. See CoreTests.swift to see this works.
        */
        vDSP.subtract(lhs.elements, from: rhs.elements, count: dim, result: &outElem)
        
        return DSPComplexVector(elements: outElem, in: lhs.space)
    }
    
}


extension DSPComplexVector: Equatable {
    
    public static func == (lhs: DSPComplexVector, rhs: DSPComplexVector) -> Bool {
        
        assert(lhs.space == rhs.space, "Cannot equate two vectors in different spaces")
        
        let dim = lhs.space.dimension
        
        for i in 0..<dim {
            
            if lhs.elements.realp[i] != rhs.elements.realp[i] { return false }
            if lhs.elements.imagp[i] != rhs.elements.imagp[i] { return false }
            
        }
        
        return true
    }
}

extension DSPComplexVector: Addable {}

extension DSPComplexVector: CustomStringConvertible {
    public var description: String {
        var output = "Vector in space \(space.description) with identifier: \(space.identifier)\n"
        for i in 0..<space.dimension {
            output.append("\(elements.realp[i]) + \(elements.imagp[i]) i\n")
        }
        return output
    }
}

