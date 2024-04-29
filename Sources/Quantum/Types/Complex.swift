/*
 See ComplexNumber protocol for much of the functionality.
 
 The idea of using protocols was to enable more future flexibility.
 This may have been unnecessary on reflection.
 
 Note: struct goes directly on the stack so as well as been data
 type (avoiding sharing of state issues) so it avoids the heap.
 */
import Foundation


public struct Complex: ComplexNumber, Hashable {
    
    
    public init(_ r: Double) {
        self.init(real: r)
    }
    
    public init(real: Double) {
        self.real = real
        self.imag = 0.0
    }

    public var real: Double
    public var imag: Double
        
    public init(real: Double, imag: Double) {
        self.real = real
        self.imag = imag
    }
    
    public init (modulus: Double, argument: Double) {
        self.init(real: modulus * cos(argument),
                  imag: modulus * sin(argument))

    }
}




extension Complex: CustomStringConvertible {
    public var description: String {
        return "\(real) + \(imag) i"
    }
}

//sort this shit later

//extension Complex: Has_Exp {
//    public static func exp(_ x: Self) -> Self {
//        let realExp: Double = exp(x.real)
//        return( Self(real: realExp * cos(x.imag),
//                     imag: realExp * sin(x.imag)) )
//    }
//}
//extension Complex: Has_Sqrt where T: Has_Exp & Has_Cos & Has_Sin & Has_Atan & Has_Sqrt {
//    public static func sqrt(_ x: Self) -> Self {
//        return Self(modulus: T.sqrt(x.modulus), argument: x.argument / T(2) )
//    }
//}
extension Complex {
    public static func abs(_ x: Self) -> Double {
        return x.modulus
    }
}
//  Created by M J Everitt on 16/01/2022.
