/*
 Contains code for generating the states needed to test and produce production
 figures for Jaynes Cummings model.
 
 In a production library would expand to contain a catalogue of standard states.
*/
 import Foundation

extension StateSpace {
    public func makeCoherentState(alpha: Complex) -> Vector {
        
        var output = Vector(in: self)
        
        let prefactor = exp( -(alpha.modulus * alpha.modulus / 2) )

        output[0] = Complex(real: prefactor)
        for n in 1 ..< self.dimension {
            // avoid computing n! as this can get problematic
            output[n] = output[n-1] * alpha / sqrt(Double(n))
        }
        return output
    }
    
    public func makeNumberState(_ n: Int) -> Vector {
        var output = Vector(in: self)
        output[n] = Complex(1)
        return output
    }
    
    public func makeVector(from input: [Double]) -> Vector {
        let output = input.map( { Complex($0) } )
        return Vector(elements: output, in: self)
    }
    
    public func makeVector(from input: [Complex]) -> Vector {
        return Vector(elements: input, in: self)
    }

}
//  Created by M J Everitt on 21/01/2022.

