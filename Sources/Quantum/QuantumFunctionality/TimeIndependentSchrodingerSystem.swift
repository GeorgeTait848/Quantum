// An example of a class that can solve the time dependent Schrödinger equation for a time independent Hamiltonian (for sparse or dense cases).

import Foundation
open class TimeIndependentSchrodingerSystem {
    public var Psi: Vector
    public var minus_i_H: Matrix
    public var time: Double
    
    public init(initialstate: Vector,
                hamiltonian: Matrix) {
        let minusi = Complex( real: 0.0, imag: -1.0 )
        time = 0.0
        Psi = initialstate
        minus_i_H = minusi * hamiltonian
    }
    
    open func schrodingerEquation(psi: Vector, time: Double)-> Vector {
            
        if diagonalSparse { return diagSparse_minus_i_H! * psi }
        
        if sparse { return sparse_minus_i_H! * psi }
        
        return minus_i_H * psi
            
    }
    
    
    public var diagonalSparse = false
    public var diagSparse_minus_i_H: DiagonalSparseMatrix?
    
    
    public var sparse = false
    public var sparse_minus_i_H: SparseMatrix?
    
    
    open func useSparseAlgebra() {
        diagonalSparse = false
        
        sparse = true
        sparse_minus_i_H = SparseMatrix(from: minus_i_H)
    }
    
    open func useDiagonalSparseAlgebra() {
        diagonalSparse = true
        diagSparse_minus_i_H = DiagonalSparseMatrix(from: minus_i_H)
    }
    
    
    public func useNonSparseAlgebra() {
        sparse = false
        diagonalSparse = false
    }

    public func evolve(by dt: Double) {
        
//        multiStepIvpIntegrator(from: time,
//                               to: time + dt,
//                               first_try_of_stepsize: dt,
//                               smallest_allowed_value_of_stepsize: 1.0e-8,
//                               accuracy: 10e-6,
//                               y: &Psi,
//                               derivative_function: schrodingerEquation)
//
//
//
//        time += dt
        
        adaptiveRungeKuttaOverRange(f: schrodingerEquation, y: &Psi, from: time, to: time + dt, h: dt, relativeTol: 1e-6)
        
        time += dt
        
        /*
            More sophisticated versions of this class might allow one to choose which
            integrator to use as we could replace the above with e.g.
         
                  Psi = doRungeKuttaStep(t: time,
                                         h: dt,
                                         y: Psi,
                                         derivative_function: schrodingerEquation)
        */
    }
}
//  Created by M J Everitt on 21/01/2022.


/*Quick hack. making equivalent structures for each of the types conforming to OperatorType so far. That is, Matrix, SparseMatrix and DiagonalSparseMatrix.
 Note, this is NOT good practice.
*/


open class TimeIndependentSchrodingerSparseSystem {
    public var Psi: Vector
    public var minus_i_H: SparseMatrix
    public var time: Double
    
    public init(initialstate: Vector,
                hamiltonian: SparseMatrix) {
        let minusi = Complex( real: 0.0, imag: -1.0 )
        time = 0.0
        Psi = initialstate
        minus_i_H = minusi * hamiltonian
    }
    
    open func schrodingerEquation(psi: Vector, time: Double)-> Vector {
            
      return minus_i_H * psi
        
      
            
    }

    public func evolve(by dt: Double) {
        
//        multiStepIvpIntegrator(from: time,
//                               to: time + dt,
//                               first_try_of_stepsize: dt,
//                               smallest_allowed_value_of_stepsize: 1.0e-8,
//                               accuracy: 10e-6,
//                               y: &Psi,
//                               derivative_function: schrodingerEquation)
//
//        time += dt
        
        adaptiveRungeKuttaOverRange(f: schrodingerEquation, y: &Psi, from: time, to: time + dt, h: dt, relativeTol: 1e-6)
        
        time += dt
        
        /*
            More sophisticated versions of this class might allow one to choose which
            integrator to use as we could replace the above with e.g.
         
                  Psi = doRungeKuttaStep(t: time,
                                         h: dt,
                                         y: Psi,
                                         derivative_function: schrodingerEquation)
        */
    }
    
    
    public func simulate(times: [Double], takeExpectationOf qoperators: SparseMatrix...) -> [[Complex]] {
        
        var expectationValues = [[Complex]](repeating: [Complex](repeating: Complex(real: 0), count: times.count), count: qoperators.count)
        
        for i in 0..<times.count-1 {
            
            for j in 0..<qoperators.count {
                expectationValues[j][i] = qoperators[j].expectationValue(of: Psi)
            }
            
            evolve(by: times[i+1]-times[i])
            
        }
        
        for j in 0..<qoperators.count {
            expectationValues[j][times.count-1] = qoperators[j].expectationValue(of: Psi)
            
        }
        
        return expectationValues
        
        
    }
}


open class TimeIndependentSchrodingerDiagonalSparseSystem {
    public var Psi: Vector
    public var minus_i_H: DiagonalSparseMatrix
    public var time: Double
    
    public init(initialstate: Vector,
                hamiltonian: DiagonalSparseMatrix) {
        let minusi = Complex(real: 0.0, imag: -1.0 )
        time = 0.0
        Psi = initialstate
        minus_i_H = minusi * hamiltonian
    }
    
    open func schrodingerEquation(psi: Vector, time: Double)-> Vector {
            
      return minus_i_H * psi
        
      
            
    }

    public func evolve(by dt: Double) {
        
//        multiStepIvpIntegrator(from: time,
//                               to: time + dt,
//                               first_try_of_stepsize: dt,
//                               smallest_allowed_value_of_stepsize: 1.0e-8,
//                               accuracy: 10e-6,
//                               y: &Psi,
//                               derivative_function: schrodingerEquation)
        
        adaptiveRungeKuttaOverRange(f: schrodingerEquation, y: &Psi, from: time, to: time + dt, h: dt, relativeTol: 1e-6)
        
        time += dt
        
        /*
            More sophisticated versions of this class might allow one to choose which
            integrator to use as we could replace the above with e.g.
         
                  Psi = doRungeKuttaStep(t: time,
                                         h: dt,
                                         y: Psi,
                                         derivative_function: schrodingerEquation)
        */
    }
    
    
    public func simulate(times: [Double], takeExpectationOf qoperators: DiagonalSparseMatrix...) -> [[Complex]] {
        
        var expectationValues = [[Complex]](repeating: [Complex](repeating: Complex(real: 0), count: times.count), count: qoperators.count)
        
        for i in 0..<times.count-1 {
            
            for j in 0..<qoperators.count {
                expectationValues[j][i] = qoperators[j].expectationValue(of: Psi)
            }
            
            evolve(by: times[i+1]-times[i])
            
        }
        
        for j in 0..<qoperators.count {
            expectationValues[j][times.count-1] = qoperators[j].expectationValue(of: Psi)
            
        }
        
        return expectationValues
        
        
    }
}

