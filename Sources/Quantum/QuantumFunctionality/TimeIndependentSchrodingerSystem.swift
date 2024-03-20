// An example of a class that can solve the time dependent Schrödinger equation for a time independent Hamiltonian (for sparse or dense cases).

import Foundation
open class TimeIndependentSchrodingerSystem {
    public var Psi: StateVector
    public var minus_i_H: MatrixOperator
    public var time: Real
    
    public init(initialstate: StateVector,
                hamiltonian: MatrixOperator) {
        let minusi = ComplexReal( real: 0.0, imag: -1.0 )
        time = 0.0
        Psi = initialstate
        minus_i_H = minusi * hamiltonian
    }
    
    open func schrodingerEquation(time: Real, psi: StateVector)-> StateVector {
            
        if diagonalSparse { return diagSparse_minus_i_H! * psi }
        
        if sparse { return sparse_minus_i_H! * psi }
        
        return minus_i_H * psi
            
    }
    
    
    public var diagonalSparse = false
    public var diagSparse_minus_i_H: DiagonalSparseMatrix<ComplexReal>?
    
    
    public var sparse = false
    public var sparse_minus_i_H: SparseMatrix<ComplexReal>?
    
    
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

    public func evolve(by dt: Real) {
        
        multiStepIvpIntegrator(from: time,
                               to: time + dt,
                               first_try_of_stepsize: dt,
                               smallest_allowed_value_of_stepsize: 1.0e-8,
                               accuracy: 10e-6,
                               y: &Psi,
                               derivative_function: schrodingerEquation)
        
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
    public var Psi: StateVector
    public var minus_i_H: SparseMatrix<ComplexReal>
    public var time: Real
    
    public init(initialstate: StateVector,
                hamiltonian: SparseMatrix<ComplexReal>) {
        let minusi = ComplexReal( real: 0.0, imag: -1.0 )
        time = 0.0
        Psi = initialstate
        minus_i_H = minusi * hamiltonian
    }
    
    open func schrodingerEquation(time: Real, psi: StateVector)-> StateVector {
            
      return minus_i_H * psi
        
      
            
    }

    public func evolve(by dt: Real) {
        
        multiStepIvpIntegrator(from: time,
                               to: time + dt,
                               first_try_of_stepsize: dt,
                               smallest_allowed_value_of_stepsize: 1.0e-8,
                               accuracy: 10e-6,
                               y: &Psi,
                               derivative_function: schrodingerEquation)
        
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
    
    
    public func simulate(times: [Real], takeExpectationOf qoperators: SparseMatrix<ComplexReal>...) -> [[ComplexReal]] {
        
        var expectationValues = [[ComplexReal]](repeating: [ComplexReal](repeating: ComplexReal(0), count: times.count), count: qoperators.count)
        
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
    public var Psi: StateVector
    public var minus_i_H: DiagonalSparseMatrix<ComplexReal>
    public var time: Real
    
    public init(initialstate: StateVector,
                hamiltonian: DiagonalSparseMatrix<ComplexReal>) {
        let minusi = ComplexReal( real: 0.0, imag: -1.0 )
        time = 0.0
        Psi = initialstate
        minus_i_H = minusi * hamiltonian
    }
    
    open func schrodingerEquation(time: Real, psi: StateVector)-> StateVector {
            
      return minus_i_H * psi
        
      
            
    }

    public func evolve(by dt: Real) {
        
        multiStepIvpIntegrator(from: time,
                               to: time + dt,
                               first_try_of_stepsize: dt,
                               smallest_allowed_value_of_stepsize: 1.0e-8,
                               accuracy: 10e-6,
                               y: &Psi,
                               derivative_function: schrodingerEquation)
        
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
    
    
    public func simulate(times: [Real], takeExpectationOf qoperators: DiagonalSparseMatrix<ComplexReal>...) -> [[ComplexReal]] {
        
        var expectationValues = [[ComplexReal]](repeating: [ComplexReal](repeating: ComplexReal(0), count: times.count), count: qoperators.count)
        
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

