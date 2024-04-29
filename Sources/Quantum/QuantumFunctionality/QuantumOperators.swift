/*
 Contains code for generating the operators needed to test and produce production
 figures for Jaynes Cummings model.
 
 In a production library would expand to contain a catalogue of standard states.
 
 In a violation of TDD some other operators such as Parity and Angular momentum
 were added as it seemed sensible to do this while in the correct frame of mind.
 
 There is some duplication of capability as general angular momentum operators can
 also make Pauli operators. Is this an allowable exception or should we remove
 this from future versions of the code.
*/

import Foundation

public class StateSpace: VectorSpace {
    public var numberOperator: Matrix { return makeNumberOperator() }
    public var nullOperator: Matrix { return makeNullOperator() }
    private var I = Complex(real: 0, imag: 1)

    private func makeNullOperator() -> Matrix {
        var output = Matrix(in: self)
        for i in 0 ..< dimension {
            output[i, i] = Complex(real: 0.0)
        }
        return output
    }

    private func makeNumberOperator() -> Matrix {
        var output = Matrix(in: self)
        for i in 0 ..< dimension {
            output[i, i] = Complex(real: Double(i))
        }
        return output
    }
    public var annihilationOperator: Matrix {
        return makeAnnihilationOperator()
    }

    private func makeAnnihilationOperator() -> Matrix {
        var output = Matrix(in: self)
        for i in 0 ..< output.space.dimension - 1 {
            output[i, i+1] = Complex(real: sqrt( Double(i+1) ) , imag: 0.0)
        }
        return output
    }
    
    public var creationOperator: Matrix {
        return makeCreationOperator()
    }
    
    private func makeCreationOperator() -> Matrix {
        var output = Matrix(in: self)
        for i in 0 ..< output.space.dimension - 1 {
            output[i+1, i] = Complex(real: sqrt( Double(i+1) ) , imag: 0.0)
        }
        return output
    }
    
    public var parityOperator: Matrix { return makeParityOperator() }
    private func makeParityOperator() -> Matrix {
        var output = Matrix(in: self)
        var alternatingSignOne = 1 // -1^0
        for n in 0 ..< output.space.dimension {
            output[n, n] = Complex(real: Double(alternatingSignOne)) // -1^n
            alternatingSignOne = -alternatingSignOne
        }
        return output
    }
    // MARK: - Pauli operators
    // Use these so often its worth having them as a special case
    public var sigmaX: Matrix { return makePauliSpinX() }
    public var sigmaY: Matrix { return makePauliSpinY() }
    public var sigmaZ: Matrix { return makePauliSpinZ() }
    public var sigmaPlus: Matrix { return makeSpinRaising() }
    public var sigmaMinus: Matrix { return makeSpinLowering() }

    private func makePauliSpinX () -> Matrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
        var output = Matrix(in: self)

        output[0, 1] = Complex(real: 1)
        output[1, 0] = Complex(real: 1)
        return output
    }
    private func makePauliSpinY () -> Matrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
    
        var output = Matrix(in: self)

        output[0, 1] = -I
        output[1, 0] = I
        return output
    }
    private func makePauliSpinZ () -> Matrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")

        var output = Matrix(in: self)

        output[0, 0] = Complex(real: 1)
        output[1, 1] = Complex(real: -1)
        return output
    }
    private func makeSpinRaising () -> Matrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")

        var output = Matrix(in: self)

        output[0, 1] = Complex(real: 1)
        return output
    }
    private func makeSpinLowering () -> Matrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")

        var output = Matrix(in: self)

        output[1, 0] = Complex(real: 1)
        return output
    }
// MARK: - Now gneral Anular monetum opertaors in the standard baisis
    public var Jx: Matrix { return makeAnguarMonetumInStandrardBasisJx() }
    public var Jy: Matrix { return makeAnguarMonetumInStandrardBasisJy() }
    public var Jz: Matrix { return makeAnguarMonetumInStandrardBasisJz() }
    public var Jplus: Matrix{ return makeAnguarMonetumInStandrardBasisJplus() }
    public var Jminus: Matrix { return makeAnguarMonetumInStandrardBasisJminus() }
    
    
    public func makeAnguarMonetumInStandrardBasisJx () -> Matrix {
        var output = Matrix(in: self)

        let s = Double(dimension - 1)/2.0
        
        for i in 0 ..< dimension {
            for j in 0 ..< dimension {
                let m = s - Double(i)
                let mp = s - Double(j)
                let temp =  Double( delta(i, j+1) + delta(i+1,j) ) * 0.5 * Double.sqrt( (s * (s + 1)) - (m * mp) )
                output[i,j] = Complex(real: temp)
            }
        }
        return output
    }
    public func makeAnguarMonetumInStandrardBasisJy () -> Matrix {
        var output = Matrix(in: self)

        let s = Double(dimension - 1)/2.0

        for i in 0 ..< dimension {
            for j in 0 ..< dimension {
                let m = s - Double(i)
                let mp = s - Double(j)
                // standard basis means i,j index the oposite way to m and m'.
                let temp =  Double(delta(i+1,j) - delta(i, j+1)) * 0.5 * Double.sqrt( (s * (s + 1)) - Double(m * mp) )
                output[j,i] = Complex(real: temp) * I
            }
        }
        return output
    }
    public func makeAnguarMonetumInStandrardBasisJz () -> Matrix {
        var output = Matrix(in: self)

        let s = Double(dimension - 1)/2.0

        for i in 0 ..< dimension {
                let m = s - Double(i)
                output[i,i] = Complex(real: m)
        }
        return output
    }
    public func makeAnguarMonetumInStandrardBasisJplus () -> Matrix {
        var output = Matrix(in: self)

        let s = Double(dimension - 1)/2.0

        for i in 0 ..< dimension {
            for j in 0 ..< dimension {
                let m = s - Double(i)
                let mp = s - Double(j)
                // standard basis means i,j index the oposite way to m and m'.
                let temp =  Double(delta(i+1, j)) *  Double.sqrt( (s * (s + 1)) - Double(m * mp) )
                output[i,j] = Complex(real: temp)
            }
        }
        return output
    }
    public func makeAnguarMonetumInStandrardBasisJminus () -> Matrix {
        var output = Matrix(in: self)
        // TODO: Double does fix precision here so could be improved.
        let s = Double(dimension - 1)/2.0

        for i in 0 ..< dimension {
            for j in 0 ..< dimension {
                let m = s - Double(i)
                let mp = s - Double(j)
                // standard basis means i,j index the opposite way to m and m'.
                let temp =  Double(delta(i, j+1)) *  Double.sqrt( (s * (s + 1)) - Double(m * mp) )
                output[i,j] = Complex(real: temp)
            }
        }
        return output
    }

    public func exponentialOfScaledAnnihilationOperator(scaleFactor alpha: Complex) -> Matrix {
        var output = Matrix(in: self)
        // TODO: Double does fix precision here so could be improved.
        for m in 0 ..< dimension {
            output[m,m] = Complex(real: 1.0)
            for n in m + 1 ..< dimension {
                output[m,n] = output[m,n-1] * alpha * Complex(real: Double.sqrt(Double(n))/Double(n-m))
            }
        }
        return output
    }
    
    
    public func exponentialOfScaledCreationOperator(scaleFactor alpha: Complex) -> Matrix {
        exponentialOfScaledAnnihilationOperator(scaleFactor: alpha).transpose()
    }
    
    public func pauliVectorOperatorExponentialEulerFormula(x_Component ax: Double,
                                                          y_Component ay: Double,
                                                          z_Component az: Double ) -> Matrix {
        let theta =  sqrt(ax * ax + ay * ay + az * az)
        let normalised_x = ax / theta
        let normalised_y = ay / theta
        let normalised_z = az / theta
        
        let na = normalised_x * self.sigmaX + normalised_y * self.sigmaY + normalised_z * self.sigmaZ
        
        return self.identityOperator * cos(theta) + I * na * sin(theta)
    }
    
    // might use these a lot - so worth the duplication even though not DRY
    public func exponentialEulerFormulaSigmaX(_ theta: Double) -> Matrix {
        return self.identityOperator * cos(theta) + I * self.sigmaX * sin(theta)
    }
    
    public func exponentialEulerFormulaSigmaY(_ theta: Double) -> Matrix {
       
        return self.identityOperator * cos(theta) + I * self.sigmaY * sin(theta)
    }
    
    public func exponentialEulerFormulaSigmaZ(_ theta: Double) -> Matrix {
       
        return self.identityOperator * cos(theta) + I * self.sigmaZ * sin(theta)
    }

    
    // see problem 2.1 Nielsen & Chuang, Quantum Computation and Quantum Information 10th anniversary edition
    public func functionOfPauliVectorOperator(x_Component ax: Double,
                                              y_Component ay: Double,
                                              z_Component az: Double,
                                              function f: (Double) -> Complex ) -> Matrix {
        let theta =  sqrt(ax * ax + ay * ay + az * az)
        let normalised_x = ax / theta
        let normalised_y = ay / theta
        let normalised_z = az / theta
        
        let na = normalised_x * self.sigmaX + normalised_y * self.sigmaY + normalised_z * self.sigmaZ
        
        let half = 0.5
        
        return half * (self.identityOperator * ( f(theta) + f(-theta) ) + na * ( f(theta) - f(-theta) ) )
    }
    
    
    public var identityOperator_diagonalSparse: DiagonalSparseMatrix {return makeIdentityMatrixDiagonalSparse()}
    public var SigmaZ_diagonalSparse: DiagonalSparseMatrix {return makeSigmaZDiagonalSparse()}
    public var sigmaPlus_diagonalSparse: DiagonalSparseMatrix { return makeSpinRaisingDiagonalSparse() }
    public var sigmaMinus_diagonalSparse: DiagonalSparseMatrix { return makeSpinLoweringDiagonalSparse() }
    public var numberOperator_diagonalSparse: DiagonalSparseMatrix { return makeNumberOperatorDiagonalSparse() }
    public var annhilationOperator_diagonalSparse: DiagonalSparseMatrix { return makeAnnhilationOperatorDiagonalSparse() }
    public var creationOperator_diagonalSparse: DiagonalSparseMatrix { return makeCreationOperatorDiagonalSparse() }
    
    
    public func makeIdentityMatrixDiagonalSparse() -> DiagonalSparseMatrix {
        var mainDiag = MatrixDiagonal(dimension: dimension, diagIdx: 0, elements: [:])
        
        for i in 0..<dimension {
            mainDiag[i] = Complex(real: 1.0)
        }
        
        return DiagonalSparseMatrix(in: self, diagonals: [0: mainDiag])
    }
    
    public func makeSigmaZDiagonalSparse() -> DiagonalSparseMatrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
            
        let Sz = DiagonalSparseMatrix(in: self,
                                                   diagonals: [0: MatrixDiagonal(dimension: 2,
                                                                                            diagIdx: 0,
                                                                                              elements: [0:Complex(real: 1.0),
                                                                                                         1:Complex(real: -1.0)])])
        
        return Sz
    }
    
    public func makeSpinRaisingDiagonalSparse() -> DiagonalSparseMatrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
        
        let S_plus = DiagonalSparseMatrix(in: self,
                                                       diagonals: [1:
                                                                    MatrixDiagonal(dimension: 2, diagIdx: 1,
                                                                                   elements: [0:Complex(real: 1.0)])])
        
        return S_plus
    }
    
    public func makeSpinLoweringDiagonalSparse() -> DiagonalSparseMatrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
        
        let S_minus = DiagonalSparseMatrix(in: self, diagonals: [-1 : MatrixDiagonal(dimension: 2, diagIdx: -1, elements: [1:Complex(real: 1.0)])])
        
        return S_minus
    }
    
    
    public func makeNumberOperatorDiagonalSparse() -> DiagonalSparseMatrix {
        
        var mainDiag = MatrixDiagonal(dimension: dimension, diagIdx: 0, elements: [:])
        
        for i in 1..<dimension {
            mainDiag[i] = Complex(real: Double(i))
        }
        
        return DiagonalSparseMatrix(in: self, diagonals: [0: mainDiag])
    }
    
    public func makeAnnhilationOperatorDiagonalSparse() -> DiagonalSparseMatrix {
        
        var oneDiag = MatrixDiagonal(dimension: dimension, diagIdx: 1, elements: [:])
        
        for i in 0..<dimension-1 {
            oneDiag[i] = Complex(real: sqrt(Double(i)))
        }
        
        return DiagonalSparseMatrix(in: self, diagonals: [1: oneDiag])
    }
    
    
    public func makeCreationOperatorDiagonalSparse() -> DiagonalSparseMatrix {
        
        var minusOneDiag = MatrixDiagonal(dimension: dimension, diagIdx: -1, elements: [:])
        
        for i in 1..<dimension {
            minusOneDiag[i] = Complex(real: sqrt(Double(i)))
        }
        
        return DiagonalSparseMatrix(in: self, diagonals: [-1: minusOneDiag])
    }
    
    
    public var identityOperator_sparse: SparseMatrix {return makeIdentityMatrixSparse()}
    public var SigmaZ_sparse: SparseMatrix {return makeSigmaZSparse()}
    public var sigmaPlus_sparse: SparseMatrix { return makeSpinRaisingSparse() }
    public var sigmaMinus_sparse: SparseMatrix { return makeSpinLoweringSparse() }
    public var numberOperator_sparse: SparseMatrix { return makeNumberOperatorSparse() }
    public var annhilationOperator_sparse: SparseMatrix { return makeAnnhilationOperatorSparse() }
    public var creationOperator_sparse: SparseMatrix { return makeCreationOperatorSparse() }
    
    
    public func makeIdentityMatrixSparse() -> SparseMatrix {
        var output = SparseMatrix(in: self)
        for row in 0..<dimension {
            output.values.append(CoordinateStorage(value: Complex(real: 1.0), row: row, col: row))
            output.nonzero_elements_per_row[row] = 1
            output.row_first_element_offsets[row] = row
        }
        return output
    }
    
    public func makeSigmaZSparse() -> SparseMatrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
            
        var Sz = SparseMatrix(values: [CoordinateStorage(value: Complex(real: 1.0), row: 0, col: 0),
                                       CoordinateStorage(value: Complex(real: -1.0), row: 1, col: 1)],
                                           in: self)

        
        return Sz
    }
    
    public func makeSpinRaisingSparse() -> SparseMatrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
        
        var S_plus = SparseMatrix(values: [CoordinateStorage(value: Complex(real: 1.0), row: 0, col: 1)], in: self)

        return S_plus
    }
    
    public func makeSpinLoweringSparse() -> SparseMatrix {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
        
        var S_minus = SparseMatrix(values: [CoordinateStorage(value: Complex(real: 1.0), row: 1, col: 0)], in: self )
        return S_minus
    }
    
    
    public func makeNumberOperatorSparse() -> SparseMatrix {
        
        var output = SparseMatrix(in: self)
        for i in 1..<dimension {
            output.values.append(CoordinateStorage(value: Complex(real: Double(i)), row: i, col: i))
            output.nonzero_elements_per_row[i] = 1
            output.row_first_element_offsets[i] = i-1
        }
        return output
    }
    
    public func makeAnnhilationOperatorSparse() -> SparseMatrix {
        
        var output = SparseMatrix(in: self)
        
        for i in 0..<dimension-1 {
            output.values.append(CoordinateStorage(value: Complex(real: sqrt(Double(i))), row: i, col: i+1))
            output.nonzero_elements_per_row[i] = 1
            output.row_first_element_offsets[i] = i
        }
        
        return output
    }
    
    
    public func makeCreationOperatorSparse() -> SparseMatrix {
        
        var output = SparseMatrix(in: self)
        
        for i in 1..<dimension {
            output.values.append(CoordinateStorage(value: Complex(real: sqrt(Double(i))), row: i, col: i-1))
            output.nonzero_elements_per_row[i] = 1
            output.row_first_element_offsets[i] = i-1
        }
        
        return output
    }
    
    

}
//  Created by M J Everitt on 21/01/2022.
