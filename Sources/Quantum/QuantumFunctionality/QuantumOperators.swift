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

public class StateSpace: VectorSpace<Complex<Real>> {
    public var numberOperator: MatrixOperator { return makeNumberOperator() }
    public var nullOperator: MatrixOperator { return makeNullOperator() }
    private var I = ScalarField(real: 0, imag: 1)

    private func makeNullOperator() -> MatrixOperator {
        var output = MatrixOperator(in: self)
        for i in 0 ..< dimension {
            output[i, i] = ScalarField(0)
        }
        return output
    }

    private func makeNumberOperator() -> MatrixOperator {
        var output = MatrixOperator(in: self)
        for i in 0 ..< dimension {
            output[i, i] = ScalarField(i)
        }
        return output
    }
    public var annihilationOperator: MatrixOperator {
        return makeAnnihilationOperator()
    }

    private func makeAnnihilationOperator() -> MatrixOperator {
        var output = MatrixOperator(in: self)
        for i in 0 ..< output.space.dimension - 1 {
            output[i, i+1] = ScalarField( real: Real.sqrt( Real(i+1) ) , imag: Real(0.0) )
        }
        return output
    }
    
    public var creationOperator: MatrixOperator {
        return makeCreationOperator()
    }
    
    private func makeCreationOperator() -> MatrixOperator {
        var output = MatrixOperator(in: self)
        for i in 0 ..< output.space.dimension - 1 {
            output[i+1, i] = ScalarField( real: Real.sqrt( Real(i+1) ) , imag: Real(0.0) )
        }
        return output
    }
    
    public var parityOperator: MatrixOperator { return makeParityOperator() }
    private func makeParityOperator() -> MatrixOperator {
        var output = MatrixOperator(in: self)
        var alternatingSignOne = 1 // -1^0
        for n in 0 ..< output.space.dimension {
            output[n, n] = ScalarField(alternatingSignOne) // -1^n
            alternatingSignOne = -alternatingSignOne
        }
        return output
    }
    // MARK: - Pauli operators
    // Use these so often its worth having them as a special case
    public var sigmaX: MatrixOperator { return makePauliSpinX() }
    public var sigmaY: MatrixOperator { return makePauliSpinY() }
    public var sigmaZ: MatrixOperator { return makePauliSpinZ() }
    public var sigmaPlus: MatrixOperator { return makeSpinRaising() }
    public var sigmaMinus: MatrixOperator { return makeSpinLowering() }

    private func makePauliSpinX () -> MatrixOperator {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
        var output = MatrixOperator(in: self)

        output[0, 1] = ScalarField(1)
        output[1, 0] = ScalarField(1)
        return output
    }
    private func makePauliSpinY () -> MatrixOperator {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
    
        var output = MatrixOperator(in: self)

        output[0, 1] = -I
        output[1, 0] = I
        return output
    }
    private func makePauliSpinZ () -> MatrixOperator {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")

        var output = MatrixOperator(in: self)

        output[0, 0] = ScalarField(1)
        output[1, 1] = ScalarField(-1)
        return output
    }
    private func makeSpinRaising () -> MatrixOperator {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")

        var output = MatrixOperator(in: self)

        output[0, 1] = ScalarField(1)
        return output
    }
    private func makeSpinLowering () -> MatrixOperator {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")

        var output = MatrixOperator(in: self)

        output[1, 0] = ScalarField(1)
        return output
    }
// MARK: - Now gneral Anular monetum opertaors in the standard baisis
    public var Jx: MatrixOperator { return makeAnguarMonetumInStandrardBasisJx() }
    public var Jy: MatrixOperator { return makeAnguarMonetumInStandrardBasisJy() }
    public var Jz: MatrixOperator { return makeAnguarMonetumInStandrardBasisJz() }
    public var Jplus: MatrixOperator { return makeAnguarMonetumInStandrardBasisJplus() }
    public var Jminus: MatrixOperator { return makeAnguarMonetumInStandrardBasisJminus() }
    
    
    public func makeAnguarMonetumInStandrardBasisJx () -> MatrixOperator {
        var output = MatrixOperator(in: self)

        let s = Double(dimension - 1)/2.0
        
        for i in 0 ..< dimension {
            for j in 0 ..< dimension {
                let m = s - Double(i)
                let mp = s - Double(j)
                let temp =  Double( delta(i, j+1) + delta(i+1,j) ) * 0.5 * Double.sqrt( (s * (s + 1)) - (m * mp) )
                output[i,j] = ScalarField(temp)
            }
        }
        return output
    }
    public func makeAnguarMonetumInStandrardBasisJy () -> MatrixOperator {
        var output = MatrixOperator(in: self)

        let s = Double(dimension - 1)/2.0

        for i in 0 ..< dimension {
            for j in 0 ..< dimension {
                let m = s - Double(i)
                let mp = s - Double(j)
                // standard basis means i,j index the oposite way to m and m'.
                let temp =  Double(delta(i+1,j) - delta(i, j+1)) * 0.5 * Double.sqrt( (s * (s + 1)) - Double(m * mp) )
                output[j,i] = ScalarField(temp) * I
            }
        }
        return output
    }
    public func makeAnguarMonetumInStandrardBasisJz () -> MatrixOperator {
        var output = MatrixOperator(in: self)

        let s = Double(dimension - 1)/2.0

        for i in 0 ..< dimension {
                let m = s - Double(i)
                output[i,i] = ScalarField(m)
        }
        return output
    }
    public func makeAnguarMonetumInStandrardBasisJplus () -> MatrixOperator {
        var output = MatrixOperator(in: self)

        let s = Double(dimension - 1)/2.0

        for i in 0 ..< dimension {
            for j in 0 ..< dimension {
                let m = s - Double(i)
                let mp = s - Double(j)
                // standard basis means i,j index the oposite way to m and m'.
                let temp =  Double(delta(i+1, j)) *  Double.sqrt( (s * (s + 1)) - Double(m * mp) )
                output[i,j] = ScalarField(temp)
            }
        }
        return output
    }
    public func makeAnguarMonetumInStandrardBasisJminus () -> MatrixOperator {
        var output = MatrixOperator(in: self)
        // TODO: Double does fix precision here so could be improved.
        let s = Double(dimension - 1)/2.0

        for i in 0 ..< dimension {
            for j in 0 ..< dimension {
                let m = s - Double(i)
                let mp = s - Double(j)
                // standard basis means i,j index the opposite way to m and m'.
                let temp =  Double(delta(i, j+1)) *  Double.sqrt( (s * (s + 1)) - Double(m * mp) )
                output[i,j] = ScalarField(temp)
            }
        }
        return output
    }

    public func exponentialOfScaledAnnihilationOperator(scaleFactor alpha: ScalarField) -> MatrixOperator {
        var output = MatrixOperator(in: self)
        // TODO: Double does fix precision here so could be improved.
        for m in 0 ..< dimension {
            output[m,m] = ScalarField(real: 1.0)
            for n in m + 1 ..< dimension {
                output[m,n] = output[m,n-1] * alpha * ScalarField(Double.sqrt(Double(n))/Double(n-m))
            }
        }
        return output
    }
    public func exponentialOfScaledCreationOperator(scaleFactor alpha: ScalarField) -> MatrixOperator {
        exponentialOfScaledAnnihilationOperator(scaleFactor: alpha).transpose()
    }
    
    public func pauliVectorOperatorExponentialEulerFormula(x_Component ax: Real,
                                                          y_Component ay: Real,
                                                          z_Component az: Real ) -> MatrixOperator {
        let theta =  Real.sqrt(ax * ax + ay * ay + az * az)
        let normalised_x = ax / theta
        let normalised_y = ay / theta
        let normalised_z = az / theta
        
        let na = normalised_x * self.sigmaX + normalised_y * self.sigmaY + normalised_z * self.sigmaZ
        
        let i = ScalarField( real: Real(0.0), imag: Real(1.0) )
        
        return self.identityOperator * Real.cos(theta) + i * na * Real.sin(theta)
    }
    
    // might use these a lot - so worth the duplication even though not DRY
    public func exponentialEulerFormulaSigmaX(_ theta: Real) -> MatrixOperator {
        let i = ScalarField( real: Real(0.0), imag: Real(1.0) )
        return self.identityOperator * Real.cos(theta) + i * self.sigmaX * Real.sin(theta)
    }
    
    public func exponentialEulerFormulaSigmaY(_ theta: Real) -> MatrixOperator {
        let i = ScalarField(real: Real(0.0), imag: Real(1.0) )
        return self.identityOperator * Real.cos(theta) + i * self.sigmaY * Real.sin(theta)
    }
    
    public func exponentialEulerFormulaSigmaZ(_ theta: Real) -> MatrixOperator {
        let i = ScalarField( real: Real(0.0), imag: Real(1.0) )
        return self.identityOperator * Real.cos(theta) + i * self.sigmaZ * Real.sin(theta)
    }

    
    // see problem 2.1 Nielsen & Chuang, Quantum Computation and Quantum Information 10th anniversary edition
    public func functionOfPauliVectorOperator(x_Component ax: Real,
                                              y_Component ay: Real,
                                              z_Component az: Real,
                                              function f: (Real) -> Complex<Real> ) -> MatrixOperator {
        let theta =  Real.sqrt(ax * ax + ay * ay + az * az)
        let normalised_x = ax / theta
        let normalised_y = ay / theta
        let normalised_z = az / theta
        
        let na = normalised_x * self.sigmaX + normalised_y * self.sigmaY + normalised_z * self.sigmaZ
        
        let half = Real(1)/Real(2)
        
        return half * (self.identityOperator * ( f(theta) + f(-theta) ) + na * ( f(theta) - f(-theta) ) )
    }
    
    
    public var identityOperator_diagonalSparse: DiagonalSparseMatrix<ComplexReal> {return makeIdentityMatrixDiagonalSparse()}
    public var SigmaZ_diagonalSparse: DiagonalSparseMatrix<ComplexReal> {return makeSigmaZDiagonalSparse()}
    public var sigmaPlus_diagonalSparse: DiagonalSparseMatrix<ComplexReal> { return makeSpinRaisingDiagonalSparse() }
    public var sigmaMinus_diagonalSparse: DiagonalSparseMatrix<ComplexReal> { return makeSpinLoweringDiagonalSparse() }
    public var numberOperator_diagonalSparse: DiagonalSparseMatrix<ComplexReal> { return makeNumberOperatorDiagonalSparse() }
    public var annhilationOperator_diagonalSparse: DiagonalSparseMatrix<ComplexReal> { return makeAnnhilationOperatorDiagonalSparse() }
    public var creationOperator_diagonalSparse: DiagonalSparseMatrix<ComplexReal> { return makeCreationOperatorDiagonalSparse() }
    
    
    public func makeIdentityMatrixDiagonalSparse() -> DiagonalSparseMatrix<ComplexReal> {
        var mainDiag = MatrixDiagonal<ComplexReal>(dimension: dimension, diagIdx: 0, elements: [:])
        
        for i in 0..<dimension {
            mainDiag[i] = ComplexReal(real: 1)
        }
        
        return DiagonalSparseMatrix(in: self, diagonals: [0: mainDiag])
    }
    
    public func makeSigmaZDiagonalSparse() -> DiagonalSparseMatrix<ComplexReal> {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
            
        let Sz = DiagonalSparseMatrix<ComplexReal>(in: self,
                                                   diagonals: [0: MatrixDiagonal<ComplexReal>(dimension: 2,
                                                                                            diagIdx: 0,
                                                                                            elements: [0:ComplexReal(real: 1),
                                                                                                       1:ComplexReal(real: -1)])])
        
        return Sz
    }
    
    public func makeSpinRaisingDiagonalSparse() -> DiagonalSparseMatrix<ComplexReal> {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
        
        let S_plus = DiagonalSparseMatrix<ComplexReal>(in: self,
                                                       diagonals: [1:
                                                                    MatrixDiagonal<ComplexReal>(dimension: 2, diagIdx: 1,
                                                                                                elements: [0:ComplexReal(real: 1)])])
        
        return S_plus
    }
    
    public func makeSpinLoweringDiagonalSparse() -> DiagonalSparseMatrix<ComplexReal> {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
        
        let S_minus = DiagonalSparseMatrix<ComplexReal>(in: self, diagonals: [-1 : MatrixDiagonal<ComplexReal>(dimension: 2, diagIdx: -1, elements: [1:ComplexReal(real: 1)])])
        
        return S_minus
    }
    
    
    public func makeNumberOperatorDiagonalSparse() -> DiagonalSparseMatrix<ComplexReal> {
        
        var mainDiag = MatrixDiagonal<ComplexReal>(dimension: dimension, diagIdx: 0, elements: [:])
        
        for i in 1..<dimension {
            mainDiag[i] = ComplexReal(real: Real(i))
        }
        
        return DiagonalSparseMatrix<ComplexReal>(in: self, diagonals: [0: mainDiag])
    }
    
    public func makeAnnhilationOperatorDiagonalSparse() -> DiagonalSparseMatrix<ComplexReal> {
        
        var oneDiag = MatrixDiagonal<ComplexReal>(dimension: dimension, diagIdx: 1, elements: [:])
        
        for i in 0..<dimension-1 {
            oneDiag[i] = ComplexReal(real: Real.sqrt(Real(i)))
        }
        
        return DiagonalSparseMatrix<ComplexReal>(in: self, diagonals: [1: oneDiag])
    }
    
    
    public func makeCreationOperatorDiagonalSparse() -> DiagonalSparseMatrix<ComplexReal> {
        
        var minusOneDiag = MatrixDiagonal<ComplexReal>(dimension: dimension, diagIdx: -1, elements: [:])
        
        for i in 1..<dimension {
            minusOneDiag[i] = ComplexReal(real: Real.sqrt(Real(i)))
        }
        
        return DiagonalSparseMatrix<ComplexReal>(in: self, diagonals: [-1: minusOneDiag])
    }
    
    
    public var identityOperator_sparse: SparseMatrix<ComplexReal> {return makeIdentityMatrixSparse()}
    public var SigmaZ_sparse: SparseMatrix<ComplexReal> {return makeSigmaZSparse()}
    public var sigmaPlus_sparse: SparseMatrix<ComplexReal> { return makeSpinRaisingSparse() }
    public var sigmaMinus_sparse: SparseMatrix<ComplexReal> { return makeSpinLoweringSparse() }
    public var numberOperator_sparse: SparseMatrix<ComplexReal> { return makeNumberOperatorSparse() }
    public var annhilationOperator_sparse: SparseMatrix<ComplexReal> { return makeAnnhilationOperatorSparse() }
    public var creationOperator_sparse: SparseMatrix<ComplexReal> { return makeCreationOperatorSparse() }
    
    
    public func makeIdentityMatrixSparse() -> SparseMatrix<ComplexReal> {
        var output = SparseMatrix(in: self)
        for row in 0..<dimension {
            output.values.append(CoordinateStorage(value: ComplexReal(real: 1), row: row, col: row))
            output.nonzero_elements_per_row[row] = 1
            output.row_first_element_offsets[row] = row
        }
        return output
    }
    
    public func makeSigmaZSparse() -> SparseMatrix<ComplexReal> {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
            
        var Sz = SparseMatrix<ComplexReal>(values: [CoordinateStorage(value: ComplexReal(real: 1), row: 0, col: 0),
                                                    CoordinateStorage(value: ComplexReal(real: -1), row: 1, col: 1)],
                                           in: self)

        
        return Sz
    }
    
    public func makeSpinRaisingSparse() -> SparseMatrix<ComplexReal> {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
        
        var S_plus = SparseMatrix(values: [CoordinateStorage(value: ComplexReal(real: 1), row: 0, col: 1)], in: self)

        return S_plus
    }
    
    public func makeSpinLoweringSparse() -> SparseMatrix<ComplexReal> {
        assert(self.dimension == 2, "Dimension should be 2 but is \(self.dimension)")
        
        var S_minus = SparseMatrix(values: [CoordinateStorage(value: ComplexReal(real: 1), row: 1, col: 0)], in: self )
        return S_minus
    }
    
    
    public func makeNumberOperatorSparse() -> SparseMatrix<ComplexReal> {
        
        var output = SparseMatrix(in: self)
        for i in 1..<dimension {
            output.values.append(CoordinateStorage(value: ComplexReal(real: Double(i)), row: i, col: i))
            output.nonzero_elements_per_row[i] = 1
            output.row_first_element_offsets[i] = i-1
        }
        return output
    }
    
    public func makeAnnhilationOperatorSparse() -> SparseMatrix<ComplexReal> {
        
        var output = SparseMatrix(in: self)
        
        for i in 0..<dimension-1 {
            output.values.append(CoordinateStorage(value: ComplexReal(real: Real.sqrt(Real(i))), row: i, col: i+1))
            output.nonzero_elements_per_row[i] = 1
            output.row_first_element_offsets[i] = i
        }
        
        return output
    }
    
    
    public func makeCreationOperatorSparse() -> SparseMatrix<ComplexReal> {
        
        var output = SparseMatrix(in: self)
        
        for i in 1..<dimension {
            output.values.append(CoordinateStorage(value: ComplexReal(real: Real.sqrt(Real(i))), row: i, col: i-1))
            output.nonzero_elements_per_row[i] = 1
            output.row_first_element_offsets[i] = i-1
        }
        
        return output
    }
    
    

}
//  Created by M J Everitt on 21/01/2022.
