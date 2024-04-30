// Ad hoc "Unit Tests" note that this style is not a good example of TDD and a result of
// trying to code a library in only one week (please excuse any unfixed spelling mistakes).

// Too many of the tests focus on implementation and not behaviours.
// One might consider the code that is used to generate the images for each chapter are better tests of the library than those presented here.

import Foundation
import Accelerate
import XCTest
@testable import Quantum

final class CoreTests: XCTestCase {
    typealias C = Complex
    typealias V = Vector
    
    override func setUpWithError() throws {
        // Put setup code here. This method is called before the invocation of each test method in the class.
    }
    
    override func tearDownWithError() throws {
        // Put teardown code here. This method is called after the invocation of each test method in the class.
    }
    
    func test_Complex_Basic() throws {
        
        let i = C(real: 0.0, imag: 1.0)
        let zero = C(0)
        let three = C(3)
        let threeDbl = C(3.0)
        let threeI = C(real: 0.0, imag: 3.0)
        let threeplusthreei = C(real: 3.0, imag: 3.0)
        XCTAssertEqual( three - threeDbl , zero )
        XCTAssertEqual( i * three , threeI )
        XCTAssertEqual( i * i * i, -i )
        XCTAssertEqual( three - threeDbl , zero )
        XCTAssertEqual(threeplusthreei / three , C(real: 1.0, imag: 1.0) )
        XCTAssertEqual(threeplusthreei / 3.0 , C(real: 1.0, imag: 1.0) )
        XCTAssertEqual( 3 + 3.0 * i , threeplusthreei )
        XCTAssertEqual( i * 3.0 + 3.0 , threeplusthreei )
        XCTAssertEqual( -3.0 - 3.0 * i + threeplusthreei, zero )
        XCTAssertEqual( i * -3.0 - 3.0 + threeplusthreei, zero )
        XCTAssertEqual( (3.0 + 4.0 * i).norm, 25.0 )
        XCTAssertEqual( (3.0 + 4.0 * i).modulus, 5.0 )
        
    }
    
    func test_Has_protocolForFunctions() throws {
        let a = C(real: 0.0, imag: 2.0)
        let b = C(modulus: 2.0, argument: Double.pi/2.0)
        XCTAssertEqual( (( a * b ) - C(real: -4.0, imag: 0.0)).modulus , 0.0 , accuracy: 1.0e-14 )
        XCTAssertEqual(a.argument, Double.pi/2.0, accuracy: 1.0e-14)
        XCTAssertEqual((-a).argument, -Double.pi/2.0, accuracy: 1.0e-14)
        
        XCTAssertEqual(b.modulus, 2.0, accuracy: 1.0e-14)
    }
    
    func test_VectorBasics() throws {
        let i = C( real: 0.0, imag: 1.0 )
        
        let space1 = VectorSpace(dimension: 2, label: "test space 1")
        let u1 = V( elements: [ 3.0 + 2.0 * i , -1.0 + 4.0 * i ] , in: space1)
        let v1 = V( elements: [ -3.0 - 2.0 * i , 1.0 - 4.0 * i ] , in: space1)
        let zero1 = V(elements: [ C(real: 0), C(real: 0) ] , in: space1)
        
        let space2 = VectorSpace(dimension: 2, label: "test space 2")
        let u2 = V( elements: [ 3.0 + 2.0 * i , -1.0 + 4.0 * i ] , in: space2)
        
        XCTAssertEqual( u1 + v1 , zero1 )
        XCTAssertTrue( u1.space == v1.space )
        XCTAssertFalse( u1.space == u2.space )
        XCTAssertEqual( (u1 + v1)[1] , zero1[1] )
        XCTAssertNotEqual( u1[1] , v1[1])
        XCTAssertEqual( -u1 , v1 )
        XCTAssertNotEqual( -u1 , -v1 )
        XCTAssertEqual( u1 - v1 , u1 + u1 )
        XCTAssertEqual( u1 - u1 , zero1 )
        XCTAssertNotEqual( u1 - v1 , u1 - u1 )
        
        XCTAssertEqual( sum(u1,v1) , zero1 )
        XCTAssertEqual( sum(u1,v1,u1) , u1 )
        XCTAssertEqual( u1 + u1 , ( C(4.0) * u1 ) / C(2.0) )
        XCTAssertEqual( u1 + u1 , ( u1 * C(4.0) ) / C(2.0) )
        XCTAssertNotEqual( u1 + u1 , ( u1 * C(4.0) ) / C(3.0) )
        XCTAssertEqual( u1 + u1 , ( 4.0 * u1 ) / 2 ) // also sching int and double
        XCTAssertEqual( u1 + u1 , ( u1 * 4.0 ) / 2.0 )
        XCTAssertNotEqual( u1 + u1 , ( u1 * 4.0 ) / 3.0 )
        XCTAssertEqual( u1 + u1 , ( 4.0 * u1 ) / 2.0 )
        XCTAssertEqual( u1 + u1 , ( u1 * 4.0 ) / 2.0 )
        XCTAssertNotEqual( u1 + u1 , ( u1 * 4.0 ) / 3.0 )
        
    }
    
    func test_identity_operator() throws {
        let space = VectorSpace(dimension: 3, label: "test Space")
        
        var elem = [Complex](repeating: Complex(real: 0.0), count: space.dimension * space.dimension)
        
        for i in 0..<space.dimension {
            elem[i * space.dimension + i] = Complex(real: 1.0)
        }
        XCTAssert(
            space.identityOperator.elements == elem
        )
    }
    
    func test_matrixSubscript() throws {
        let space = VectorSpace(dimension: 2, label: "test Space")
        var matrix = Matrix(in: space)
        matrix[0,0] = Complex(1.0)
        matrix[0,1] = Complex(2.0)
        matrix[1,0] = Complex(3.0)
        matrix[1,1] = Complex(4.0)
        
        let elems = [Complex(1.0), Complex(2.0), Complex(3.0), Complex(4.0)]
        XCTAssertTrue( matrix.elements == elems)
    }
    
    func test_matrixTimesVector() throws {
        let space = VectorSpace(dimension: 2, label: "test Space")
        var matrix = Matrix(in: space)
        matrix[0,0] = Complex(1.0)
        matrix[0,1] = Complex(2.0)
        matrix[1,0] = Complex(3.0)
        matrix[1,1] = Complex(4.0)
        var vector = Vector(in: space)
        vector[0] = Complex(5.0)
        vector[1] = Complex(6.0)
        
        XCTAssertTrue( (matrix * vector).elements == [Complex(17.0),Complex(39.0)] )
    }
    
    
    func test_square_matrix_multiplication() throws {
        
        let space = VectorSpace(dimension: 2, label: "MatMult test space")
        
        let A = Matrix(elements: [Complex(1.0), Complex(2.0), Complex(3.0), Complex(4.0)], in: space)
        let B = Matrix(elements: [Complex(5.0), Complex(6.0), Complex(7.0), Complex(8.0)], in: space)
        
        let C = A*B
        
        XCTAssertEqual(C.elements, [Complex(19),Complex(22),Complex(43),Complex(50)])
        
        
        
        
    }
    
    func test_matrix_approximateEquals() throws {
        let i = C( real: 0.0, imag: 1.0 )
        let space = VectorSpace(dimension: 2, label: "test")
        let TINY = 1e-6
        let A = Matrix(elements: [Complex(1.0), Complex(2.0), Complex(3.0), Complex(-4.0)], in: space)
        let B = Matrix(elements: [Complex(1.0),Complex(2.0+0.5 * TINY),Complex(3.0),Complex(-4.0)], in: space)
        let C = Matrix(elements: [Complex(1.0),Complex(2.0 + TINY),Complex(3.0),Complex(-4.0)], in: space)
        XCTAssert( (A =~= B) == true)
        XCTAssert( (A =~= C) == false)
        
        let spaceC = VectorSpace(dimension: 2, label: "test")
        let complexOne = Complex(1)
        let complexZero = Complex(0)
        
        let Ac = Matrix(elements: [i,complexOne,complexZero,complexOne],
                        in: spaceC)
        
        let Bc = Matrix(elements: [i, complexOne + 0.5 * TINY, complexZero, complexOne],
                        in: spaceC)
        
        let Cc = Matrix(elements: [i, complexOne + 2.0 * TINY, complexZero, complexOne],
                        in: spaceC)
        XCTAssert( (Ac =~= Bc) == true)
        XCTAssert( (Ac =~= Cc) == false)
    }
    
    func test_matrix_TransposeAndAdjoint() throws {
        let i = C( real: 0.0, imag: 1.0 )
        let spaceC = VectorSpace(dimension: 2, label: "test")
        let A = Matrix(elements: [1 + i, 2 + 2 * i,
                                  3 - i, i ],
                       in: spaceC)
        let AT = Matrix(elements: [1 + i, 3 - i,
                                   2 + 2 * i, i ],
                        in: spaceC)
        let AD = Matrix(elements: [1 - i, 3 + i,
                                   2 - 2 * i, -i ],
                        in: spaceC)
        XCTAssertEqual(A.transpose(), AT)
        XCTAssertEqual(A.hermitianAdjoint(), AD)
    }
    
    // MARK: - Sparse Tests (adding spartse to make integration faster)
    func test_sparse_creation_and_equals() throws {
        let cs1 = CoordinateStorage(value: 3.0, row: 1, col: 2)
        let cs2 = CoordinateStorage(value: 6.0, row: 1, col: 2)
        XCTAssert(cs1 != cs2)
        let cs3 = CoordinateStorage(value: 3.0, row: 2, col: 2)
        XCTAssert(cs1 != cs3)
        let cs4 = CoordinateStorage(value: 3.0, row: 1, col: 1)
        XCTAssert(cs1 != cs4)
        let cs5 = CoordinateStorage(value: 3.0, row: 1, col: 2)
        XCTAssert(cs1 == cs5)
        XCTAssert( (cs1 * 2.0 == cs2), "not ClosedUnderScalarFieldMultiplication" )
        XCTAssert(cs1 < cs3)
        XCTAssert(cs3 > cs1)
    }
    
    func test_makingSparseMatrix() throws {
        let space = VectorSpace(dimension: 3, label: "test Space")
        let n = Matrix(elements: [0.0, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 2.0], in: space)
        let sparse_n = SparseMatrix(from: n)
        
        XCTAssert(sparse_n.values.count == 2)
        XCTAssert(sparse_n.values[0] == CoordinateStorage(value: 1.0, row: 1, col: 1))
        XCTAssert(sparse_n.values[1] == CoordinateStorage(value: 2.0, row: 2, col: 2))
        XCTAssert(Matrix(fromSparse: sparse_n) == n)
        let a = Matrix(elements: [0.0, 1.0,  0.0,
                                          0.0 , 0.0, sqrt(2.0),
                                          0.0 , 0.0, 0.0], in: space)
        let sparse_a = SparseMatrix(from: a)
        
        let apn = a + n
        
        let sparse_apn = sparse_a + sparse_n
        XCTAssert(Matrix(fromSparse: sparse_apn) == apn)
        
        let v = Vector(elements: [1.0,2.0,3.0], in: space)
        
        XCTAssert( apn * v == sparse_apn * v )
        XCTAssert( Matrix(fromSparse: 2.0 * sparse_apn) == 2.0 * apn)
        XCTAssert( Matrix(fromSparse: sparse_apn * 2.0) == 2.0 * apn)
        XCTAssert( Matrix(fromSparse: sparse_apn * 2.0) == apn * 2.0)
        XCTAssert( Matrix(fromSparse: 2.0 * sparse_apn) == apn * 2.0)
        
        
        let spaceC = StateSpace(dimension: 3, label: "test complex Space")
        let nC = spaceC.numberOperator
        let sparse_nC = SparseMatrix(from: nC)
        let aC = spaceC.annihilationOperator
        let sparse_aC = SparseMatrix(from: aC)
        
        let apnC = aC + nC
        let sparse_apnC = sparse_aC + sparse_nC
        XCTAssert(Matrix(fromSparse: sparse_apnC) == apnC)
        
        // check correct inner product is used if Has_Conjugate
        XCTAssert(v.innerProduct(dualVector: v) == Complex(14))
        let vC = Vector(elements: [1.0,2.0,3.0], in: spaceC)
        XCTAssert(Complex(real: 14.0, imag: 0.0) == vC.innerProduct(dualVector: vC) )
    }
    
    
    func test_tensorProduct() throws {
        let space1 = VectorSpace(dimension: 2, label: "test Space 1")
        let matrix1 = Matrix(elements: [1.0,2.0,3.0,4.0], in: space1)
        let space2 = VectorSpace(dimension: 2, label: "test Space 2")
        let matrix2 = Matrix(elements: [0.0,5.0,6.0,7.0], in: space2)
        let totalSpace = VectorSpace(tensorProductOf: [space1, space2], label: "TP space")
        let bigMatrix = totalSpace.tensorProduct(of: matrix1, with: matrix2)
        
        XCTAssertTrue(bigMatrix.elements == [ 0.0, 5.0 , 0.0  , 10.0 ,
                                              6.0, 7.0 , 12.0 , 14.0,
                                              0.0, 15.0, 0.0  , 20.0,
                                              18.0, 21.0,24.0  , 28.0].map{Complex($0)} )
    }
    
    
    func testParallelWrite() throws {
        var array = Array(repeating: 0, count: 100)

        DispatchQueue.concurrentPerform(iterations: array.count) { index in
            // Write to the array element at index
            array[index] = index
        }

        print(array)

    }
    
    func test_subspaceOfTotalSpaceInit() throws {
        
        let twoQubitLeftSpace = VectorSpace(dimension: 2, label: "two qubit left")
        let twoQubitRightSpace = VectorSpace(dimension: 2, label: "two qubit right")
        
        let twoQubitTotalSpace = VectorSpace(tensorProductOf: [twoQubitLeftSpace, twoQubitRightSpace], label: "")
        
        let twoQubitSubspace1 = twoQubitTotalSpace.createSubspace(spacesToKeep: twoQubitLeftSpace, label: "")
        XCTAssert(twoQubitSubspace1 == twoQubitLeftSpace)
        XCTAssert(twoQubitSubspace1.dimension == 2)
        
        let twoQubitSubspace2 = twoQubitTotalSpace.createSubspace(spacesToKeep: twoQubitRightSpace, label: "")
        print(twoQubitSubspace2.identifier, twoQubitRightSpace.identifier)
        XCTAssert(twoQubitSubspace2 == twoQubitRightSpace)
        XCTAssert(twoQubitSubspace2.dimension == 2)
        
        
        let JCSpinSpace = VectorSpace(dimension: 2, label: "JC Spin Space")
        let JCFieldSpace = VectorSpace(dimension: 30, label: "JC Field Space")
        
        let JCTotalSpace = VectorSpace(tensorProductOf: [JCSpinSpace, JCFieldSpace], label: "JC Total Space")
        
        let JCSubSpace1 = JCTotalSpace.createSubspace(spacesToKeep: JCSpinSpace, label: "")
        XCTAssert(JCSubSpace1 == JCSpinSpace)
        XCTAssert(JCSubSpace1.dimension == 2)
        
        let JCSubSpace2 = JCTotalSpace.createSubspace(spacesToKeep: JCFieldSpace, label: "")
        XCTAssert(JCSubSpace2 == JCFieldSpace)
        XCTAssert(JCSubSpace2.dimension == 30)
        
        let quinpartiteSpace1 = VectorSpace(dimension: 2, label: "Quinpartite Space 1")
        let quinpartiteSpace2 = VectorSpace(dimension: 3, label: "Quinpartite Space 2")
        let quinpartiteSpace3 = VectorSpace(dimension: 4, label: "Quinpartite Space 3")
        let quinpartiteSpace4 = VectorSpace(dimension: 5, label: "Quinpartite Space 4")
        let quinpartiteSpace5 = VectorSpace(dimension: 6, label: "Quinpartite Space 5")
        
        
        let quinpartiteTotalSpace = VectorSpace(tensorProductOf: [quinpartiteSpace1, quinpartiteSpace2, quinpartiteSpace3, quinpartiteSpace4, quinpartiteSpace5], label: "Quinpartite Total Space")
        
        let firstThirdFifthSpaces = quinpartiteTotalSpace.createSubspace(spacesToKeep: quinpartiteSpace1, quinpartiteSpace3, quinpartiteSpace5, label: "First, Third and Fifth Spaces of Quinpartite.")
        
        XCTAssert(firstThirdFifthSpaces.setofSpaces.count == 3)
        XCTAssert(firstThirdFifthSpaces.dimension == 48)
        XCTAssert(firstThirdFifthSpaces.setofSpaces[0].dimension == 2)
        XCTAssert(firstThirdFifthSpaces.setofSpaces[1].dimension == 4)
        XCTAssert(firstThirdFifthSpaces.setofSpaces[2].dimension == 6)
        
        let secondSpace = quinpartiteTotalSpace.createSubspace(spacesToKeep: quinpartiteSpace2, label: "")
        XCTAssert(secondSpace == quinpartiteSpace2)
        XCTAssert(secondSpace.dimension == 3)
        
        let firstAndFifthSpaces = quinpartiteTotalSpace.createSubspace(spacesToKeep: quinpartiteSpace1, quinpartiteSpace5, label: "")
        
        XCTAssert(firstAndFifthSpaces.setofSpaces.count == 2)
        XCTAssert(firstAndFifthSpaces.dimension == 12)
        XCTAssert(firstAndFifthSpaces.setofSpaces[0].dimension == 2)
        XCTAssert(firstAndFifthSpaces.setofSpaces[1].dimension == 6)
        
    }
    
    func test_makingDiagonalSparseMatrix() throws {
        
        let testSpace = VectorSpace(dimension: 4, label: "")
        
        let denseMatrix = Matrix(elements: [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16], in: testSpace)
        
        let diagonalSparse = DiagonalSparseMatrix(from: denseMatrix)
        
        XCTAssertEqual(diagonalSparse.space, testSpace)
        XCTAssertEqual(diagonalSparse.diagonals[-3]!.elements, [3:Complex(13)])
        XCTAssertEqual(diagonalSparse.diagonals[-2]!.elements, [2:Complex(9), 3: Complex(14)])
        XCTAssertEqual(diagonalSparse.diagonals[-1]!.elements, [1:5, 2:10, 3:15].mapValues{Complex($0)})
        XCTAssertEqual(diagonalSparse.diagonals[0]!.elements, [0:1, 1:6, 2:11, 3:16].mapValues{Complex($0)})
        XCTAssertEqual(diagonalSparse.diagonals[1]!.elements, [0:2, 1:7, 2:12].mapValues{Complex($0)})
        XCTAssertEqual(diagonalSparse.diagonals[2]!.elements, [0:3, 1:8].mapValues{Complex($0)})
        XCTAssertEqual(diagonalSparse.diagonals[3]!.elements, [0:4].mapValues{Complex($0)})
        
    }
    
    func test_diagonalSparseAdd() throws {
        
        let space = VectorSpace(dimension: 3, label: "test diag sparse add space")
        let lhs = DiagonalSparseMatrix(in: space, diagonals: [0:MatrixDiagonal(dimension: space.dimension, diagIdx: 0, elements: [0:1, 1:1, 2:1].mapValues{Complex($0)})])
        
        let rhs = DiagonalSparseMatrix(in: space, diagonals: [1:MatrixDiagonal(dimension: space.dimension, diagIdx: 1, elements: [0:1, 1:1].mapValues{Complex($0)})])
       
        XCTAssertEqual((lhs + rhs).diagonals[0]!.elements, [0:1, 1:1, 2:1].mapValues{Complex($0)})
        XCTAssertEqual((lhs + rhs).diagonals[1]!.elements, [0:1, 1:1].mapValues{Complex($0)})
    }
    
    func test_diagonalSparseSubtract() throws {
        
        let space = VectorSpace(dimension: 3, label: "test diag sparse add space")
        let lhs = DiagonalSparseMatrix(in: space, diagonals: [0:MatrixDiagonal(dimension: space.dimension, diagIdx: 0, elements: [0:1, 1:1, 2:1].mapValues{Complex($0)})])
        
        let rhs = DiagonalSparseMatrix(in: space, diagonals: [1:MatrixDiagonal(dimension: space.dimension, diagIdx: 1, elements: [0:1, 1:1].mapValues{Complex($0)})])
        
        XCTAssertEqual((lhs - rhs).diagonals[0]!.elements, [0:1, 1:1, 2:1].mapValues{Complex($0)})
        XCTAssertEqual((lhs - rhs).diagonals[1]!.elements, [0:-1, 1:-1].mapValues{Complex($0)})
    }
    
    func test_diagonalSparseVecMult() throws {
        
        let space = VectorSpace(dimension: 2, label: "test Space")
        
        let zeroDiag = MatrixDiagonal(dimension: 2, diagIdx: 0, elements: [0: 1.0, 1:4.0].mapValues{Complex($0)})
        let oneDiag = MatrixDiagonal(dimension: 2, diagIdx: 1, elements: [0:2.0].mapValues{Complex($0)})
        let minusOneDiag = MatrixDiagonal(dimension: 2, diagIdx: -1, elements: [1:3.0].mapValues{Complex($0)})
        
        let matrix = DiagonalSparseMatrix(in: space,
                                          diagonals: [0: zeroDiag, 1: oneDiag, -1: minusOneDiag])
        
        let vector = Vector(elements: [5.0, 6.0], in: space)
        
        XCTAssertTrue( (matrix * vector).elements == [17.0,39.0].map{Complex($0)} )
    }
    
    func test_diagonalSparseMatMult() throws {
        
        let testDiagMMSpace = VectorSpace(dimension: 4, label: "")
        
        let A = DiagonalSparseMatrix(in: testDiagMMSpace, diagonals:
                                        [0: MatrixDiagonal(dimension: 4, diagIdx: 0, elements: [0:1, 1:3, 2:5, 3:7].mapValues{Complex($0)}),
                                         -1: MatrixDiagonal(dimension: 4, diagIdx: -1, elements: [1: 2, 2:4, 3: 6].mapValues{Complex($0)})])
        
        
        let B = DiagonalSparseMatrix(in: testDiagMMSpace, diagonals:
                                        [-1: MatrixDiagonal(dimension: 4, diagIdx: -1, elements: [1: 1, 2:1, 3:1].mapValues{Complex($0)}),
                                          1: MatrixDiagonal(dimension: 4, diagIdx: 1, elements: [0:1, 1:1, 2:1].mapValues{Complex($0)})])
        
        let C = A*B
        
        XCTAssertEqual(C.diagonals[-2]!.elements, [2:4, 3:6].mapValues{Complex($0)})
        XCTAssertEqual(C.diagonals[-1]!.elements, [1:3, 2:5, 3:7].mapValues{Complex($0)})
        XCTAssertEqual(C.diagonals[0]!.elements, [1:2, 2:4, 3:6].mapValues{Complex($0)})
        XCTAssertEqual(C.diagonals[1]!.elements, [0:1, 1:3, 2:5].mapValues{Complex($0)})
        
        XCTAssertNil(C.diagonals[-3])
        XCTAssertNil(C.diagonals[2])
        XCTAssertNil(C.diagonals[3])
        
        
    }
    
    func test_diagonalSparseTensorProduct() throws {
        
        
        let leftSubSpace = VectorSpace(dimension: 2, label: "")
        let rightSubSpace = VectorSpace(dimension: 3, label: "")
        let fullSpace = VectorSpace(tensorProductOf: [leftSubSpace, rightSubSpace],
                                            label: "")
        let A = DiagonalSparseMatrix(in: leftSubSpace, diagonals: [0: MatrixDiagonal(dimension: 2, diagIdx: 0, elements: [0:Complex(1)]),
                                                                   1: MatrixDiagonal(dimension: 2, diagIdx: 1, elements: [0:Complex(1)])])
        
        let B = DiagonalSparseMatrix(in: rightSubSpace, diagonals: [-1: MatrixDiagonal(dimension: 3, diagIdx: -1, elements: [1:1, 2:2].mapValues{Complex($0)}),
                                                                     2: MatrixDiagonal(dimension: 3, diagIdx: 2, elements: [0:Complex(1)])])
        
        let C = fullSpace.tensorProduct(of: A, with: B)
        
        XCTAssertEqual(C.diagonals[2]!.elements, [0:1, 1:1, 2:2].mapValues{Complex($0)})
        XCTAssertEqual(C.diagonals[5]!.elements, [0:1].mapValues{Complex($0)})
        XCTAssertEqual(C.diagonals[-1]!.elements, [1:1, 2:2].mapValues{Complex($0)})
        
        for diagIdx in [-5,-4,-3,-2,0,1,3,4] {
            XCTAssertNil(C.diagonals[diagIdx])
        }
        
    }
    
    func test_matrixFromDiagonalSparse() throws {
        
        let testSpace = VectorSpace(dimension: 4, label: "")
        
        let denseMatrix = Matrix(elements: [1.0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16], in: testSpace)
        
        let diagonalSparse = DiagonalSparseMatrix(from: denseMatrix)
        
        let reproducedDenseMatrix = Matrix(from: diagonalSparse)
        
        XCTAssertEqual(reproducedDenseMatrix, denseMatrix)
        
        
    }
    
    func test_diagonalSparseIndexing() throws {
        
        let testSpace = VectorSpace(dimension: 4, label: "")
        
        let denseMatrix = Matrix(elements: [1.0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16], in: testSpace)
        
        let diagonalSparse = DiagonalSparseMatrix(from: denseMatrix)
        
        for row in 0..<4 {
            for col in 0..<4 {
                
                XCTAssertEqual(denseMatrix[row,col], diagonalSparse[row,col])
            }
        }
    }
    
    func test_sparseIndexing_get() throws {
        
        let testSpace = VectorSpace(dimension: 4, label: "")
        
        let denseMatrix = Matrix(elements: [1.0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16], in: testSpace)
        
        let sparse_mat = SparseMatrix(from: denseMatrix)
        
        for row in 0..<4 {
            for col in 0..<4 {
                
                XCTAssertEqual(denseMatrix[row,col], sparse_mat[row,col])
            }
        }
    }
    
    
    func testConvertIntToBasisState() throws {
        
        /*
         Case 1: First subsystem of dimension 3, second of dimension 2
         Testing all 6 integers possible
         */
        XCTAssertEqual(convertIntToBasisState(intToConvert: 0, basisStateDimensions: [3,2]), [0,0])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 1, basisStateDimensions: [3,2]), [0,1])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 2, basisStateDimensions: [3,2]), [1,0])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 3, basisStateDimensions: [3,2]), [1,1])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 4, basisStateDimensions: [3,2]), [2,0])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 5, basisStateDimensions: [3,2]), [2,1])
        
        
        /*
         Case 2: 3 subsystems, first of dimension 2, second of dim 3, third of dim 4
         Testing for numbers 0,2,4,6,8,10,12,14,16, 19, 23
         */
        
        
        XCTAssertEqual(convertIntToBasisState(intToConvert: 0, basisStateDimensions: [2,3,4]), [0,0,0])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 2, basisStateDimensions: [2,3,4]), [0,0,2])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 4, basisStateDimensions: [2,3,4]), [0,1,0])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 6, basisStateDimensions: [2,3,4]), [0,1,2])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 8, basisStateDimensions: [2,3,4]), [0,2,0])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 10, basisStateDimensions: [2,3,4]), [0,2,2])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 12, basisStateDimensions: [2,3,4]), [1,0,0])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 14, basisStateDimensions: [2,3,4]), [1,0,2])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 16, basisStateDimensions: [2,3,4]), [1,1,0])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 19, basisStateDimensions: [2,3,4]), [1,1,3])
        XCTAssertEqual(convertIntToBasisState(intToConvert: 23, basisStateDimensions: [2,3,4]), [1,2,3])
        
        
        
    }
   
    
    func test_SEPC_tensorProdBasisVecWithLeftIdentity() throws {
        let identityDim = 3
        let basisDim = 4
        
        XCTAssertEqual(SEPC_tensorProdBasisVecWithLeftIdentity(identityDimension: identityDim, basisStateDimension: basisDim, basisStateIdx: 0), [0,4,8])
        XCTAssertEqual(SEPC_tensorProdBasisVecWithLeftIdentity(identityDimension: identityDim, basisStateDimension: basisDim, basisStateIdx: 1), [1,5,9])
        XCTAssertEqual(SEPC_tensorProdBasisVecWithLeftIdentity(identityDimension: identityDim, basisStateDimension: basisDim, basisStateIdx: 2), [2,6,10])
        XCTAssertEqual(SEPC_tensorProdBasisVecWithLeftIdentity(identityDimension: identityDim, basisStateDimension: basisDim, basisStateIdx: 3), [3,7,11])
    }
    
    func test_SEPC_nonSquareTensorProdWithRightIdentity() throws {
        
        let lhsSEPC1 = [1,4]
        let identityDim1 = 2
        XCTAssertEqual(SEPC_nonSquareTensorProdWithRightIdentity(identityDimension: identityDim1, lhs: lhsSEPC1), [2,3,8,9])
        
        let lhsSEPC2 = [0,2,7,4]
        let identityDim2 = 3
        XCTAssertEqual(SEPC_nonSquareTensorProdWithRightIdentity(identityDimension: identityDim2, lhs: lhsSEPC2), [0,1,2,6,7,8,21,22,23,12,13,14])
    }
    
    func testLinearFitting() throws {
        
        let x1Data = [0.0,1.0,2.0,3.0,4.0,5.0]
        let y1Data = x1Data.map{$0 * 2}
        
        
        let (slope1, intercept1) = getLinearFitCoefficientsFromLeastSquaresMethod(x1Data, y1Data)
        
        let y2Data = x1Data.map{$0 * 0.5 + 1}
        
        let (slope2, intercept2) = getLinearFitCoefficientsFromLeastSquaresMethod(x1Data, y2Data)
        
        XCTAssertEqual(slope1, 2.0)
        XCTAssertEqual(intercept1, 0.0)
        
        XCTAssertEqual(slope2, 0.5)
        XCTAssertEqual(intercept2, 1)
    }
    
    func test_partialTraceMatrixTwoQubit() throws {
        
        let leftSpace = VectorSpace(dimension: 2, label: "")
        let rightSpace = VectorSpace(dimension: 2, label: "")
        
        let totalSpace = VectorSpace(tensorProductOf: [leftSpace, rightSpace], label: "")
        
        let densityMatrix1 = Matrix(elements:[Double](repeating: 0.25, count: totalSpace.dimension*totalSpace.dimension),in: totalSpace)
        
        let reducedDensityMatrixLeft1 = densityMatrix1.traceInto(subSpace: leftSpace)
        let reducedDensityMatrixRight1 = densityMatrix1.traceInto(subSpace: rightSpace)
        
        XCTAssert(reducedDensityMatrixLeft1.elements == [0.5,0.5,0.5,0.5].map{Complex($0)})
        XCTAssert(reducedDensityMatrixRight1.elements == [0.5,0.5,0.5,0.5].map{Complex($0)})
        
        var densityMatrix2 = Matrix(in: totalSpace)
        
        densityMatrix2[0,0] = Complex(0.5)
        densityMatrix2[0,3] = Complex(0.5)
        densityMatrix2[3,0] = Complex(0.5)
        densityMatrix2[3,3] = Complex(0.5)
        
        let reducedDensityMatrixLeft2 = densityMatrix2.traceInto(subSpace: leftSpace)
        let reducedDensityMatrixRight2 = densityMatrix2.traceInto(subSpace: rightSpace)
       
        XCTAssert(reducedDensityMatrixLeft2.elements == [0.5, 0.0, 0.0, 0.5].map{Complex($0)})
        XCTAssert(reducedDensityMatrixRight2.elements == [0.5, 0.0, 0.0, 0.5].map{Complex($0)})
        
        
    }
    
    func test_partialTraceMatrixThreeQubit() throws {
        
        let spaceA = VectorSpace(dimension: 2, label: "")
        let spaceB = VectorSpace(dimension: 2, label: "")
        let spaceC = VectorSpace(dimension: 2, label: "")
        let totalSpace = VectorSpace(tensorProductOf: [spaceA, spaceB, spaceC], label: "")
        
        var densityMatrixABC = Matrix(in: totalSpace)
        
        //rho = (1/2)*(|GHZ><GHZ| + |W><W|)
        // |W> = (1/sqrt(3))(|001> + |010> + |100>)
        //|GHZ> = (1/sqrt(2))(|000> + |111>)
        densityMatrixABC[0,0] = Complex(1.0/4.0)
        densityMatrixABC[0,7] = Complex(1.0/4)
        
        densityMatrixABC[1,1] = Complex(1.0/6)
        densityMatrixABC[1,2] = Complex(1.0/6)
        densityMatrixABC[1,6] = Complex(1.0/6)
        
        densityMatrixABC[2,1] = Complex(1.0/6)
        densityMatrixABC[2,2] = Complex(1.0/6)
        densityMatrixABC[2,6] = Complex(1.0/6)
        
        densityMatrixABC[6,1] = Complex(1.0/6)
        densityMatrixABC[6,2] = Complex(1.0/6)
        densityMatrixABC[6,6] = Complex(1.0/6)
        
        densityMatrixABC[7,0] = Complex(1.0/4)
        densityMatrixABC[7,7] = Complex(1.0/4)
        
        //tracing out middle space
        let spacesAC = totalSpace.createSubspace(spacesToKeep: spaceA, spaceC, label: "")
        let reducedDensityAC = densityMatrixABC.traceInto(subSpace: spacesAC)
        
        let correctElements = [5.0/12, 0.0, 1.0/6, 0, 0, 1.0/6, 0, 0, 1.0/6, 0, 1.0/6, 0, 0, 0, 0, 1.0/4]
        
        for i in 0..<correctElements.count {
            XCTAssertEqual(reducedDensityAC.elements[i].real, correctElements[i], accuracy: 0.00000001)
            
        }
    }
    
    
    func testPartialTraceMultiQubit() throws {
        // Testing for 6 qubits tracing out three. Product state of each qubit in |𝜓> = (1/sqrt(2))(|0> + |1>)
        //After trace out should have an 8x8 matrix where every element is 1/8.
        let numQubits = 6
        let indicesToTraceOut = [0,3,5]
        
        var dm_val = 1.0
        var totalSpaceDim = 1
        
        for _ in 0..<numQubits {
            dm_val /= 2.0
            totalSpaceDim *= 2
        }
        
        let dm_elements = [Double](repeating: dm_val, count: totalSpaceDim*totalSpaceDim)
        
        var subspaces: [VectorSpace] = []
        
        for i in 0..<numQubits {
            subspaces.append(VectorSpace(dimension: 2, label: "Qubit \(i) subspace"))
        }
        
        let totalSpace = VectorSpace(tensorProductOf: subspaces, label: "total qubit space")
        
        let densityMatrix = Matrix(elements: dm_elements, in: totalSpace)
        
        let subspaceToKeep = totalSpace.createSubspace(removeSpacesAtIndices: indicesToTraceOut, label: "subspace to keep")
        
        let reduced_dm = densityMatrix.traceInto(subSpace: subspaceToKeep)
        
        for elem in reduced_dm.elements {
            XCTAssertEqual(elem.real, 1.0/8.0, accuracy: 0.000001)
            XCTAssertEqual(elem.imag, 0.0, accuracy: 0.000001)
        }
    }
    
    func testPartialTraceTwoQubitSparse() throws {
        
        /*
         rho =  (0.5*|0><0| + 0.5 * |1><1|)_A ⊗ (0.5 * |0><0| + |1><1|)_B
         */
        let leftSpace = VectorSpace(dimension: 2, label: "")
        let rightSpace = VectorSpace(dimension: 2, label: "")
        
        let totalSpace = VectorSpace(tensorProductOf: [leftSpace, rightSpace], label: "")
        
        var elems: [CoordinateStorage] = []
        
        for i in 0..<totalSpace.dimension {
            elems.append(CoordinateStorage(value: 1.0/Double(totalSpace.dimension), row: i, col: i))
        }
        
        let densityMatrix1 = SparseMatrix(values: elems, in: totalSpace)
        
        let reducedDensityMatrixLeft1 = densityMatrix1.traceInto(subSpace: leftSpace)
        let reducedDensityMatrixRight1 = densityMatrix1.traceInto(subSpace: rightSpace)
        
        print("\nleft: ", reducedDensityMatrixLeft1)
        print("\nright: ",reducedDensityMatrixRight1, "\n")
 
    
        
    }

    
    func test_accelerate_mmul() throws {
        
        let space = VectorSpace(dimension: 2, label: "Accelerate mmul test space")
        
        let A = Matrix(elements: [Complex(1.0), Complex(2.0), Complex(3.0), Complex(4.0)], in: space)
        let B = Matrix(elements: [Complex(5.0), Complex(6.0), Complex(7.0), Complex(8.0)], in: space)
        
        let C = A.accelerateMatrixMultiplication(B)
        
        XCTAssertEqual(C.elements, [Complex(19),Complex(22),Complex(43),Complex(50)])
        
        
    }
    
    
    func test_Accelerate_matrixTimesVector() throws {
        let space = VectorSpace(dimension: 2, label: "test Space")
        var matrix = Matrix(in: space)
        matrix[0,0] = Complex(1.0)
        matrix[0,1] = Complex(2.0)
        matrix[1,0] = Complex(3.0)
        matrix[1,1] = Complex(4.0)
        var vector = Vector(in: space)
        vector[0] = Complex(5.0)
        vector[1] = Complex(6.0)
        
        let res = matrix.accelerateVectorMult(vector)
        
        
        XCTAssertTrue(res.elements == [Complex(17.0),Complex(39.0)] )
    }
    
    func test_Accelerate_matrixTimesVector_new_rep() throws {
        let space = VectorSpace(dimension: 2, label: "test Space")
        
        let lhsReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: space.dimension*space.dimension)
        _ = lhsReal.initialize(from: [1,2,3,4])
        
        let lhsImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: space.dimension*space.dimension)
        _ = lhsImag.initialize(from: [0,0,0,0])
        
        let lhsElem = DSPDoubleSplitComplex(realp: lhsReal.baseAddress!, imagp: lhsImag.baseAddress!)
        
        
        
        var matrix = DSPComplexMatrix(elements: lhsElem, in: space)
        var vector = Vector(in: space)
        vector[0] = Complex(5.0)
        vector[1] = Complex(6.0)
        
        let res = matrix.accelerateVectorMult(vector)
        
        
        XCTAssertTrue(res.elements == [Complex(17.0),Complex(39.0)] )
    }
    
    func test_DSP_matrix_add() throws {
        let space = VectorSpace(dimension: 2, label: "test Space")
        
        let lhsReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: space.dimension*space.dimension)
        _ = lhsReal.initialize(from: [1,2,3,4])
        
        let lhsImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: space.dimension*space.dimension)
        _ = lhsImag.initialize(from: [5,6,7,8])
        
        let lhsElem = DSPDoubleSplitComplex(realp: lhsReal.baseAddress!, imagp: lhsImag.baseAddress!)
        
        var lhs = DSPComplexMatrix(elements: lhsElem, in: space)
        
        let rhsReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: space.dimension*space.dimension)
        _ = rhsReal.initialize(from: [1,2,3,4])
        
        let rhsImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: space.dimension*space.dimension)
        _ = rhsImag.initialize(from: [5,6,7,8])
        
        let rhsElem = DSPDoubleSplitComplex(realp: rhsReal.baseAddress!, imagp: rhsImag.baseAddress!)
        
        var rhs = DSPComplexMatrix(elements: rhsElem, in: space)
        
        let res = lhs + rhs
        
        let res_real = [2.0,4,6,8]
        let res_imag = [10.0, 12, 14, 16]
        
        for i in 0..<space.dimension * space.dimension {
            XCTAssertEqual(res.elements.realp[i], res_real[i])
            XCTAssertEqual(res.elements.imagp[i], res_imag[i])
        }
    }
    
    func test_DSP_matrix_subtract() throws {
        let space = VectorSpace(dimension: 2, label: "test Space")
        
        let lhsReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: space.dimension*space.dimension)
        _ = lhsReal.initialize(from: [1,2,3,4])
        
        let lhsImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: space.dimension*space.dimension)
        _ = lhsImag.initialize(from: [5,6,7,8])
        
        let lhsElem = DSPDoubleSplitComplex(realp: lhsReal.baseAddress!, imagp: lhsImag.baseAddress!)
        
        var lhs = DSPComplexMatrix(elements: lhsElem, in: space)
        
        let rhsReal = UnsafeMutableBufferPointer<Double>.allocate(capacity: space.dimension*space.dimension)
        _ = rhsReal.initialize(from: [1,2,3,4])
        
        let rhsImag = UnsafeMutableBufferPointer<Double>.allocate(capacity: space.dimension*space.dimension)
        _ = rhsImag.initialize(from: [5,6,7,8])
        
        let rhsElem = DSPDoubleSplitComplex(realp: rhsReal.baseAddress!, imagp: rhsImag.baseAddress!)
        
        var rhs = DSPComplexMatrix(elements: rhsElem, in: space)
        
        let res = lhs - rhs
        
        for i in 0..<space.dimension * space.dimension {
            XCTAssertEqual(res.elements.realp[i], 0.0)
            XCTAssertEqual(res.elements.imagp[i], 0.0)
        }
    }
    
}




//  Created by  M J Everitt on 17/01/2022.

