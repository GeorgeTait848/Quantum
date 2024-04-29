/*
 Standard utility functions used by different structs and classes (to make code DRY).
 
 e.g. scalarBinaryOperation allows us to multiply, divide etc operators and vectors element wise by scalar
 
 TODO: Optimise
 It is not clear that map and zip are as fast as one might expect.
 Test performance against for loop and consider adding an automated way to switch
 depending on e.g. size of arrays or platform.
 
 TODO: Nested Kronecker products
 For
    1) tensor products of multiple spaces
    2) tensor products of tensor product spaces
 
 TODO: Consider if arrays could be replaced by Collection in `elementwise' funcs to make more general
 */
import Foundation

public func scalarBinaryOperation<T>(array: [T],by value: T, operation: (T,T)->T) -> [T] {
    return array.map( { (element: T) -> T in return operation(element, value) } )
}

public func elementWisePrefixOperation<T>(array: [T], operation: (T)->T) -> [T] {
    return array.map( { (element: T) -> T in return operation(element) } )
}

public func elementWiseBinaryOperation<T>(thisArry: [T], thatArray: [T], operation: (T,T)->T ) -> [T] {
    assert(thisArry.count == thatArray.count)
    return zip(thisArry,thatArray).map(operation)
}

func checkInSameSpace<U: livesInAVectorSpace,V: livesInAVectorSpace>(_ lhs: U,_ rhs: V) {
    assert(lhs.space.identifier == rhs.space.identifier,
           """
              Incompatable spaces:
                  - \(lhs.space.description)
                  - \(rhs.space.description)
              """)
}
public func repeatedly<T>(apply binaryFunction: (T,T)->T, _ items: [T]) ->T {
    assert(items.count > 1, "need more than one item")
    var output = items[0]
    for i in 1 ..< items.count {
        output = binaryFunction( output , items[i] )
    }
    return output
}
public func repeatedly<T>(apply binaryFunction: (T,T)->T, _ items: T ...) ->T {
    return repeatedly(apply: binaryFunction, items)
}
public func sum<A: Addable>(_ vectors: A...) -> A {
    repeatedly(apply: +, vectors)
}

// For dereferencing a one-d arras as if its a matrix.
//      row * dim + col   is column major
//      row  + col * dim  is row major
// currently format is not made clear in self-documenting code.
// TODO: consider adding a boolean "column major" argument to make this explicit.
// (maybe with default value true)
func atIndex(row: Int, column: Int, nColumns: Int) -> Int { 
    return (row * nColumns) + column
} 

// kronecker delta
// TODO: Consider renaming
func delta(_ n: Int,_ m: Int) -> Int {
    var out = 0
    if n == m {
        out = 1
    }
    return out
}

// https://en.wikipedia.org/wiki/Kronecker_product
// [accessed: 12/01/2022 - see Definition]
public func kronekerProduct  (
    A: [Complex], rowsA: Int, colsA :Int,
    B: [Complex], rowsB: Int, colsB :Int)
-> [Complex]  {
    assert(A.count == rowsA * colsA, "dimension of A bad: dim(A) = \(A.count), rowsA = \(rowsA), colsB = \(colsA)")
  assert(B.count == rowsB * colsB, "dimension of A bad: dim(A) = \(B.count), rowsA = \(rowsB), colsB = \(colsB)")
    var C = Array(repeating: Complex(real: 0.0),
                count: rowsA*rowsB*colsA*colsB)
    
  let colsC = colsA * colsB
  let p = rowsB // so notation is same as wikipedia
  let q = colsB

  let index = { (_ R: Int,_ C: Int,_ N: Int) -> Int in
        return atIndex(row: R, column: C, nColumns: N) }

  for r in 0 ..< rowsA {
    for s in 0 ..< colsA {
      let A_rs = A[index(r,s, colsA)]
      for v in 0 ..< rowsB {
        for w in 0 ..< colsB {
          C[index(p * r + v, q * s + w, colsC)] =
                            A_rs * B[index(v,w, colsB)]
        }
      }
    }
  }
  return C
}


public func getRowLimits(dim: Int, diagIdx: Int) -> (Int, Int) {
    
    if diagIdx <= 0 {
        return (-diagIdx, dim - 1)
    }
    
    return (0, dim - diagIdx - 1)
    
}


/*
 The following functions are subroutines for the partial trace algorithm.
 The partial trace algorithm is a generalisation of the mathematics described in
 https://www.ryanlarose.com/uploads/1/1/5/8/115879647/quic06-states-trace.pdf
 Eqns. (6.23),(6.24).
 
 Which is specific to two subsystems. This algorithm works for an arbitrary number of subsystems,
 with any number of subspaces being traced out, at any position in the total space.
 
 The partial trace algorithm uses the fact that the matrices that pre and post multiply the density matrix have
 the following properties:
 
 - The postmultiplying matrix has a single non-zero element per column, and this element always has value 1. This property
 will be referred to in the code as SEPC (Single Element Per Column). This property ensures that the postmultiplication of
 the density matrix by SEPC matrix will output a matrix in which each element will only have a single contribution from the density matrix.
 That is, it removes the summation over the column of the rhs matrix of matrix multiplication.
 
 - The premultiplying matrix is the transpose of the post multiplying one, hence it will have a single element per row.
 This allows for the highly efficicent method of "index selection".
 
 Index selection means that each basis state will "select" exactly one element of the density matrix for an element of the
 reduced density matrix. Hence, denoting the basis state of the total space being traced out as 𝛼 and its corresponding contribution
 to the reduced density matrix as 𝜌_𝛼. Then, the reduced density matix, ρ_r is simply:
 
 ρ_r = ∑_𝛼 ρ_𝛼
 
 Where, ρ_𝛼 is defined as follows. Given that any basis state 𝛼 will result in a SEPC matrix which postmultiplies ρ, and whose transpose
 premultiplies ρ. We can store the SEPC as a list of indices. Each element refers to the row at which the nonzero element is found for each
 column in the SEPC. Denote this list as L. Then, one can show quite simply with linear algebra that:
 
 ρ_𝛼[i,j] = ρ[L[i], L[j]]
 
 So, the entirety of the partial trace algorithm is determining the list L. One can show that there are only 4 algorithms required to determine L.
 These are:
 
 1) Tensor product of a basis vector with identity matrix from the left using SEPC. This is used when the first space to trace out is not the first space
 in the total space. Let the identity be of dimension n and the basis state of dimension m, let the index of the nonzero element of the basis state be i.
 
 L[j] = i + jm  ,  j ∈ [0,n-1]
 
 
 2) Tensor product of nonsquare SEPC matrix with identiy from right. This is used when there are two subspaces being traced out, which are NOT consecutive
 in the total space. Hence, we need to tensor product with all identities between these two spaces. Note that the case where the first space is traced out
 but the second is not is a special case of this algorithm and can be handled with the same logic (this would be the tensor product of basis vector with
 identity from left).
 
 Let the SEPC (lhs) be of dimension nm x m, and the identity (rhs) be of dimension p x p. Then:
 
 L[j] = lhs[j/p] * p + j % p    ,   j ∈ [0, mp-1]
 
 3) Tensor product of nonsquare SEPC (lhs) with basis state (rhs). Used when tracing out a basis state in the space which is not the first space to be traced out.
 Note that tensor product of consecutive basis states, where the first of which is the first space to be traced out and is the first space in the total space, is
 a special case of this algorithm and follows the exact same logic.
 
 Let the SEPC (lhs) be of dimensino nm x m and the basis state (rhs) be of dimension p and index i. Then,
 
 L[j] = lhs[j] * p + i  ,   j ∈ [0,m-1]
 
 4) Tensor product of identity matrices in SEPC logic. This is used to tensor product the identities between nonconsecutive spaces to be traced out.
 For any identity matrix L[j] = j.
 For the tensor product of any two: I_n ⊗ I_m = I_(n x m)
 So, one only needs to track the change in dimension here.
 
 
 I recommend that anybody struggling to understand these algorithms does some concrete examples of bipartite systems of different dimension hilbert spaces to see
 that these actually work and how they are used to select the elements.
 */

public func convertIntToBasisState(intToConvert n: Int, basisStateDimensions b: [Int]) -> [Int] {
    /*
     Function to be used in partial trace algorithm.
     When tracing out an arbitrary number of spaces, each with arbitrary dimension hilbert space,
     it is necessary to ensure all possible combinations of basis states being traced out are accounted for.
     
     A set of spaces {H_i} will have total dimension d = 𝛱_i(dim(H_i)), hence d basis states to trace out.
     
     This function takes any n ∈ [0,d-1] and returns the basis state index for each subspace.
     Adapted from the algorithm for decimal to binary conversion at:
     
     https://www.geeksforgeeks.org/python-decimal-to-binary-list-conversion/
     
     Specifically, "Method 3: Using While"
     */
    var n_copy = n
    var i = b.count-1
    var output = [Int](repeating: 0, count: b.count)
    
    while n_copy > 0 {
        output[i] = n_copy % b[i]
        n_copy /= b[i]
        i -= 1
    }
    
    return output
}


public func SEPC_tensorProdBasisVecWithLeftIdentity(identityDimension: Int, basisStateDimension: Int, basisStateIdx: Int) -> [Int] {
    
    var output = [Int](repeating: 0, count: identityDimension)
    
    for j in 0..<output.count {
        output[j] = j * basisStateDimension + basisStateIdx
    }
    
    return output
}

public func SEPC_nonSquareTensorProdWithRightIdentity(identityDimension: Int, lhs: [Int]) -> [Int] {
    
    let lhsNumCols = lhs.count
    var output = [Int](repeating: 0, count: identityDimension * lhsNumCols)
    
    for j in 0..<output.count {
        output[j] = lhs[j/identityDimension] * identityDimension + j % identityDimension
    }
    
    return output
}


public func SEPC_matrixTensorProdWithRightBasisState(basisStateDimension: Int, basisStateIdx: Int, lhs: [Int]) -> [Int] {
    
    return lhs.map{$0 * basisStateDimension + basisStateIdx}
}

//  Created by M J Everitt on 17/01/2022.

public func getBasisStateContributionToPartialTrace<T: OperatorType> (traceInto subspace: VectorSpace,
                                                                      fullMatrix: T,
                                                                      selectionIndices: [Int]) -> T {
    
    var output = T(in: subspace)
    
    for row in 0..<subspace.dimension {
        for col in 0..<subspace.dimension {
            output[row,col] = fullMatrix[selectionIndices[row], selectionIndices[col]]
        }
    }
    
    return output
    
}

public func getBasisStatePartialTraceSelectionIndices(totalSpaceDimensions: [Int],
                                                      traceOutSpacesAtIdx: [Int],
                                                      currentBasisState: [Int]) -> [Int] {
    
    let traceOutSpacesDimensions = traceOutSpacesAtIdx.map{totalSpaceDimensions[$0]}
    var SEPC_indices: [Int]? = nil
    var identityDimension = 1
    var traceOutNextIdx = 0
    
    for i in 0..<totalSpaceDimensions.count {
        
        //all required spaces traced out, include any tensor prods of remaining spaces.
        if traceOutNextIdx == traceOutSpacesAtIdx.count {
            for j in i..<totalSpaceDimensions.count {
                identityDimension *= totalSpaceDimensions[j]
            }
            SEPC_indices = SEPC_nonSquareTensorProdWithRightIdentity(identityDimension: identityDimension, lhs: SEPC_indices!)
            break
        }
        
        //dont need to trace out this space, move to next
        if i != traceOutSpacesAtIdx[traceOutNextIdx] {
            identityDimension *= totalSpaceDimensions[i]
            continue
        }
        
        //do need to trace out this space, and it is the first space we need to trace out.
        if SEPC_indices == nil {
            
            let basisDimension = traceOutSpacesDimensions[traceOutNextIdx]
            SEPC_indices = SEPC_tensorProdBasisVecWithLeftIdentity(identityDimension: identityDimension,
                                                                   basisStateDimension: basisDimension ,
                                                                   basisStateIdx: currentBasisState[traceOutNextIdx])
            identityDimension = 1
            traceOutNextIdx += 1
            continue
        }
        
        
        //need to trace out this space and it is not the first space we are tracing out.
        
        if identityDimension > 1 {
            //i.e we have nonconsecutive spaces to trace out
            SEPC_indices = SEPC_nonSquareTensorProdWithRightIdentity(identityDimension: identityDimension, lhs: SEPC_indices!)
        }
        
        
        
        let basisDimension = traceOutSpacesDimensions[traceOutNextIdx]
        SEPC_indices = SEPC_matrixTensorProdWithRightBasisState(basisStateDimension: basisDimension,
                                                                basisStateIdx: currentBasisState[traceOutNextIdx],
                                                                lhs: SEPC_indices!)
        traceOutNextIdx += 1
        identityDimension = 1
        continue
    }
    
    return SEPC_indices!
}
