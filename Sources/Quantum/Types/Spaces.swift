/*
 
 See Adding Quantum Structure to the Code: Spaces (section 8.4.3.1 of the book)
 The purpose of this code is very much to be able to map the mathematical structure
 of problems we are trying to solve into code and to prevent human error.
 
 So for example, thinking of any matrix as an operator in a vector space we do not
 want e.g. to be allowed to add two matrices from separate vector spaces together or
 have a matrix from one vector space act on a vector in another space even if they
 are appropriately dimensioned.
 
 In addition, if we make tensor products of spaces we want the code to be able to check
 that the user can consistently work with operators in that space.
 
 TODO: Optimise
 TODO: Nested Kronecker products
 For
    1) tensor products of multiple spaces
    2) tensor products of tensor product spaces
 TODO: Add logging to spaces and its members, such as vectors
 */
import Foundation

fileprivate var space_counter = 0
// Note typealias ScalarField is important rather than <T> to make clear
// closure under scalar arithmetic as we have several different types working together.

public class VectorSpace: Hashable  {
    public let dimension: Int
    public let description: String
    internal let identifier: Int
    internal var setofSpaces: [VectorSpace]


//    public typealias SpaceVector = Vector<T>

    public init(dimension: Int, label: String) {
        self.dimension = dimension
        self.description = label
        self.identifier = space_counter
        space_counter += 1
        setofSpaces = []
        self.setofSpaces.append(self)
    }

    public static func == (lhs: VectorSpace, rhs: VectorSpace) -> Bool {
        return lhs.identifier == rhs.identifier
    }
    public func hash(into hasher: inout Hasher) {
        hasher.combine(identifier)
        hasher.combine(setofSpaces)
    }
    
    public var identityOperator: Matrix {
        return  makeIdentityMatrix()
    }
    
    public func makeIdentityMatrix() -> Matrix {
        var output = Matrix(in: self)
        for i in 0 ..< dimension {
            output[i, i] = Complex(1)
        }
        return output
    }
    
    
    public func createSubspace(spacesToKeep: VectorSpace..., label: String, sorted: Bool = false) -> VectorSpace {
        
        var sortedSpacesToKeep = spacesToKeep
        if !sorted {
            sortedSpacesToKeep = spacesToKeep.sorted(by: {$0.identifier < $1.identifier})
        }
        
        assert(sortedSpacesToKeep != setofSpaces, "Cannot create a subspace which is the total space.")
        
        //TODO: check each space in spacesToKeep is a subspace of the totalSpace. For now, assume that it is.
        if sortedSpacesToKeep.count == 1 {
            return sortedSpacesToKeep[0]
        }
        
//        TODO: somehow check if the space which is exactly a combination of these spaces already exists, if it does,
//        return that space. Not sure how this would be done without storing every space which has been defined.
//        Difficult problem due to the generics, cannot predefine the types of elements of a generic list.
        return VectorSpace(tensorProductOf: sortedSpacesToKeep, label: label, sorted: true)

    }
    
    public func createSubspace(removeSpacesAtIndices: [Int], label: String) -> VectorSpace {
        
        var spacesToKeep: [VectorSpace] = []
        
        for i in 0..<setofSpaces.count {
            if removeSpacesAtIndices.contains(i) { continue }
            spacesToKeep.append(setofSpaces[i])
        }
        
        return VectorSpace(tensorProductOf: spacesToKeep, label: label, sorted: true)
        
        
    }
}

// MARK: - Add Tensor product capability
extension VectorSpace {
    convenience public init(tensorProductOf spaces: [VectorSpace],
                            label: String, sorted: Bool = false) {
        assert(spaces.count > 1, "Cannot make a tensor product of one space")
        
        var tempDimension = 1
        for space in spaces {
            tempDimension *= space.dimension;
        }
        
        var sortedSpaces = spaces
        if !sorted {
            sortedSpaces = spaces.sorted(by: {$0.identifier < $1.identifier} )
        }
        
       for i in 0 ..< sortedSpaces.count-1 {
            assert(sortedSpaces[i] != sortedSpaces[i+1],
                   "Error: The same space included more than once in tensor product space")
        }

        var componentSpaceLabels = "\n"
        for space in sortedSpaces {
            componentSpaceLabels.append(contentsOf: space.description)
        }

        self.init(dimension: tempDimension, label: label + " (tensor product space)")
        setofSpaces = sortedSpaces
    }
    
    
    
    
    public func tensorProduct(of A: Matrix, with B: Matrix) -> Matrix {
        
        assert (setofSpaces.count == 2, "Tensor product currently only implemented for two operators not in tensor product spaces")

        let temp = [A,B].sorted(by: { $0.space.identifier < $1.space.identifier} )

        for i in 0 ..< setofSpaces.count {
            assert (temp[i].space == setofSpaces[i], "operator in wrong space " + A.space.description)
        }
        assert(temp[0].elements.count == (temp[0].space.dimension * temp[0].space.dimension), "\(temp[0].elements.count) != \(temp[0].space.dimension)^2" )
        assert(temp[1].elements.count == (temp[1].space.dimension * temp[1].space.dimension), "\(temp[1].elements.count) != \(temp[1].space.dimension)^2" )
        let C = kronekerProduct(A: temp[0].elements,
                                rowsA: temp[0].space.dimension,
                                colsA: temp[0].space.dimension,
                                B: temp[1].elements,
                                rowsB: temp[1].space.dimension,
                                colsB: temp[1].space.dimension)
        
        return Matrix(elements: C, in: self)
    }

    public func tensorProduct(of A: DiagonalSparseMatrix,
                              with B: DiagonalSparseMatrix)
    -> DiagonalSparseMatrix {
        
        assert (setofSpaces.count == 2, "Tensor product currently only implemented for two operators not in tensor product spaces")

        let temp = [A,B].sorted(by: { $0.space.identifier < $1.space.identifier} )

        for i in 0 ..< setofSpaces.count {
            assert (temp[i].space == setofSpaces[i], "operator in wrong space " + A.space.description)
        }
        
        var output  = DiagonalSparseMatrix(in: self)
        
        for (lhsDiagIdx, lhsDiag) in A.diagonals {
            for (rhsDiagIdx, rhsDiag) in B.diagonals {
                
                let currOutputDiagIdx = B.space.dimension * lhsDiagIdx + rhsDiagIdx
                let currOutputDiag = lhsDiag.tensorProductWith(rhsDiag)
                
                
                if let _ = output[currOutputDiagIdx] {
                    
                    output[currOutputDiagIdx]! += currOutputDiag
                }
                
                else {
                    output[currOutputDiagIdx] = currOutputDiag
                }
            }
        }
        return output
        
        
    }
    
    public func tensorProduct(of A: Vector,
                              with B: Vector)
    -> Vector {
        assert (setofSpaces.count == 2, "Tensor product currently only implemented for two operators not in tensor product spaces")

        let temp = [A,B].sorted(by: { $0.space.identifier < $1.space.identifier} )

        for i in 0 ..< setofSpaces.count {
            assert (temp[i].space == setofSpaces[i], "operator in wrong space " + A.space.description)
        }
        
        assert(temp[0].elements.count == temp[0].space.dimension, "\(temp[0].elements.count) != \(temp[0].space.dimension)")
        assert(temp[1].elements.count == temp[1].space.dimension, "\(temp[1].elements.count) != \(temp[1].space.dimension)")
        
        let C = kronekerProduct(A: temp[0].elements,
                                rowsA: temp[0].space.dimension,
                                colsA: 1,
                                B: temp[1].elements,
                                rowsB: temp[1].space.dimension,
                                colsB: 1)
        
        return Vector(elements: C, in: self)
    }
}
//  Created by M J Everitt on 18/01/2022.
