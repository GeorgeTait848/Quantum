//
//  File.swift
//
//
//  Created by George Tait on 05/10/2023.
//

import Foundation

public func benchmarkDiagonalSparseMatMult(dims: [Int], writeDataToFileAt: String) {
    
    let diagonalSparseMM: (DiagonalSparseMatrix, DiagonalSparseMatrix) -> DiagonalSparseMatrix  =  {A,B in A*B}
    let lc = 50
    
    var lowDensities = [Double](repeating: 0.0, count: lc)
    
    for i in 0..<lowDensities.count {
        lowDensities[i] += 0.1*Double(i)/Double(lc)
    }
    
    
    let hc = 18
    var highDensities = [Double](repeating: 0.0, count: hc)
    
    
    for i in 0..<hc {
        highDensities[i] += 0.1 + 0.9*Double(i)/Double(hc)
    }
    
    let densities = lowDensities + highDensities
}


public func makeMatrixByDiagonalDensity(in space: VectorSpace, densityOfDiagonals: Double) -> DiagonalSparseMatrix {
    
    var output = DiagonalSparseMatrix(in: space)
    
    let numDiags = Int(Double(2*space.dimension-1)*densityOfDiagonals)
    
    var n = 0; var currDiagIdx = 0
    
    while n < numDiags {
        let (rowLowerLim, rowUpperLim) = getRowLimits(dim: space.dimension, diagIdx: currDiagIdx)
        var elem: [Int:Complex] = [:]
        
        for row in rowLowerLim...rowUpperLim {
            elem[row] = Complex(real: 1, imag: 0)
        }
        
        let curr = MatrixDiagonal(dimension: space.dimension, diagIdx: currDiagIdx, elements: elem)
        output[currDiagIdx] = curr
        
        n += 1
        
        if currDiagIdx == 0 {
            currDiagIdx = 1
            continue
        }
        
        if currDiagIdx < 0 {
            currDiagIdx = 1-currDiagIdx
            continue
        }
        
        currDiagIdx *= -1
    }
    
    return output
    
}

public func estimateTimeComplexity(codeToExecute: (Int) -> Void, dims: [Int]) -> (Double, Double, Double, Double) {
    
    let elapsedTimes = getElapsedTimeData(codeToExecute: codeToExecute, dims: dims)
    let totalTime = elapsedTimes.reduce(0, +)
    let logElapsedTimes = elapsedTimes.map{ log($0) }
    
    let logDims = dims.map{ log(Double($0)) }

    let (omega, intercept) = getLinearFitCoefficientsFromLeastSquaresMethod(logDims, logElapsedTimes)
    
    let r_sq = getRSquaredCoefficient(logDims, logElapsedTimes)
    return (omega, r_sq, intercept, totalTime)
    
}


public func getElapsedTimeData(codeToExecute: (Int) -> Void, dims: [Int]) -> [Double]{
    
    var times = [Double](repeating: 0, count: dims.count)
    for i in 0..<dims.count {
        
        let startTime = clock()
        codeToExecute(dims[i])
        let endTime = clock()
        
        let elapsedTime = Double((endTime - startTime))/Double(CLOCKS_PER_SEC/1_000)
        
        times[i] = elapsedTime
    }
    
    
    return times
}




public func storeElapsedTimeDataToFile(codeToExecute: (Int) -> Void, dims: [Int], pathFromHomeToFile: String) {
    
    let outputText = convertElapsedTimeDataToWriteableFormat(codeToExecute: codeToExecute, dims: dims)
    
    let writeFilename = FileManager.default.homeDirectoryForCurrentUser.path + pathFromHomeToFile
    
    
    _ = FileManager.default.createFile(atPath: writeFilename, contents: nil, attributes: nil)
    
    do {
    
    try  outputText.write(toFile: writeFilename, atomically: false, encoding: .utf8) }
    
    catch { errorStream.write("Can not write to output \n")}
    
    
}

public func convertElapsedTimeDataToWriteableFormat(codeToExecute: (Int) -> Void, dims: [Int]) -> String {
    
    let times = getElapsedTimeData(codeToExecute: codeToExecute, dims: dims)
    var outputText = ""
    
    for i in 0..<dims.count {
        
        outputText += "\(dims[i]) \t \(times[i]) \n"
        
    }
    
    return outputText
}


public func getRSquaredCoefficient(_ x: [Double], _ y: [Double]) -> Double {
    
//  Method from https://en.wikipedia.org/wiki/Coefficient_of_determination
    
    let (estimatingSlope, estimatingIntercept) = getLinearFitCoefficientsFromLeastSquaresMethod(x, y)
    let ybar = y.reduce(0.0, +)/Double(y.count)
    
    var sumOfResSquares = 0.0
    var totalSumOfSquares = 0.0
    
    for i in 0..<x.count {
        
        let f_i = x[i]*estimatingSlope + estimatingIntercept
        
        sumOfResSquares += (y[i] - f_i) * (y[i] - f_i)
        totalSumOfSquares += (y[i] - ybar) * (y[i] - ybar)
        
    }
    
    return 1.0 - sumOfResSquares/totalSumOfSquares
    
}

public func getLinearFitCoefficientsFromLeastSquaresMethod(_ x: [Double], _ y: [Double]) -> (Double, Double) {
    
//    Method from: https://www.varsitytutors.com/hotmath/hotmath_help/topics/line-of-best-fit

    
    let xbar = x.reduce(0.0, +)/Double(x.count)
    let ybar = y.reduce(0.0, +)/Double(y.count)
    
    var numerator = 0.0
    var denom = 0.0
    
    for i in 0..<x.count {
        
        numerator += (x[i] - xbar)*(y[i] - ybar)
        denom += (x[i] - xbar)*(x[i] - xbar)
    }
    
    let slope = numerator/denom
    let intercept = ybar - slope*xbar
    
    return (slope,intercept)
    
}


