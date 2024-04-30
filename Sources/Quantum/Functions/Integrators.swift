/*
 The idea for this approach is that the same integrators can be used for very many
 different types. This might include collections of Doubles, Complex<Double> or any
 type that is a Scalar but also more exotic types as well such as quaternions over
 a rational number type, arbitrary precision arithmetic or even business model types
 in e.g. stock market modelling.
 
 This is an example of applying The Dependency Inversion Principle to scientific coding.
 
 This leads to a bit of a coding overhead at this stage but each algorithm should
 only need to be written once.
 
 Type safety is used here to make sure that the independent variable is of the right
 kind (we would not want to accidentally cast e.g. time to a complex type but as some do
 look at complex time we should still allow the possibility).
 
 Some of the material is based off Numerical Recipes in C versions of algorithms
 written by others. Numerical recipes state:
 
    If you analyse the ideas contained in a program, and then express those ideas
    in your own completely different implementation, then that new program implementation
    belongs to you. That is what we have done for those programs in this book
    that are not entirely of our own devising. When programs in this book are said to be
    “based” on programs published in copyright sources, we mean that the ideas are the
    same. The expression of these ideas as source code is our own. We believe that no
    material in this book infringes on an existing copyright.
 
 Re-expressing that code in terms of generics but preserving the traceability by using
 standard labels found in the literature is in the spirit of the above statement and
 so I believe that this code does not infringe on any existing copyright. This is because
 my implementation using generics to satisfy the The Dependency Inversion Principle makes this
 work as fundamentally different from that of Numerical Recipes as theirs was from the sources
 they used as a basis for their code. I also hope our use encourages sales of the Numerical
 Recipes books which I believe are a rather useful resource coving aspects of scientific
 coding that we have not covered.
 
 Also note that there is some odd behaviour in my implementation of odeint for non-linear
 systems. In the one week time frame it was not possible to determine if this was simply
 due to numerical instability or my implementation of the algorithm itself.

 TODO: Add more advanced integrators.
 TODO: Much more testing.
*/

import Foundation
import CoreLocation


public protocol AdaptiveSteppable: Integrable, EquatableUnderRelativeTolerance {
    
}

public protocol EquatableUnderRelativeTolerance {
    
    func  checkEquivalenceUnderRelativeTolerance (_ other: Self, relativeTol: Double) -> Bool
}

public protocol Integrable:
                  Addable & ClosedUnderScalarFieldMultiplication{}


public func adaptiveRungeKuttaOverRange <T: AdaptiveSteppable> (f: (T, Double) -> T,
                                                                y: inout T,  from x_0: Double,
                                                                to x_f: Double, h: Double,
                                                                relativeTol: Double) {
    
    var current_h = h
    var current_x = x_0
    
    
    while current_x < x_f {

        current_h = adaptStep(f: f, y: y, x: current_x, h: current_h, relativeTol: relativeTol)
        y = rungeKuttaFourthOrder(f: f, y: y, x: current_x, h: current_h)
        current_x += current_h
        
    }
        
    current_h = x_f - current_x
    y = rungeKuttaFourthOrder(f: f, y: y, x: current_x, h: current_h)
    
}




public func adaptiveRungeKuttaFourthOrder <T: AdaptiveSteppable> (f: (T, Double) -> T,
                                                                  y: T, x: Double,
                                                                  h: Double,
                                                                  relativeTol: Double) -> T {
    
    let new_h = adaptStep(f: f, y: y, x: x, h: h, relativeTol: relativeTol)
    let new_y = rungeKuttaFourthOrder(f: f, y: y, x: x, h: new_h)
    return new_y
    
    
    
}



public func adaptStep <T: AdaptiveSteppable> (f: (T, Double) -> T, y: T, x: Double, h: Double, relativeTol: Double) -> Double {
    
    let hh = h / 2.0
    let dh = h * 2.0
    let min_h = 0.0005
    
    let new_y_h = rungeKuttaFourthOrder(f: f, y: y, x: x, h: h)
    let new_y_hh = rungeKuttaFourthOrder(f: f, y: y, x: x, h: hh)
    let new_y_dh = rungeKuttaFourthOrder(f: f, y: y, x: x, h: dh)
    
    if h < min_h {
        return min_h
    }
 
  else if new_y_hh.checkEquivalenceUnderRelativeTolerance(new_y_h, relativeTol: relativeTol) == false {
        return adaptStep(f: f, y: y, x: x, h: hh, relativeTol: relativeTol)
    }
    
    else if new_y_h.checkEquivalenceUnderRelativeTolerance(new_y_dh, relativeTol: relativeTol) == false {
        return h
    }
    
    else {
        
        return adaptStep(f: f, y: y, x: x, h: dh, relativeTol: relativeTol)
    }
    
}

public func rungeKuttaFourthOrder <T: Integrable> (f: (T, Double) -> T, y: T, x: Double, h: Double) -> T {
    
    let K_1 = f(y , x)
    let K_2 = f(y + (h * K_1 * 0.5), x + h/2.0)
    let K_3 = f(y + (h * K_2 * 0.5), x + h/2.0)
    let K_4 = f(y + (h * K_3), x + h)
    
    let tot_slope = K_1 + 2.0 * K_2 + 2.0 * K_3 + K_4
    
    let total_change = (h/6)*tot_slope
    let new_y = y + total_change
    
    return new_y
    
}




public func doRungeKuttaStep
    < IndependentVariable: Has_IntegerInitializer &
                          Addable &
                          Dividable,
      IntegrandType > (
    t: IndependentVariable,
    h: IndependentVariable,
    y: IntegrandType,
    derivative_function derivs: (IndependentVariable,
                                 IntegrandType)
                                 -> IntegrandType,
    add: (IntegrandType, IntegrandType)
         -> IntegrandType,
    times: (IndependentVariable, IntegrandType)
           -> IntegrandType,
    dydx inputdydx: IntegrandType! = nil )
-> IntegrandType
{
    let dydx: IntegrandType
    if inputdydx != nil {
        dydx = inputdydx
    } else {
        dydx = derivs( t , y )
    }

    let hh = h / IndependentVariable(2)
    let h6 = h / IndependentVariable(6)
    let xh = t + hh
    var yt = add( y, times( hh, dydx ) )
    var dyt = derivs(xh,yt)
    yt = add( y, times( hh, dyt ) )
    var dym = derivs(xh, yt)
    yt = add( y, times( h, dym ) )
    dym = add( dym , dyt )
    dyt = derivs( t + h, yt )
    let sum1 = add( dydx , dyt )
    let sum2 = add( dym , dym )
    return add(y , times( h6, add( sum1, sum2 ) ) )
}

//public func doRungeKuttaStep <T: OdeIntegrable> (
//    t: Complex,
//    h: Complex,
//    y: T,
//    derivative_function return_derivatives:
//    (Complex, T) -> T,
//    dydx inputdydx: T! = nil
//) -> T
//where T.OdeScalar: Has_IntegerInitializer & Addable & Dividable
//{
//    if inputdydx != nil {
//        return doRungeKuttaStep(t: t,
//                                h: h,
//                                y: y,
//                                derivative_function: return_derivatives,
//                                add: T.odeAdd,
//                                times: T.odeMultiply,
//                                dydx: inputdydx)
//    } else {
//        return doRungeKuttaStep(t: t,
//                                h: h,
//                                y: y,
//                                derivative_function: return_derivatives,
//                                add: T.odeAdd,
//                                times: T.odeMultiply)
//    }
//}



// MARK: - Multistep method
// modified from odeint – does not store intermediate values (for simplified usage).
// TODO: Find out why this behaves as expected for Jyannes-Cummings but cannot integrate over long intervals for duffing oscillator (step size blows up)
//public func multiStepIvpIntegrator< T: OdeMultStepIntegrable >
//(
//    from x1: Double,
//    to x2: Double,
//    first_try_of_stepsize h1: Double,
//    smallest_allowed_value_of_stepsize hmin: Double,
//    accuracy eps: Double,
//    y ystart: inout T,
//    derivative_function derivs: (Double, T) -> T
//) {
//
//    let MAXSTP = 10000
//    let TINY = 1e-30
//    var x = x1
//    var h = h1
//    var y = ystart
//
//    var yscal = Array(repeating: 0.0, count: y.count)
//    let arrayTindicies = zip(yscal.indices, y.indices)
//
//    for _ in 0 ..< MAXSTP {
//        let dydx = derivs(x,y)
//        for (iscal, iT) in arrayTindicies {
//            yscal[iscal] = T.ode_abs(of: y[iT]) + ( T.ode_abs(of: dydx[iT]) * h ) + TINY
//        }
//        if ( (x + h - x2 )*( x + h - x1) > 0.0 ) {
//            h = x2  - x
//        }
//        var hdid = 0.0
//        var hnext = 0.0
//        rkqs(y: &y,
//             dydx: dydx,
//             x: &x,
//             htry: h,
//             eps: eps,
//             yscal: yscal,
//             hdid: &hdid,
//             hnext: &hnext,
//             derivative_function: derivs)
//        if ( (x - x2) * (x2 - x1) >= 0.0 ) {
//            ystart = y
//            return
//        }
//        if (hnext * hnext <= hmin * hmin) {
//            errorStream.write("Step size too small in odeint")
//        }
//        h=hnext;
//    }
//    errorStream.write("Too many steps in routine odeint")
//}
//
//// MARK: internal rk45 step with error estimate based of NRC rkqs
//func rkqs< T: OdeMultStepIntegrable >(y: inout T,
//                                      dydx: T,
//                                      x: inout Double,
//                                      htry: Double,
//                                      eps: Double,
//                                      yscal: [Double],
//                                      hdid: inout Double,
//                                      hnext: inout Double,
//                                      derivative_function derivs: (Double, T) -> T)
//{
//
//    let SAFETY = 0.9
//    let PGROW  = -0.2
//    let PSHRNK = -0.25
//    let ERRCON = 189.0/1_000_000
//    var ytemp = T.odeMultiply(scalar: 0.0, integrand: y)
//    var yerr = T.odeMultiply(scalar: 0.0, integrand: y)
//    var h = htry
//    let arrayTindicies = zip(yscal.indices, yerr.indices)
//
//    while (true) {
//        fifthOrderCashKarpRungeKutta(t: x,
//                                     h: h,
//                                     y: y,
//                                     yout: &ytemp,
//                                     yerr: &yerr,
//                                     derivative_function: derivs, dydx: dydx)
//
//        var errmax = 0.0
//        for (iscal, ierr) in arrayTindicies {
//            errmax = max(errmax , T.ode_abs(of: yerr[ierr]) / yscal[iscal] )
//        }
//        errmax = errmax / eps
//        if (errmax > 1.0) {
//            let htemp = SAFETY * h * pow(errmax,PSHRNK)
//            if (h >= 0.0 ) {
//                h = max(htemp, h / 10.0 )
//            } else {
//                h = min(htemp, h / 10.0 )
//            }
//            let xnew = (x) + h // only used here so bad name - maybe remove & change the if statement
//            if (xnew == x) {
//                errorStream.write("step-size underflow in rkqs")
//            }
//            continue
//        } else {
//            if (errmax > ERRCON) {
//                hnext = SAFETY * h * pow(errmax, PGROW)
//            } else {
//                hnext = 5.0 * h
//            }
//            hdid = h
//            x = x + hdid
//            y = ytemp
//            break
//        }
//    }
//}
//// MARK: - Fifth Order Cash Karp Runge Kutta
//public func fifthOrderCashKarpRungeKutta <T: OdeIntegrable> (
//    t: Double,
//    h: Double,
//    y: T,
//    yout: inout T,
//    yerr: inout T,
//    derivative_function return_derivatives:
//    (Double, T) -> T,
//    dydx inputdydx: T
//)
//
//{
//    fifthOrderCashKarpRungeKutta(t: t,
//                                 h: h,
//                                 y: y,
//                                 yout: &yout,
//                                 yerr: &yerr,
//                                 derivative_function: return_derivatives,
//                                 add: T.odeAdd,
//                                 times: T.odeMultiply,
//                                 dydx: inputdydx)
//}
//
//
//// based on rkck in numerical.recipes/book/book.html
//// Cash, J.R., and Karp, A.H. 1990, ACM Transactions on Mathematical Software, vol. 16, pp. 201– 222. [2]
//private func fifthOrderCashKarpRungeKutta
//    < IndependentVariable: Has_IntegerInitializer &
//                           Addable &
//                           Dividable &
//                           Multipliable,
//      IntegrandType>
//    (
//    t x: IndependentVariable,
//    h: IndependentVariable,
//    y: IntegrandType,
//    yout: inout IntegrandType,
//    yerr: inout IntegrandType,
//    derivative_function derivs: (IndependentVariable,
//                                 IntegrandType)
//                                 -> IntegrandType,
//    add: (IntegrandType, IntegrandType)
//         -> IntegrandType,
//    times: (IndependentVariable, IntegrandType)
//    -> IntegrandType ,
//    dydx: IntegrandType
//    )
//{
//    typealias R = IndependentVariable
//    // capitalised to avoid confusion with Addable sum in utility functions
//
//    func Sum(_ values: IntegrandType ...) -> IntegrandType {
//        return repeatedly (apply: add , values)
//    }
//
//    let a2 = R(2)/R(10),       a3  = R(3)/R(10),        a4 = R(6)/R(10),
//        a5 = R(1),             a6  = R(875)/R(1000),
//        b21 = R(2)/R(10),      b31 = R(3)/R(40),         b32 = R(9)/R(40),
//        b41 = R(3)/R(10),      b42 = R(-9)/R(10),        b43 = R(12)/R(10),
//        b51 = R(-11) / R(54),  b52 = R(25)/R(10),        b53 = R(-70)/R(27),
//        b54 = R(35)/R(27),     b61 = R(1631)/R(55296),   b62 = R(175)/R(512),
//        b63 = R(575)/R(13824), b64 = R(44275)/R(110592), b65 = R(253) / R(4096),
//        c1  = R(37)/R(378),    c3  = R(250)/R(621),      c4  = R(125)/R(594),
//        c6  = R(512)/R(1771),
//        dc5 = R(-277) / R(14336),
//        dc1 = c1 + ( R(-2825)  / R(27648)),
//        dc3 = c3 + ( R(-18575) / R(48384)),
//        dc4 = c4 + ( R(-13525) / R(55296) ),
//        dc6 = c6 + ( R( -25 )  / R(100) )
//
//    var ytemp = add (y , times (b21 * h , dydx) )  // First Step
//
//    let ak2 = derivs(x + a2 * h, ytemp) // Second Step
//
//    ytemp =  Sum( y ,
//                  times( h * b31, dydx) ,
//                  times( h * b32, ak2) )
//
//    let ak3 = derivs(x + a3 * h, ytemp) // Third Step
//
//    ytemp = Sum( y ,
//                 times( h*b41, dydx),
//                 times( h*b42, ak2) ,
//                 times( h*b43, ak3) )
//
//    let ak4 = derivs(x + a4 * h, ytemp) // Fourth Step
//
//    ytemp = Sum( y ,
//                 times( h*b51, dydx),
//                 times( h*b52, ak2) ,
//                 times( h*b53, ak3) ,
//                 times( h*b54, ak4) )
//
//    let ak5 = derivs(x + a5 * h, ytemp) // Fifth step.
//
//    ytemp = Sum( y ,
//                 times( h*b61, dydx) ,
//                 times( h*b62, ak2 ) ,
//                 times( h*b63, ak3 ) ,
//                 times( h*b64, ak4 ) ,
//                 times( h*b65, ak5 ) )
//    let ak6 = derivs(x+a6*h, ytemp)  // Sixth step.
//    // Accumulate increments with proper weights.
//    yout = Sum( y ,
//                times( h*c1, dydx) ,
//                times( h*c3, ak3 ) ,
//                times( h*c4, ak4 ) ,
//                times( h*c6, ak6 ))
//
//    //Estimate error as difference between fourth and fifth order methods.
//    yerr = Sum(times( h*dc1, dydx),
//               times( h*dc3, ak3 ),
//               times( h*dc4, ak4 ),
//               times( h*dc5, ak5 ),
//               times( h*dc6, ak6 ))
//}
//
//
//// Created by M J Everitt on 21/01/2022.
//
