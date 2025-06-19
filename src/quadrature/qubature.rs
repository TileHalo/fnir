use itertools::izip;

use crate::geom::{GeomCell, Point, Triangle};


/// General N-dimensional quadrature (N > 1)
pub trait Qubature<A: GeomCell<M,  D>, I, O, const D: usize, const M: usize> {

    fn nint<F>(&self, func: F, cell: A) -> crate::Result<O> where F: Fn(Point<D>) -> O;

}

pub struct GaussTriQuadrature<const D: usize, const Q: usize> {
    pub a: [f64; D],
    pub b: [f64; D],
    pub c: [f64; D],
    pub xi: Vec<f64>,
    pub eta: Vec<f64>,
    pub nu: Vec<f64>,
}

pub struct GaussTetQuadrature<const D: usize> {
    pub a: [f64; D],
    pub b: [f64; D],
    pub c: [f64; D],
    pub xi: Vec<f64>,
    pub eta: Vec<f64>,
    pub nu: Vec<f64>,
}

impl <const D: usize, const Q: usize> Qubature<Triangle<D>, f64, f64, D, 2> for GaussTriQuadrature<D, Q> {

    fn nint<F>(&self, func: F, cell: Triangle<D>) -> crate::Result<f64> where F: Fn(Point<D>) -> f64 {
        let jac = cell.jacobian_meas();
        let res: f64 = izip!(self.xi.iter(), self.eta.iter(), self.nu.iter()).map(|(xi, eta, nu)| nu*func(cell.map_reference(Point::new([*xi, *eta])))).sum();
        Ok(jac*res)
    }
}

// These two should be implemented
// Function gausstriangle(Q):
//     Input:
//         Q – number of Gauss-Legendre points for the 1D base quadrature
//     Output:
//         xi, eta – coordinates of integration points in reference triangle
//         nu      – corresponding integration weights
//
//     1. Set P = Q  // Number of base 1D quadrature points
//
//     2. If P == 1:
//         - Use centroid of the triangle:
//             xi  = [1/3]
//             eta = [1/3]
//             nu  = [0.5]  // Area of the triangle
//
//     3. Else:
//         - Compute 1D Gauss-Legendre points and weights on [0,1]:
//             x, w1 = GaussLegendre(P)     // x: nodes, w1: weights
//             Transpose x to row vector
//
//         - Initialize:
//             - First point:
//                 xi_temp  = 1 - x[0]
//                 eta_temp = 0.5 * x[0]
//                 nu_temp  = x[0] * w1[0]
//
//         - For j = 1 to P-1 (i.e., second to last point in x):
//             a. Estimate number of vertical points:
//                 Pj = max(2, ceil(x[j] / x[P-1] * P))
//
//             b. Compute vertical Gauss-Legendre rule:
//                 yj, wj = GaussLegendre(Pj)
//
//             c. Scale yj and weights to triangle geometry:
//                 y_mapped = x[j] * yj      // eta direction
//                 w2 = x[j] * wj            // Jacobian scaled weights
//
//             d. For each y in y_mapped:
//                 xi_temp.append(1 - x[j])
//                 eta_temp.append(y)
//                 nu_temp.append(w1[j] * corresponding w2)
//
//         - Final outputs:
//             xi  = xi_temp
//             eta = eta_temp
//             nu  = nu_temp


// Function gausstetrahedron(Q):
//     Input:
//         Q – number of base Gauss-Legendre points (1D quadrature order)
//     Output:
//         xi, eta, zeta – integration points in the reference tetrahedron
//         nu            – corresponding integration weights
//
//     1. If Q == 1:
//         - Use centroid of the tetrahedron:
//             xi   = [1/4]
//             eta  = [1/4]
//             zeta = [1/4]
//             nu   = [1/6]  // Volume of reference tetrahedron
//
//     2. Else:
//         - Compute 1D Gauss-Legendre nodes and weights on [0,1]:
//             x, w1 = GaussLegendre(Q)
//
//         - Initialize empty lists:
//             xi_list, eta_list, zeta_list, weight_list = []
//
//         - Loop over first dimension (outer):
//             For i = 1 to Q:
//                 xi_i = x[i]
//                 wi   = w1[i]
//
//                 // Determine number of points in second dimension
//                 Qi = max(2, ceil(xi_i / x[Q-1] * Q))
//                 y, w2 = GaussLegendre(Qi)
//
//                 For j = 1 to Qi:
//                     eta_j = xi_i * y[j]
//                     wj = w2[j]
//
//                     // Determine number of points in third dimension
//                     Qij = max(2, ceil(eta_j / x[Q-1] * Q))
//                     z, w3 = GaussLegendre(Qij)
//
//                     For k = 1 to Qij:
//                         zeta_k = eta_j * z[k]
//                         wk = w3[k]
//
//                         // Map to reference tetrahedron:
//                         xi    = 1 - xi_i
//                         eta   = eta_j
//                         zeta  = zeta_k
//                         weight = wi * wj * wk * xi_i * eta_j  // Jacobian scaling
//
//                         Append xi, eta, zeta, weight to lists
//
//         - Return:
//             xi = xi_list
//             eta = eta_list
//             zeta = zeta_list
//             nu = weight_list
