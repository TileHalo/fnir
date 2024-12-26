//! Numerical quadrature rules.
//! To use a (generic, non-existent) rule `Rule`, one must first call `Rule::new()`,
//! and then either use `Rule::nint` or `Rule::integrate` depending if rule modifications are
//! desired or not.
use std::ops;

use crate::algebra as alg;

/// Trait that all quadrature rules implement. All quadratures blanket implement the (`Integral`)[`super::Integral`]
/// trait
pub trait Quadrature<I, O> {
    /// Default number of iterations for quadrature
    const DEFAULTN: usize;
    fn nint<F>(&self, func: F, a: Option<I>, b: Option<I>, n: usize) -> crate::Result<O>
    where
        F: Fn(I) -> O;
}

impl<Q: Quadrature<I, O>, I, O> super::Integral<I, O> for Q {
    fn integrate<F>(&self, func: F, a: Option<I>, b: Option<I>) -> crate::Result<O>
    where
        F: Fn(I) -> O,
    {
        self.nint(func, a, b, Self::DEFAULTN)
    }
}

/// Standard trapezoidal rule for finite intervals.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Trapezoidal(usize);
/// Michalski-Mosig version of tanh-sinh/exp-sinh/sinh-sinh double ended quadrature with improvements
/// from Dr. Robert  van Engelen from Genivia Labs [qtsh](https://www.genivia.com/files/qthsh.pdf)
/// This quadrature accepts finite, semi-infinite, and infinite intervals, and uses
/// tanh-sinh, exp-sinh, and sinh-sinh rules respectively
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct DEQuad {
    eta: usize,
    nmax: usize,
    kappa: f64,
    maxlev: usize,
    eps: f64,
}

impl<
        I: alg::FloatField + ops::Mul<O, Output = O>,
        O: alg::VectorSpace<I, f64, O> + alg::CommutingAddSub<O> + ops::Div<f64, Output = O>,
    > Quadrature<I, O> for Trapezoidal
where
    f64: alg::CommutingAll<I> + alg::CommutingMul<O>,
{
    const DEFAULTN: usize = 20;
    fn nint<F>(&self, func: F, start: Option<I>, end: Option<I>, n: usize) -> crate::Result<O>
    where
        F: Fn(I) -> O,
    {
        match (start, end) {
            (Some(a), Some(b)) => {
                let l: I = b - a;
                let dx: I = l / n as f64;
                let mut res = O::zero();
                for i in 0..n {
                    res += (dx
                        * ((func(a + (i as f64) * dx)) + func(a + ((i as f64) + I::one()) * dx)))
                        / 2.0;
                }
                Ok(res)
            }
            (None, Some(_)) | (Some(_), None) => Err(crate::IntegrationError::SemiInfiniteIntegral),
            (None, None) => Err(crate::IntegrationError::InfiniteIntegral),
        }
    }
}

impl DEQuad {
    pub fn new() -> Self {
        Self {
            eta: 1,
            nmax: 24,
            kappa: 1.0e-15,
            maxlev: 5,
            eps: 1.0e-7,
        }
    }
    fn tanhsinh<
        O: num_complex::ComplexFloat<Real = I> + alg::VectorSpace<I, f64, O>,
        I: alg::FloatField + ops::Mul<O, Output = O>,
        F,
    >(
        &self,
        func: F,
        a: I,
        b: I,
        nmax: usize,
    ) -> crate::Result<O>
    where
        F: Fn(I) -> O,
        f64: alg::CommutingAll<I> + alg::CommutingMul<O>,
    {
        let nm: usize = if nmax == 0 { 24 } else { nmax };
        let maxlev = if self.maxlev == 0 { 6 } else { self.maxlev };
        let c = (a + b) / 2.0;
        let d = (b - a) / 2.0;
        let mut h = 2.0;

        let mut s = func(c);

        let (mut fp, mut fm) = (O::zero(), O::zero());

        let eps: f64 = if self.eps == 0.0 { 1.0e-9 } else { self.eps };
        let tol: f64 = 10.0 * eps;

        for k in 0..maxlev {
            let mut p = O::zero();
            h /= 2.0;
            let mut eh = f64::exp(h);
            let mut t = eh * I::one();
            if k > 0 {
                eh *= eh;
            }
            for _ in 0..nm {
                let u = I::exp(1.0 / t - t);
                let r = 2.0 * u / (1.0 + u);
                let w = (t + 1.0 / t) * r / (1.0 + u);
                let x = d * r;
                if a + x > a {
                    let y = func(a + x);
                    if y.is_finite() {
                        fp = y;
                    }
                }
                if b + x > b {
                    let y = func(b - x);
                    if y.is_finite() {
                        fm = y;
                    }
                }
                let q = w * (fp + fm);
                p += q;
                t *= eh;
                if q.abs() <= eps * p.abs() {
                    break;
                }
            }

            let v = s - p;
            s += p;
            if v.abs() <= (tol * s).abs() {
                break;
            }
        }

        Ok(d * h * s)
    }
    fn expsinh<
        O: num_complex::ComplexFloat<Real = I> + alg::VectorSpace<I, f64, O>,
        I: alg::FloatField + ops::Mul<O, Output = O>,
        F,
    >(
        &self,
        func: F,
        a: I,
        nmax: usize,
    ) -> crate::Result<O>
    where
        F: Fn(I) -> O,
        f64: alg::CommutingAll<I> + alg::CommutingMul<O>,
    {
        let nm: usize = if nmax == 0 { 24 } else { nmax };
        let maxlev = if self.maxlev == 0 { 6 } else { self.maxlev };
        let c = a;
        let d = 1.0;
        let mut h = 2.0;

        let mut s = func(a + d);

        let eps: f64 = if self.eps == 0.0 { 1.0e-9 } else { self.eps };
        let tol: f64 = 10.0 * eps;

        for k in 0..maxlev {
            let mut q = O::zero();
            let mut p = q;
            h /= 2.0;
            let mut eh = f64::exp(h);
            let mut t = eh * 0.5 * I::one();
            if k > 0 {
                eh *= eh;
            }
            for _ in 0..nm {
                q = O::zero();
                let r = I::exp(t - 0.25 / t);
                let w = r;
                let (x1, x2) = (c + d / r, c + d * r);
                if x1 == c { break; }
                let (y1, y2) = (func(x1), func(x2));
                if y1.is_finite() {
                    q += y1 * (1.0 / w);
                }
                if y2.is_finite() {
                    q += y2 *  w;
                }
                q *= t + 0.25 / t;
                p += q;
                t *= eh;
                if q.abs() <= eps * p.abs() {
                    break;
                }
            }

            let v = s - p;
            s += p;
            if v.abs() <= (tol * s).abs() {
                break;
            }
        }

        Ok(d * h * s)
    }

    fn sinhsinh<
        O: num_complex::ComplexFloat<Real = I> + alg::VectorSpace<I, f64, O>,
        I: alg::FloatField + ops::Mul<O, Output = O>,
        F,
    >(
        &self,
        func: F,
        nmax: usize,
    ) -> crate::Result<O>
    where
        F: Fn(I) -> O,
        f64: alg::CommutingAll<I> + alg::CommutingMul<O>,
    {
        let nm: usize = if nmax == 0 { 24 } else { nmax };
        let maxlev = if self.maxlev == 0 { 6 } else { self.maxlev };
        let c = 0.0;
        let d = 1.0;
        let mut h = 2.0;

        let mut s = func(I::zero());

        let eps: f64 = if self.eps == 0.0 { 1.0e-9 } else { self.eps };
        let tol: f64 = 10.0 * eps;

        for k in 0..maxlev {
            let mut q = O::zero();
            let mut p = q;
            h /= 2.0;
            let mut eh = f64::exp(h);
            let mut t = eh * 0.5 * I::one();
            if k > 0 {
                eh *= eh;
            }
            for _ in 0..nm {
                q = O::zero();
                let r = (I::exp(t - 0.25 / t) - 1.0 / I::exp(t - 0.25 / t)) / 2.0;
                let w = (I::exp(t - 0.25 / t) + 1.0 / I::exp(t - 0.25 / t)) / 2.0;
                let (x1, x2) = (c - d * r, c + d * r);
                let (y1, y2): (O, O) = (func(x1), func(x2));
                if y1.is_finite() {
                    q += y1 * w;
                }
                if y2.is_finite() {
                    q += y2 * w;
                }
                q *= t + 0.25 / t;
                p += q;
                t *= eh;
                if q.abs() <= eps * p.abs() {
                    break;
                }
            }

            let v = s - p;
            s += p;
            if v.abs() <= (tol * s).abs() {
                break;
            }
        }

        Ok(d * h * s)
    }
}

impl Default for DEQuad {
    fn default() -> Self {
        Self::new()
    }
}

impl<
        O: num_complex::ComplexFloat<Real = I> + alg::VectorSpace<I, f64, O>,
        I: alg::FloatField + ops::Mul<O, Output = O>,
    > Quadrature<I, O> for DEQuad
where
    f64: alg::CommutingAll<I> + alg::CommutingMul<O>,
{
    const DEFAULTN: usize = 24;
    fn nint<F>(&self, func: F, start: Option<I>, end: Option<I>, nmax: usize) -> crate::Result<O>
    where
        F: Fn(I) -> O,
    {
        match (start, end) {
            (Some(a), Some(b)) => self.tanhsinh(func, a, b, nmax),
            (Some(a), None) => self.expsinh(func, a, nmax),
            (None, Some(a)) => Ok(-self.expsinh(func, a, nmax)?),
            (None, None) => self.sinhsinh(func, nmax)
        }
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;

    use super::*;

    #[test]
    fn tanhsinh_test() {
        let rule = DEQuad::new();
        // Test various integrals from 0 to 1
        let ints: &[fn(f64) -> f64] = &[
            |x: f64| 1.0 / (1.0 + x),
            |x| 4.0 / (1.0 + x * x),
            |x| f64::acos(x),
            |x| f64::sin(x) / x,
            |x| f64::sqrt(x / (1.0 - x.powi(2))),
            |x| 1.0 / f64::sqrt(x),
            |x| 1.0 / f64::sqrt(1.0 - x),
            |x| f64::powf(x, -0.8),
            |x| 1.0 / f64::sqrt(f64::sin(std::f64::consts::PI * x)),
            |x| 1.0 / f64::sqrt(-f64::log10(x)),
        ];

        let cints: &[fn(f64) -> f64] = &[
            |x| 4.0 / (1.0 + x * x),
        ];

        let resc: &[f64] = &[2.0*f64::atan(33.0/56.0)];

        let ress: &[f64] = &[
            f64::ln(2.0),
            std::f64::consts::PI,
            1.0,
            0.946083,
            1.19814,
            2.0,
            2.0,
            5.0,
            1.6692537,
            f64::sqrt(std::f64::consts::PI * std::f64::consts::LN_10),
        ];
        for (i1, res) in ints.iter().zip(ress.iter()) {
            let ans = rule.nint(i1, Some(0.0), Some(1.0), 24).ok().unwrap();
            assert_approx_eq!(ans, res);
        }
        for (i1, res) in cints.iter().zip(resc.iter()) {
            let ans = rule.nint(i1, Some(2.0), Some(5.0), 24).ok().unwrap();
            assert_approx_eq!(ans, res);
        }
    }

    #[test]
    fn expsinh_test() {
        let rule = DEQuad::new();
        // Test various integrals from 0 to 1
        let ints: &[fn(f64) -> f64] = &[
            |x| 1.0/(1.0 + x*x),
            |x| f64::exp(-x)/f64::sqrt(x),
            |x| f64::exp(-1.0-x)/(1.0 + x),
            |x| (x*x)*f64::exp(-4.0*x),
            |x| f64::powi(f64::sqrt((x*x) + 9.0) - x, 3),
            |x| f64::exp(-3.0*x),
        ];

        let cints: &[fn(f64) -> f64] = &[
            |x| f64::exp(-x)/x,
        ];


        let ress: &[f64] = &[
            std::f64::consts::FRAC_PI_2,
            1.772453851,
            0.219383934,
            0.03125,
            30.375,
            0.333333333,
        ];

        let resc: &[f64] = &[0.219383934];
        for (i1, res) in ints.iter().zip(ress.iter()) {
            let ans = rule.nint(i1, Some(0.0), None, 24).ok().unwrap();
            assert_approx_eq!(ans, res);
        }
        for (i1, res) in cints.iter().zip(resc.iter()) {
            let ans = rule.nint(i1, Some(1.0), None, 24).ok().unwrap();
            assert_approx_eq!(ans, res);
        }
    }
    #[test]
    fn sinhsinh_test() {
        let rule = DEQuad::new();
        // Test various integrals from 0 to 1
        let ints: &[fn(f64) -> f64] = &[
            |x| x,
            |x| f64::exp(-(x.powi(2) + 2.0 * x + 1.0))
        ];



        let ress: &[f64] = &[
            0.0,
            f64::sqrt(std::f64::consts::PI),
        ];

        for (i1, res) in ints.iter().zip(ress.iter()) {
            let ans = rule.nint(i1, None, None, 24).ok().unwrap();
            assert_approx_eq!(ans, res);
        }

    }
}
