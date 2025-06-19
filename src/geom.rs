//! THis is where all geometry primitives live

use std::ops;

pub trait GeomCell<const N: usize, const M: usize> {
    type REFT;
    /// Reference cell
    fn refcell() -> Self::REFT;
    /// Jacobian measure to reference cell.
    fn jacobian_meas(self) -> f64;
    // /// Jacobian measure from arbitrary cell
    // fn jacobian_meas_arb(self, rf: Self) -> f64;
    /// Map to reference cell
    fn map_reference(self, p: Point<N>) -> Point<M>;
    // /// Map to arbitrary reference cell
    // fn map_reference_arb(self, p: Point<N>, rf: Self) -> Point<M>;
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point<const D: usize>([f64; D]);
pub type Triangle<const D: usize> = (Point<D>, Point<D>, Point<D>);

const REFERENCE_TRIANGLE: Triangle<2> = (Point([0.0, 0.0]), Point([1.0, 0.0]), Point([0.0, 1.0]));


impl<const D: usize> Point<D> {
    pub fn new(a: [f64; D]) -> Self {
        Point(a)
    }
}

impl <const D: usize> ops::Index<usize> for Point<D> {
    type Output = f64;

    fn index(&self, idx: usize) -> &f64 {
        &self.0[idx]
    }
}

impl <P: ops::Mul<f64, Output=f64> + Copy, const D: usize> ops::Mul<P> for Point<D> {
    type Output = Point<D>;

    fn mul(self, rhs: P) -> Self::Output {
        let mut p = Point([0.0; D]);

        for i in 0..D {
            p.0[i] = rhs*self.0[i];
        }
        p
    }
}

impl <const D: usize> ops::Sub<Point<D>> for Point<D> {
    type Output = Point<D>;

    fn sub(self, rhs: Point<D>) -> Self::Output {
        let mut p = Point([0.0; D]);

        for i in 0..D {
            p.0[i] =  self.0[i] - rhs.0[i];
        }
        p
    }
}

impl <const D: usize> ops::Add<Point<D>> for Point<D> {
    type Output = Point<D>;

    fn add(self, rhs: Point<D>) -> Self::Output {
        let mut p = Point([0.0; D]);

        for i in 0..D {
            p.0[i] = rhs.0[i] + self.0[i];
        }
        p
    }
}


impl<const D: usize> GeomCell<2, D> for Triangle<D> {
    /// Reference cell
    type REFT = Triangle<2>;

    fn refcell() -> Self::REFT {
        REFERENCE_TRIANGLE
    }
    /// Jacobian measure to reference cell.
    fn jacobian_meas(self) -> f64 {
        match D {
            0..=1 => panic!("This is not possible"),
            2 => 0.0,
            3 => 0.0,
            _ => panic!("Please use feature nönnönnöö"), // This should be translated into macro
            // tomfoolery
        }
    }
    /// Map to reference cell
    fn map_reference(self, p: Point<2>) -> Point<D> {
        self.0 + (self.1 - self.0)*p[0] + (self.1 - self.0)*p[1]
    }
}
