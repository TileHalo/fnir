//! Implements various abstract algebra abstractions.
//! Here, Commutative* are traits such that the output of operations is self, whereas Commuting*
//! the output is `T`.
//!
//! For constructs like magmas, monoids quasi-, semigroups, and groups for both additive and
//! multiplicative variants are present as `Additive*' and `Multiplicative*` together with their
//! commutative variants. Please note that commutative groups are named as `*AbelianGroup`
//!
//! Rings, Commutative rings and fields are implemented.
//!
//! Finally, vector spaces and modules are present.
use std::ops;

pub mod polynomial;

// Commuting lhs
/// Commuting add is helper trait for type `U` which is the result of commutative addition with `T`.
/// This type abstracts over when LHS is `T`
pub trait CommutingAdd<T: ops::Add<Self, Output = T>>:
    ops::Add<T, Output = T> + std::marker::Sized
{
}
/// CommutingSub is helper trait for type `U` which is the result of commutative subtraction with `T`.
/// This type abstracts over when LHS is `T`
pub trait CommutingSub<T: ops::Sub<Self, Output = T>>:
    ops::Sub<T, Output = T> + std::marker::Sized
{
}
/// [`CommutingMul`] is helper trait for type `U` which is the result of commutative multiplication with `T`.
/// This type abstracts over when LHS is `T`
pub trait CommutingMul<T: ops::Mul<Self, Output = T>>:
    ops::Mul<T, Output = T> + std::marker::Sized
{
}
/// [`CommutingDiv`] is helper trait for type `U` which is the result of commutative division with `T`.
/// This type abstracts over when LHS is `T`
pub trait CommutingDiv<T: ops::Div<Self, Output = T>>:
    ops::Div<T, Output = T> + std::marker::Sized
{
}

// Commutative stuff
/// Commutative addition
pub trait CommutativeAdd<T: CommutingAdd<Self>>:
    ops::Add<T, Output = Self> + ops::AddAssign<T> + std::marker::Sized
{
}
/// Commutative subtraction
pub trait CommutativeSub<T: CommutingSub<Self>>:
    ops::Sub<T, Output = Self> + ops::SubAssign<T> + std::marker::Sized
{
}
/// Commutative multiplication
pub trait CommutativeMul<T: CommutingMul<Self>>:
    ops::Mul<T, Output = Self> + ops::MulAssign<T> + std::marker::Sized
{
}

/// Commutative division
pub trait CommutativeDiv<T: CommutingDiv<Self>>:
    ops::Div<T, Output = Self> + ops::DivAssign<T> + std::marker::Sized
{
}

// Commuting trait gatherings
/// Helper trait when both addition and subtraction commute
pub trait CommutingAddSub<T: ops::Add<Self, Output = T> + ops::Sub<Self, Output = T>>:
    CommutingAdd<T> + CommutingSub<T>
{
}

/// Helper trait when both addition and multiplication commute
pub trait CommutingAddMul<T: ops::Add<Self, Output = T> + ops::Mul<Self, Output = T>>:
    CommutingAdd<T> + CommutingMul<T>
{
}

/// Helper trait when addition, subtraction and multiplication commute
pub trait CommutingRing<
    T: ops::Add<Self, Output = T> + ops::Mul<Self, Output = T> + ops::Sub<Self, Output = T>,
>: CommutingAdd<T> + CommutingMul<T> + CommutingSub<T>
{
}

/// Helper trait when all operaotions commute
pub trait CommutingAll<
    T: ops::Add<Self, Output = T>
        + ops::Mul<Self, Output = T>
        + ops::Sub<Self, Output = T>
        + ops::Div<Self, Output = T>,
>: CommutingAdd<T> + CommutingMul<T> + CommutingSub<T> + CommutingDiv<T>
{
}

// Implement commuting LHS stuff
impl<U: ops::Add<T, Output = T> + std::marker::Sized, T: ops::Add<Self, Output = T>> CommutingAdd<T>
    for U
{
}
impl<U: ops::Sub<T, Output = T> + std::marker::Sized, T: ops::Sub<Self, Output = T>> CommutingSub<T>
    for U
{
}
impl<U: ops::Mul<T, Output = T> + std::marker::Sized, T: ops::Mul<Self, Output = T>> CommutingMul<T>
    for U
{
}
impl<U: ops::Div<T, Output = T> + std::marker::Sized, T: ops::Div<Self, Output = T>> CommutingDiv<T>
    for U
{
}

// Implement commuting RHS
impl<P: CommutingAdd<T>, T: ops::Add<P, Output = T> + ops::AddAssign<P>> CommutativeAdd<P> for T {}
impl<P: CommutingSub<T>, T: ops::Sub<P, Output = T> + ops::SubAssign<P>> CommutativeSub<P> for T {}
impl<P: CommutingMul<T>, T: ops::Mul<P, Output = T> + ops::MulAssign<P>> CommutativeMul<P> for T {}
impl<P: CommutingDiv<T>, T: ops::Div<P, Output = T> + ops::DivAssign<P>> CommutativeDiv<P> for T {}

// Implement helpers
impl<
        U: CommutingAdd<T> + CommutingSub<T>,
        T: ops::Add<Self, Output = T> + ops::Sub<Self, Output = T>,
    > CommutingAddSub<T> for U
{
}

impl<
        U: CommutingAdd<T> + CommutingMul<T>,
        T: ops::Add<Self, Output = T> + ops::Mul<Self, Output = T>,
    > CommutingAddMul<T> for U
{
}

impl<
        U: CommutingAdd<T> + CommutingMul<T> + CommutingSub<T>,
        T: ops::Add<Self, Output = T> + ops::Mul<Self, Output = T> + ops::Sub<Self, Output = T>,
    > CommutingRing<T> for U
{
}

impl<
        U: CommutingAdd<T> + CommutingMul<T> + CommutingSub<T> + CommutingDiv<T>,
        T: ops::Add<Self, Output = T>
            + ops::Mul<Self, Output = T>
            + ops::Sub<Self, Output = T>
            + ops::Div<Self, Output = T>,
    > CommutingAll<T> for U
{
}

// Various additive abstract algebra constructs
pub trait AdditiveMagma<T>: ops::Add<T, Output = Self> + ops::AddAssign<T> {}
pub trait AdditiveCommutativeMagma<T: CommutingAdd<Self>>:
    AdditiveMagma<T> + CommutativeAdd<T>
{
}

pub trait AdditiveQuasigroup<T>:
    AdditiveMagma<T> + ops::Sub<T, Output = Self> + ops::SubAssign<T> + ops::Neg<Output = Self>
{
}
pub trait AdditiveCommutativeQuasigroup<T: CommutingAddSub<Self>>:
    AdditiveQuasigroup<T> + AdditiveCommutativeMagma<T> + CommutativeSub<T>
{
}

pub trait AdditiveSemigroup<T>: AdditiveMagma<T> {}
pub trait AdditiveCommutativeSemigroup<T: CommutingAdd<Self>>:
    AdditiveSemigroup<T> + AdditiveCommutativeMagma<T>
{
}

pub trait AdditiveMonoid<T>: AdditiveMagma<T> + num::Zero {}
pub trait AdditiveCommutativeMonoid<T: CommutingAdd<Self>>:
    AdditiveMonoid<T> + AdditiveCommutativeMagma<T> + num::Zero
{
}

pub trait AdditiveGroup<T>: AdditiveMonoid<T> + AdditiveQuasigroup<T> {}
pub trait AdditiveAbelianGroup<T: CommutingAddSub<Self>>:
    AdditiveGroup<T> + AdditiveCommutativeMonoid<T> + AdditiveCommutativeQuasigroup<T>
{
}

// Impl stuff
impl<T: ops::Add<P, Output = Self> + ops::AddAssign<P>, P> AdditiveMagma<P> for T {}
impl<T: AdditiveMagma<P>, P: CommutingAdd<T>> AdditiveCommutativeMagma<P> for T {}

impl<
        T: AdditiveMagma<P> + ops::Sub<P, Output = Self> + ops::SubAssign<P> + ops::Neg<Output = T>,
        P,
    > AdditiveQuasigroup<P> for T
{
}
impl<T: AdditiveQuasigroup<P>, P: CommutingAddSub<T>> AdditiveCommutativeQuasigroup<P> for T {}

impl<T: AdditiveMagma<P>, P> AdditiveSemigroup<P> for T {}
impl<T: AdditiveSemigroup<P>, P: CommutingAdd<T>> AdditiveCommutativeSemigroup<P> for T {}

impl<T: AdditiveSemigroup<P> + num::Zero, P> AdditiveMonoid<P> for T {}
impl<T: AdditiveMonoid<P>, P: CommutingAdd<T>> AdditiveCommutativeMonoid<P> for T {}

impl<
        T: AdditiveMonoid<P> + ops::Sub<P, Output = T> + ops::SubAssign<P> + ops::Neg<Output = T>,
        P,
    > AdditiveGroup<P> for T
{
}
impl<T: AdditiveGroup<P>, P: CommutingAddSub<T>> AdditiveAbelianGroup<P> for T {}

// Various multiplicative abstract algebra constructs
pub trait MultiplicativeMagma<T>: ops::Mul<T, Output = Self> + ops::MulAssign<T> {}
pub trait MultiplicativeCommutativeMagma<T: CommutingMul<Self>>:
    MultiplicativeMagma<T> + CommutativeMul<T>
{
}

pub trait MultiplicativeQuasigroup<T>:
    MultiplicativeMagma<T> + ops::Div<T, Output = Self> + ops::DivAssign<T>
{
}
pub trait MultiplicativeCommutativeQuasigroup<T: CommutingMul<Self> + CommutingDiv<Self>>:
    MultiplicativeQuasigroup<T> + MultiplicativeCommutativeMagma<T> + CommutativeDiv<T>
{
}

pub trait MultiplicativeSemigroup<T>: MultiplicativeMagma<T> {}
pub trait MultiplicativeCommutativeSemigroup<T: CommutingMul<Self>>:
    MultiplicativeSemigroup<T> + MultiplicativeCommutativeMagma<T>
{
}

pub trait MultiplicativeMonoid<T>: MultiplicativeMagma<T> + num::One {}
pub trait MultiplicativeCommutativeMonoid<T: CommutingMul<Self>>:
    MultiplicativeMonoid<T> + MultiplicativeCommutativeMagma<T> + num::One
{
}

pub trait MultiplicativeGroup<T>: MultiplicativeMonoid<T> + MultiplicativeQuasigroup<T> {}
pub trait MultiplicativeAbelianGroup<T: CommutingMul<Self> + CommutingDiv<Self>>:
    MultiplicativeGroup<T> + MultiplicativeCommutativeMonoid<T> + MultiplicativeCommutativeQuasigroup<T>
{
}

impl<T: ops::Mul<P, Output = Self> + ops::MulAssign<P>, P> MultiplicativeMagma<P> for T {}
impl<T: MultiplicativeMagma<P>, P: CommutingMul<T>> MultiplicativeCommutativeMagma<P> for T {}

impl<T: MultiplicativeMagma<P> + ops::Div<P, Output = Self> + ops::DivAssign<P>, P>
    MultiplicativeQuasigroup<P> for T
{
}
impl<T: MultiplicativeQuasigroup<P>, P: CommutingMul<T> + CommutingDiv<T>>
    MultiplicativeCommutativeQuasigroup<P> for T
{
}

impl<T: MultiplicativeMagma<P>, P> MultiplicativeSemigroup<P> for T {}
impl<T: MultiplicativeSemigroup<P>, P: CommutingMul<T>> MultiplicativeCommutativeSemigroup<P>
    for T
{
}

impl<T: MultiplicativeSemigroup<P> + num::One, P> MultiplicativeMonoid<P> for T {}
impl<T: MultiplicativeMonoid<P>, P: CommutingMul<T>> MultiplicativeCommutativeMonoid<P> for T {}

impl<T: MultiplicativeMonoid<P> + ops::Div<P, Output = T> + ops::DivAssign<P>, P>
    MultiplicativeGroup<P> for T
{
}
impl<T: MultiplicativeGroup<P>, P: CommutingMul<T> + CommutingDiv<T>> MultiplicativeAbelianGroup<P>
    for T
{
}

// Ring and field

pub trait Ring<T: CommutingSub<Self> + CommutingAdd<Self>>:
    AdditiveAbelianGroup<T> + MultiplicativeMonoid<T>
{
}
pub trait CommutativeRing<T: CommutingSub<Self> + CommutingAdd<Self> + CommutingMul<Self>>:
    Ring<T> + MultiplicativeCommutativeMonoid<T>
{
}

pub trait Field<T: CommutingSub<Self> + CommutingAdd<Self> + CommutingMul<Self>>:
    CommutativeRing<T> + MultiplicativeQuasigroup<T>
{
}

impl<
        T: AdditiveAbelianGroup<P> + MultiplicativeMonoid<P>,
        P: CommutingSub<T> + CommutingAdd<T>,
    > Ring<P> for T
{
}
impl<
        T: Ring<P> + MultiplicativeCommutativeMonoid<P>,
        P: CommutingSub<T> + CommutingAdd<T> + CommutingMul<T>,
    > CommutativeRing<P> for T
{
}
impl<
        T: CommutativeRing<P> + MultiplicativeQuasigroup<P>,
        P: CommutingSub<T> + CommutingAdd<T> + CommutingMul<T> + CommutingDiv<T>,
    > Field<P> for T
{
}

// Module and vector space

/// This is Right $R$-module $M$.
/// For further information check [`Module`] and [`VectorSpace`]
pub trait RightModule<
    R: Ring<T> + ops::Mul<Self, Output = Self>,
    T: CommutingAddSub<R> + ops::Mul<Self, Output = Self>,
    P: CommutingAddSub<Self>,
>: AdditiveAbelianGroup<P>
{
}
/// This is Left $R$-module $M$.
/// For further information check [`Module`] and [`VectorSpace`]
pub trait LeftModule<
    R: Ring<T>,
    T: CommutingAddSub<R>,
    P: CommutingAddSub<Self>,
>: AdditiveAbelianGroup<P> + ops::Mul<R, Output = Self>  + ops::Mul<T, Output = Self>
{
}

/// This is a $R$-module $M$, when $R$ is a commutative ring and thus both left and right modules
/// are the same. This is generalization of [`VectorSpace`] with rings instead of fields.
/// For more information check [`VectorSpace`]
pub trait Module<
    R: CommutativeRing<T> + CommutingMul<Self>,
    T: CommutingRing<R> + CommutingMul<Self>,
    P: CommutingAddSub<Self>,
>: AdditiveAbelianGroup<P> + CommutativeMul<R> + CommutativeMul<T>
{
}

/// Vector space where $R$ is a field
pub trait VectorSpace<
    R: Field<T> + Field<R> + CommutingMul<Self>,
    T: CommutingAll<R> + CommutingMul<Self>,
    P: CommutingAddSub<Self>,
>: AdditiveAbelianGroup<P> + CommutativeMul<R> + CommutativeMul<T>
{
}

impl<
        M: AdditiveAbelianGroup<P>,
        R: Ring<T> + ops::Mul<Self, Output = Self>,
        T: CommutingAddSub<R> + ops::Mul<Self, Output = Self>,
        P: CommutingAddSub<Self>,
    > RightModule<R, T, P> for M
{
}
impl<
        M: AdditiveAbelianGroup<P> + ops::Mul<R, Output = Self> + ops::Mul<T, Output = M>,
        R: Ring<T>,
        T: CommutingAddSub<R> + CommutingMul<Self>,
        P: CommutingAddSub<Self>,
    > LeftModule<R, T, P> for M
{
}
impl<
        M: AdditiveAbelianGroup<P> + CommutativeMul<R> + CommutativeMul<T>,
        R: CommutativeRing<T> + CommutingMul<M>,
        T: CommutingRing<R> + CommutingMul<M>,
        P: CommutingAddSub<Self>> Module<R, T, P> for M
{
}

impl<
        M: AdditiveAbelianGroup<P>  + CommutativeMul<R> + CommutativeMul<T>,
        R: Field<T> + Field<R> + CommutingMul<Self>,
        T: CommutingAll<R> + CommutingMul<Self>,
        P: CommutingAddSub<Self>,
    > VectorSpace<R, T, P> for M
{
}

/// Float helper, that is, types marked with this trait behave like floats and can be intermixed
/// with floats.
pub trait FloatField:
    num_complex::ComplexFloat + Field<Self> + Field<f64> + CommutingAll<Self> + PartialOrd + PartialEq
where
    f64: CommutingAll<Self>,
{
}

impl<
        U: num_complex::ComplexFloat
            + Field<Self>
            + Field<f64>
            + CommutingAll<Self>
            + PartialOrd
            + PartialEq,
    > FloatField for U
where
    f64: CommutingAll<Self>,
{
}
