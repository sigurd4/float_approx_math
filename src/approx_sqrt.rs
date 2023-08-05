use super::*;

#[const_trait]
pub trait ApproxSqrt
{
    /// Calculates an approximation of a square root.
    /// 
    /// The result gets iterated through the Newton-Raphson method with a given amount of iterations.
    /// 
    /// Even at 0 iterations, this algorithm is typically slower than the built-in sqrt function in the standard library,
    /// however this is possible to run at compile-time.
    /// 
    /// # Example
    /// 
    /// ```rust
    /// #![feature(const_trait_impl)]
    /// 
    /// use float_approx_math::ApproxSqrt;
    /// 
    /// const X: f32 = 2.0;
    /// const Y: f32 = X.approx_sqrt::<3>(); // Three iterations
    ///
    /// assert_eq!(Y, X.sqrt());
    /// ```
    /// 
    /// # Error
    /// 
    /// Calculating the square-root of 2:
    /// 
    /// ```rust
    /// use float_approx_math::ApproxSqrt;
    /// 
    /// const TWO: f64 = 2.0;
    /// 
    /// assert_eq!(TWO.sqrt(), 1.4142135623730951);
    /// 
    /// assert_eq!(TWO.approx_sqrt::<0>(), 1.5); // 6.07% error
    /// assert_eq!(TWO.approx_sqrt::<1>(), 1.4166666666666665); // 0.173% error
    /// assert_eq!(TWO.approx_sqrt::<2>(), 1.4142156862745097); // 0.000150% error
    /// assert_eq!(TWO.approx_sqrt::<3>(), 1.4142135623746899); // 0.000000000113% error
    /// ```
    fn approx_sqrt<const NEWTON: usize>(self) -> Self;
}

macro_rules! impl_approx_sqrt {
    ($float:ty: $bits:ty; $consts:tt) => {
        impl const ApproxSqrt for $float
        {
            fn approx_sqrt<const NEWTON: usize>(self) -> Self
            {
                let mut y = <$float>::from_bits(
                    (($consts::EXP_BIAS + 1) * (1 << (<$float>::MANTISSA_DIGITS as $bits - 2)))
                    + (<$float>::to_bits(self) >> 1)
                    - (1 << (<$float>::MANTISSA_DIGITS as $bits - 2))
                );
                let mut i = 0;
                while i < NEWTON
                {
                    y = (y + 2.0/y)/2.0;
                    i += 1;
                }
                y
            }
        }
    };
}
impl_approx_sqrt!(f32: u32; f32);
impl_approx_sqrt!(f64: u64; f64);