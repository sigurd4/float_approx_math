use super::*;

#[const_trait]
pub trait ApproxInvSqrt
{
    /// Calculates an approximation of 1/sqrt(x).
    /// 
    /// The result gets iterated through the Newton-Raphson method with a given amount of iterations.
    /// 
    /// At 5 or less iterations, it is generally faster than x.sqrt().recip() in release build.
    /// 
    /// 1 iteration was good enough for quake.
    /// 
    /// # Example
    /// 
    /// ```rust
    /// #![feature(const_trait_impl)]
    /// 
    /// use float_approx_math::ApproxInvSqrt;
    /// 
    /// const X: f32 = 2.0;
    /// const Y: f32 = X.approx_inv_sqrt::<4>(); // Three iterations
    ///
    /// assert_eq!(Y, X.sqrt().recip());
    /// ```
    fn approx_inv_sqrt<const NEWTON: usize>(self) -> Self;
}

macro_rules! impl_approx_inv_sqrt {
    ($float:ty: $bits:ty; $consts:tt) => {
        impl const ApproxInvSqrt for $float
        {
            fn approx_inv_sqrt<const NEWTON: usize>(self) -> Self
            {
                const L: $bits = 1 << (<$float>::MANTISSA_DIGITS as $bits - 1);
                const SIGMA: f64 = 0.0450466;
                const MAGIC_NUMBER: $bits = (1.5*L as f64*($consts::EXP_BIAS as f64 - SIGMA) + 0.5) as $bits;

                let mut y = <$float>::from_bits(MAGIC_NUMBER - (<$float>::to_bits(self) >> 1));

                let x2 = self*0.5;
                let mut i = 0;
                while i < NEWTON
                {
                    y *= 1.5 - x2*y*y;
                    i += 1;
                }
                y
            }
        }
    };
}

#[cfg(test)]
#[test]
fn verify_magic_number_f32()
{
    const L: u32 = 1 << (<f32>::MANTISSA_DIGITS - 1);
    const SIGMA: f64 = 0.0450466;
    const MAGIC_NUMBER: u32 = (1.5*L as f64*(f32::EXP_BIAS as f64 - SIGMA) + 0.5) as u32;

    assert_eq!(MAGIC_NUMBER, 0x5f3759df);
}

impl_approx_inv_sqrt!(f32: u32; f32);
impl_approx_inv_sqrt!(f64: u64; f64);