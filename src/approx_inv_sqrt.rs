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
    
    fn approx_inv_sqrt_unchecked<const NEWTON: usize>(self) -> Self;
}

macro_rules! impl_approx_inv_sqrt {
    ($float:ty: $bits:ty; $consts:tt) => {
        impl const ApproxInvSqrt for $float
        {
            fn approx_inv_sqrt<const NEWTON: usize>(self) -> Self
            {
                if self > 0.0
                {
                    self.approx_inv_sqrt_unchecked::<NEWTON>()
                }
                else
                {
                    <$float>::NAN
                }
            }
            fn approx_inv_sqrt_unchecked<const NEWTON: usize>(self) -> Self
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
mod test
{
    use ::test::Bencher;

    use super::*;
    use crate::tests as t;

    #[test]
    fn verify_magic_number_f32()
    {
        const L: u32 = 1 << (<f32>::MANTISSA_DIGITS - 1);
        const SIGMA: f64 = 0.0450466;
        const MAGIC_NUMBER: u32 = (1.5*L as f64*(f32::EXP_BIAS as f64 - SIGMA) + 0.5) as u32;
    
        assert_eq!(MAGIC_NUMBER, 0x5f3759df);
    }
    
    #[test]
    fn inv_sqrt()
    {
        const RANGE: f32 = 10.0;
        
        t::plot_approx("inv_sqrt_0", 0.1..RANGE, |x| x.sqrt().recip(), ApproxInvSqrt::approx_inv_sqrt::<0>);
        t::plot_approx("inv_sqrt_1", 0.1..RANGE, |x| x.sqrt().recip(), ApproxInvSqrt::approx_inv_sqrt::<1>);
        t::plot_approx("inv_sqrt_2", 0.1..RANGE, |x| x.sqrt().recip(), ApproxInvSqrt::approx_inv_sqrt::<2>);
        t::plot_approx("inv_sqrt_3", 0.1..RANGE, |x| x.sqrt().recip(), ApproxInvSqrt::approx_inv_sqrt::<3>);
    }
    
    #[bench]
    fn inv_sqrt_benchmark(_: &mut Bencher)
    {
        type F = f64;

        const N: usize = 1000;
        const S: usize = 32;

        t::plot_benchmark::<_, _, N, _>(
            "inv_sqrt",
            [
                &|x: F| x.sqrt().recip(),
                &ApproxInvSqrt::approx_inv_sqrt::<0>,
                &ApproxInvSqrt::approx_inv_sqrt::<1>,
                &ApproxInvSqrt::approx_inv_sqrt::<2>,
                &ApproxInvSqrt::approx_inv_sqrt::<3>,
                &ApproxInvSqrt::approx_inv_sqrt::<4>
            ],
            0.1..256.0,
            S
        );

        {
            const N: usize = 1000000;
            let x: Vec<F> = (0..N).map(|i| (i + 1) as F).collect();
    
            println!("sqrt dt = {:?}", t::benchmark(&x, &|x: F| x.sqrt().recip()));
            
            println!("approx_inv_sqrt::<0> dt = {:?}", t::benchmark(&x, &|x| x.approx_inv_sqrt::<0>()));
            println!("approx_inv_sqrt::<1> dt = {:?}", t::benchmark(&x, &|x| x.approx_inv_sqrt::<1>()));
            println!("approx_inv_sqrt::<2> dt = {:?}", t::benchmark(&x, &|x| x.approx_inv_sqrt::<2>()));
            println!("approx_inv_sqrt::<3> dt = {:?}", t::benchmark(&x, &|x| x.approx_inv_sqrt::<3>()));
            println!("approx_inv_sqrt::<4> dt = {:?}", t::benchmark(&x, &|x| x.approx_inv_sqrt::<4>()));
        }
    }
}

impl_approx_inv_sqrt!(f32: u32; f32);
impl_approx_inv_sqrt!(f64: u64; f64);