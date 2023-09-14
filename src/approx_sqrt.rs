use super::*;

#[const_trait]
pub trait ApproxSqrt
{
    /// Calculates an approximation of a square root.
    /// 
    /// The result gets iterated through the Newton-Raphson method with a given amount of iterations.
    /// 
    /// Even at 0 iterations, this algorithm is typically a bit slower than the built-in sqrt function in the standard library,
    /// however this method can be run at compile-time.
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

#[cfg(test)]
mod test
{
    use ::test::Bencher;

    use super::*;
    use crate::tests as t;

    #[test]
    pub fn sqrt_error()
    {
        const X: f64 = 2.0;
        const Y: f64 = X.approx_sqrt::<0>();

        println!("{}", X.sqrt());
        println!("{}", Y);
        println!("error = {}", (Y - X.sqrt())/X.sqrt());
    }
    
    #[bench]
    fn sqrt_benchmark(_: &mut Bencher)
    {
        type F = f64;

        const N: usize = 1000;
        const S: usize = 32;

        t::plot_benchmark::<_, _, N, _>(
            "sqrt",
            [
                &F::sqrt,
                &ApproxSqrt::approx_sqrt::<0>,
                &ApproxSqrt::approx_sqrt::<1>,
                &ApproxSqrt::approx_sqrt::<2>,
                &ApproxSqrt::approx_sqrt::<3>,
                &ApproxSqrt::approx_sqrt::<4>
            ],
            0.0..256.0,
            S
        );

        {
            const N: usize = 1000000;
            let x: Vec<F> = (0..N).map(|i| (i + 1) as F).collect();
    
            println!("sqrt dt = {:?}", t::benchmark(&x, &F::sqrt));
            
            println!("approx_sqrt::<0> dt = {:?}", t::benchmark(&x, &|x| x.approx_sqrt::<0>()));
            println!("approx_sqrt::<1> dt = {:?}", t::benchmark(&x, &|x| x.approx_sqrt::<1>()));
            println!("approx_sqrt::<2> dt = {:?}", t::benchmark(&x, &|x| x.approx_sqrt::<2>()));
            println!("approx_sqrt::<3> dt = {:?}", t::benchmark(&x, &|x| x.approx_sqrt::<3>()));
            println!("approx_sqrt::<4> dt = {:?}", t::benchmark(&x, &|x| x.approx_sqrt::<4>()));
        }
    }
}