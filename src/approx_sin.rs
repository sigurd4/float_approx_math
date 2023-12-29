use super::*;

use polynomial_ops::*;
use array__ops::*;

#[const_trait]
pub trait ApproxSin
{
    /// Calculates an approximation of a sine, using the Z80 / ZX Spectrum algorithm.
    /// 
    /// https://namoseley.wordpress.com/2012/09/26/arduinoavr-and-zx-spectrum-sin-routines/
    /// 
    /// This approximation is significantly slower than the built-in sin method in the standard library,
    /// however this method can be run at compile-time.
    /// 
    /// # Example
    /// 
    /// ```rust
    /// #![feature(const_trait_impl)]
    /// 
    /// use float_approx_math::ApproxSin;
    /// 
    /// const X: f32 = 2.0;
    /// let y: f32 = X.approx_sin();
    ///
    /// assert!((y - X.sin()).abs() < 0.0000005); // Less than 0.0000005 abs error
    /// ```
    fn approx_sin(self) -> Self;
}

macro_rules! impl_approx_sin {
    ($float:ty; $consts:tt) => {
        impl /*const*/ ApproxSin for $float
        {
            fn approx_sin(self) -> Self
            {
                const N: usize = 6;
                const C: [$float; N] = [
                    1.276278962,
                    -0.285261569,
                    0.009118016,
                    -0.000136587,
                    0.000001185,
                    -0.000000007
                ];
                let t: [[$float; N]; N] = ArrayOps::fill(
                    /*const*/ |n| Into::<Option<[$float; N]>>::into(ChebyshevPolynomial::new_of_first_kind(n)).unwrap()
                );
                let p: [$float; N] = t.zip(C)
                    .map2(/*const*/ |(t, c)| t.map2(const |tn| c*tn))
                    .reduce(/*const*/ |a, b| a.zip(b).map2(const |(a, b)| a + b))
                    .unwrap_or_default();

                let mut w = self*$consts::FRAC_2_PI;
                let mut i = 0;
                while i < 4
                {
                    w = (w - 1.0);
                    if i % 2 == 0 && w < 0.0
                    {
                        w = -w;
                    }
                    w %= 4.0;
                    i += 1;
                }
                let w = if w > 1.0 {2.0 - w} else if w < -1.0 {-2.0 - w} else {w};

                let z = 2.0*w*w - 1.0;

                p.evaluate_as_polynomial(z)*w
            }
        }
    };
}
impl_approx_sin!(f32; f32);
impl_approx_sin!(f64; f64);

#[cfg(test)]
mod test
{
    use ::test::Bencher;

    use super::*;
    use crate::tests as t;

    #[test]
    fn sin()
    {
        const RANGE: f32 = f32::TAU;
        t::plot_approx("sin", -RANGE..RANGE, f32::sin, ApproxSin::approx_sin)
    }
    
    #[bench]
    fn sin_benchmark(_: &mut Bencher)
    {
        type F = f64;

        const N: usize = 500;
        const S: usize = 32;

        t::plot_benchmark::<_, _, N, _>(
            "sin",
            [
                &F::sin,
                &ApproxSin::approx_sin
            ],
            -f64::TAU..f64::TAU,
            S
        )
    }
}