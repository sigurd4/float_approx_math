use super::*;

use polynomial_ops::*;
use array__ops::*;

#[const_trait]
pub trait ApproxCos
{
    /// Calculates an approximation of a cosine, using the Z80 / ZX Spectrum algorithm.
    /// 
    /// https://namoseley.wordpress.com/2012/09/26/arduinoavr-and-zx-spectrum-sin-routines/
    /// 
    /// This approximation is significantly slower than the built-in cos method in the standard library,
    /// however this method can be run at compile-time.
    /// 
    /// # Example
    /// 
    /// ```rust
    /// #![feature(const_trait_impl)]
    /// 
    /// use float_approx_math::ApproxCos;
    /// 
    /// const X: f32 = 2.0;
    /// let y: f32 = X.approx_cos();
    ///
    /// assert!((y - X.cos()).abs() < 0.0000005); // Less than 0.0000005 abs error
    /// ```
    fn approx_cos(self) -> Self;
}

macro_rules! impl_approx_cos {
    ($float:ty; $consts:tt) => {
        impl /*const*/ ApproxCos for $float
        {
            fn approx_cos(self) -> Self
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

                let mut w = self*$consts::FRAC_2_PI + 1.0;
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
impl_approx_cos!(f32; f32);
impl_approx_cos!(f64; f64);

#[cfg(test)]
mod test
{
    use ::test::Bencher;

    use super::*;
    use crate::tests as t;

    #[test]
    fn cos()
    {
        const RANGE: f32 = f32::TAU;
        t::plot_approx("cos", -RANGE..RANGE, f32::cos, ApproxCos::approx_cos)
    }
    
    #[bench]
    fn cos_benchmark(_: &mut Bencher)
    {
        type F = f64;

        const N: usize = 500;
        const S: usize = 32;

        t::plot_benchmark::<_, _, N, _>(
            "cos",
            [
                &F::cos,
                &ApproxCos::approx_cos
            ],
            -f64::TAU..f64::TAU,
            S
        )
    }
}