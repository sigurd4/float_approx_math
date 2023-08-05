use super::*;

use polynomial_ops::*;
use array_trait::*;

#[const_trait]
pub trait ApproxSin
{
    /// https://namoseley.wordpress.com/2012/09/26/arduinoavr-and-zx-spectrum-sin-routines/
    /// 
    /// 
    fn approx_sin(self) -> Self;
}

macro_rules! impl_approx_sin {
    ($float:ty; $consts:tt) => {
        impl const ApproxSin for $float
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
                const T: [[$float; N]; N] = ArrayOps::fill(
                    const |n| Into::<Option<[$float; N]>>::into(ChebyshevPolynomial::new_of_first_kind(n)).unwrap()
                );
                const P: [$float; N] = T.zip2(C)
                    .map2(const |(t, c)| t.map2(const |tn| c*tn))
                    .reduce(const |a, b| a.zip2(b).map2(const |(a, b)| a + b))
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

                P.evaluate_as_polynomial(z)*w
            }
        }
    };
}
impl_approx_sin!(f32; f32);
impl_approx_sin!(f64; f64);