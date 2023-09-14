#![cfg_attr(not(any(test/*, feature = "std"*/)), no_std)]

#![feature(test)]
#![feature(const_trait_impl)]
#![feature(const_float_bits_conv)]
#![feature(const_fn_floating_point_arithmetic)]
#![feature(const_option)]
#![feature(const_closures)]
#![feature(const_mut_refs)]
#![feature(const_option_ext)]
#![feature(array_zip)]
#![feature(generic_const_exprs)]
#![feature(generic_arg_infer)]

moddef::moddef!(
    flat(pub) mod {
        approx_sqrt,
        approx_inv_sqrt,

        approx_sin,
        approx_cos
    },
    mod {
        plot for cfg(test)
    }
);

mod f32
{
    pub use core::f32::consts::*;

    pub(crate) const EXP_BIAS: u32 = 127;
}
mod f64
{
    pub use core::f64::consts::*;

    pub(crate) const EXP_BIAS: u64 = 1023;
}

#[cfg(test)]
extern crate test;

#[cfg(test)]
mod tests {
    use std::{time::Duration, ops::RangeBounds};

    use array_trait::{ArrayOps, Array2dOps};
    use linspace::{LinspaceArray, Linspace};
    use test::Bencher;
    
    const PLOT_TARGET: &str = "plots";

    use super::*;
    
    #[test]
    fn inv_sqrt()
    {
        const RANGE: f32 = 10.0;
        
        plot_approx("inv_sqrt_0", 0.1..RANGE, |x| x.sqrt().recip(), ApproxInvSqrt::approx_inv_sqrt::<0>);
        plot_approx("inv_sqrt_1", 0.1..RANGE, |x| x.sqrt().recip(), ApproxInvSqrt::approx_inv_sqrt::<1>);
        plot_approx("inv_sqrt_2", 0.1..RANGE, |x| x.sqrt().recip(), ApproxInvSqrt::approx_inv_sqrt::<2>);
        plot_approx("inv_sqrt_3", 0.1..RANGE, |x| x.sqrt().recip(), ApproxInvSqrt::approx_inv_sqrt::<3>);
    }

    #[test]
    fn sin()
    {
        const RANGE: f32 = f32::TAU;
        plot_approx("sin", -RANGE..RANGE, f32::sin, ApproxSin::approx_sin)
    }
    
    #[bench]
    fn sin_benchmark(_: &mut Bencher)
    {
        type F = f64;

        const N: usize = 500;
        const S: usize = 32;

        plot_benchmark::<_, _, N, _>(
            "sin",
            [
                &F::sin,
                &ApproxSin::approx_sin
            ],
            -f64::TAU..f64::TAU,
            S
        )
    }
    
    #[test]
    fn cos()
    {
        const RANGE: f32 = f32::TAU;
        plot_approx("cos", -RANGE..RANGE, f32::cos, ApproxCos::approx_cos)
    }
    
    #[bench]
    fn cos_benchmark(_: &mut Bencher)
    {
        type F = f64;

        const N: usize = 500;
        const S: usize = 32;

        plot_benchmark::<_, _, N, _>(
            "cos",
            [
                &F::cos,
                &ApproxCos::approx_cos
            ],
            -f64::TAU..f64::TAU,
            S
        )
    }

    #[bench]
    fn sqrt_benchmark(_: &mut Bencher)
    {
        type F = f64;

        const N: usize = 1000;
        const S: usize = 32;

        plot_benchmark::<_, _, N, _>(
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
    
            println!("sqrt dt = {:?}", benchmark(&x, &F::sqrt));
            
            println!("approx_sqrt::<0> dt = {:?}", benchmark(&x, &|x| x.approx_sqrt::<0>()));
            println!("approx_sqrt::<1> dt = {:?}", benchmark(&x, &|x| x.approx_sqrt::<1>()));
            println!("approx_sqrt::<2> dt = {:?}", benchmark(&x, &|x| x.approx_sqrt::<2>()));
            println!("approx_sqrt::<3> dt = {:?}", benchmark(&x, &|x| x.approx_sqrt::<3>()));
            println!("approx_sqrt::<4> dt = {:?}", benchmark(&x, &|x| x.approx_sqrt::<4>()));
        }
    }

    #[bench]
    fn inv_sqrt_benchmark(_: &mut Bencher)
    {
        type F = f64;

        const N: usize = 1000;
        const S: usize = 32;

        plot_benchmark::<_, _, N, _>(
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
    
            println!("sqrt dt = {:?}", benchmark(&x, &|x: F| x.sqrt().recip()));
            
            println!("approx_inv_sqrt::<0> dt = {:?}", benchmark(&x, &|x| x.approx_inv_sqrt::<0>()));
            println!("approx_inv_sqrt::<1> dt = {:?}", benchmark(&x, &|x| x.approx_inv_sqrt::<1>()));
            println!("approx_inv_sqrt::<2> dt = {:?}", benchmark(&x, &|x| x.approx_inv_sqrt::<2>()));
            println!("approx_inv_sqrt::<3> dt = {:?}", benchmark(&x, &|x| x.approx_inv_sqrt::<3>()));
            println!("approx_inv_sqrt::<4> dt = {:?}", benchmark(&x, &|x| x.approx_inv_sqrt::<4>()));
        }
    }

    fn plot_benchmark<T, R, const N: usize, const M: usize>(
        fn_name: &str,
        func: [&dyn Fn(T) -> R; M],
        x: impl RangeBounds<T> + Linspace<T> + Clone,
        smoothing: usize
    ) -> ()
    where
        T: Clone
    {
        let n: [usize; N] = (1..=N).linspace_array();
        let t = n
            .map(move |n| {
                let x = x.clone()
                    .linspace(n);
                func.map(|f| (0..smoothing).map(|_| benchmark(&x, f).as_secs_f32())
                    .reduce(|a, b| a + b)
                    .unwrap()/smoothing as f32
                )
            }).transpose();
        
        let plot_title: &str = &format!("{fn_name}(x) benchmark");
        let plot_path: &str = &format!("{PLOT_TARGET}/{fn_name}_benchmark.png");
        
        plot::plot_curves(plot_title, plot_path, [n.map(|n| n as f32); M], t).expect("Plot error")
    }

    fn benchmark<T, R>(x: &[T], f: &dyn Fn(T) -> R) -> Duration
    where
        T: Clone
    {
        use std::time::SystemTime;

        let x = x.to_vec();
        let t0 = SystemTime::now();
        x.into_iter().for_each(|x| {f(x);});
        t0.elapsed().unwrap()
    }

    const N: usize = 2048;

    fn plot_approx<R>(
        fn_name: &str,
        range: R,
        func: impl Fn(f32) -> f32,
        approx: impl Fn(f32) -> f32
    )
    where
        R: RangeBounds<f32> + LinspaceArray<f32, N>
    {
        let x: [f32; N] = range.linspace_array();
        let y_approx = x.map(approx);

        let y = x.map(func);

        let plot_title: &str = &format!("{fn_name}(x)");
        let plot_path: &str = &format!("{PLOT_TARGET}/{fn_name}.png");

        plot::plot_curves(plot_title, plot_path, [x, x], [y_approx, y])
            .expect("Plot error");

        let (avg_error, max_abs_error) = y.zip(y_approx)
            .map(|(y, y_approx)| y - y_approx)
            .map(|y| (y, y.abs()))
            .reduce(|a, b| (a.0 + b.0, a.1.max(b.1)))
            .map(|(sum_error, max_abs_error)| (sum_error/N as f32, max_abs_error))
            .unwrap_or_default();
        println!("Average Error: {}", avg_error);
        println!("Max |Error|: {}", max_abs_error);
    }

    #[test]
    fn sqrt_error()
    {
        const X: f64 = 2.0;
        const Y: f64 = X.approx_sqrt::<0>();

        println!("{}", X.sqrt());
        println!("{}", Y);
        println!("error = {}", (Y - X.sqrt())/X.sqrt());
    }
}
