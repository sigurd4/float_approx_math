#![cfg_attr(not(any(test/*, feature = "std"*/)), no_std)]

#![feature(test)]
#![feature(const_trait_impl)]
#![feature(const_float_bits_conv)]
#![feature(const_fn_floating_point_arithmetic)]
#![feature(const_option)]
#![feature(const_mut_refs)]
#![feature(const_option_ext)]
#![feature(array_zip)]
#![feature(generic_arg_infer)]

#![feature(const_closures)]
#![feature(generic_const_exprs)]

moddef::moddef!(
    flat(pub) mod {
        approx_sqrt for cfg(feature = "sqrt"),
        approx_inv_sqrt for cfg(feature = "sqrt"),

        approx_sin for cfg(feature = "sin"),
        approx_cos for cfg(feature = "sin")
    },
    mod {
        plot for cfg(test)
    }
);

mod f32
{
    pub use core::f32::consts::*;

    #[allow(unused)]
    pub(crate) const EXP_BIAS: u32 = 127;
}
mod f64
{
    pub use core::f64::consts::*;

    #[allow(unused)]
    pub(crate) const EXP_BIAS: u64 = 1023;
}

#[cfg(test)]
extern crate test;

#[cfg(test)]
mod tests {
    use std::{time::Duration, ops::RangeBounds};

    use array_trait::{ArrayOps, Array2dOps};
    use linspace::{LinspaceArray, Linspace};
    
    const PLOT_TARGET: &str = "plots";

    use super::*;

    #[allow(unused)]
    pub fn plot_benchmark<T, R, const N: usize, const M: usize>(
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

    pub fn benchmark<T, R>(x: &[T], f: &dyn Fn(T) -> R) -> Duration
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

    #[allow(unused)]
    pub fn plot_approx<R>(
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

        plot::plot_curves(plot_title, plot_path, [x, x], [y, y_approx])
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
}
