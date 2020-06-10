// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*driver1.rs][driver1.rs:1]]
use anyhow::Result;
use lbfgsb::*;
use vecfx::*; // for calculate gradient norm

/// Compute function value f for the sample problem.
///
/// Evaluate f(x) and g(x) at current `x`.
///
/// Returns Err will cancel L-BFGS-B minimization.
fn evaluate(x: &[f64], g: &mut [f64]) -> Result<f64> {
    let n = x.len();
    let mut d1 = x[0] - 1.;
    let mut f = d1 * d1 * 0.25;
    for i in (2..=n) {
        /* Computing 2nd power */
        let d2 = x[i - 2];
        /* Computing 2nd power */
        d1 = x[i - 1] - d2 * d2;
        f += d1 * d1;
    }
    f *= 4.;

    // Compute gradient g for the sample problem.
    /* Computing 2nd power */
    let mut d1 = x[0];
    let mut t1 = x[1] - d1 * d1;
    g[0] = (x[0] - 1.) * 2. - x[0] * 16. * t1;

    for i in (2..=n - 1) {
        let t2 = t1;
        /* Computing 2nd power */
        d1 = x[i - 1];
        t1 = x[i] - d1 * d1;
        g[i - 1] = t2 * 8. - x[i - 1] * 16. * t1;
    }
    g[n - 1] = t1 * 8.;

    Ok(f)
}

#[test]
fn test_lbfgs() -> Result<()> {
    const N: usize = 25;
    let mut x = vec![0f64; N];

    // We now provide nbd which defines the bounds on the variables:
    // - l   specifies the lower bounds,
    // - u   specifies the upper bounds.
    //
    // First set bounds on the odd-numbered variables.
    let mut l = vec![0f64; N];
    let mut u = vec![0f64; N];
    for i in (1..=N).step_by(2) {
        l[i - 1] = 1.;
        u[i - 1] = 100.;
    }
    // Next set bounds on the even-numbered variables.
    for i in (2..=N).step_by(2) {
        l[i - 1] = -100.;
        u[i - 1] = 100.;
    }

    // We now define the starting point.
    for i in (1..=N) {
        x[i - 1] = 3.;
    }
    println!("     Solving sample problem (Rosenbrock test fcn).");
    println!("      (f = 0.0 at the optimal solution.)");
    let bounds: Vec<_> = l.into_iter().zip(u.into_iter()).collect();
    let opt = lbfgsb(x, &bounds, evaluate)?;
    assert!(opt.fx() <= 1e-8);
    assert!(dbg!(opt.gx().vec2norm()) < 1e-3);

    Ok(())
}
// driver1.rs:1 ends here
