// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*include][include:1]]
#![allow(nonstandard_style)]

#[allow(clippy::all)]
mod bindings {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}
// pub use bindings::*;

use bindings::{integer, logical, FG, FG_END, NEW_X, START};

// update signature to remove mut access from const parameters
extern "C" {
    pub fn setulb(
        n: *const integer,
        m: *const integer,
        x: *mut f64,
        l: *const f64,
        u: *const f64,
        nbd: *const integer,
        f: *mut f64,
        g: *mut f64,
        factr: *const f64,
        pgtol: *const f64,
        wa: *mut f64,
        iwa: *mut integer,
        task: *mut integer,
        iprint: *const integer,
        csave: *mut integer,
        lsave: *mut logical,
        isave: *mut integer,
        dsave: *mut f64,
    ) -> ::std::os::raw::c_int;
}
// include:1 ends here

// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*mods][mods:1]]
mod lbfgsb;
// mods:1 ends here

// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*util][util:1]]
/// Compute function value f for the sample problem.
///
/// Evaluate f(x) and g(x) at current `x`.
fn evaluate(x: &[f64], g: &mut [f64]) -> f64 {
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

    f
}

// #define IS_FG(x) ( ((x)>=FG) ?  ( ((x)<=FG_END) ? 1 : 0 ) : 0 )
fn is_fg(task: i64) -> bool {
    let task = task as u32;
    task >= FG && task <= FG_END
}
// util:1 ends here

// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*param][param:1]]
struct LbfgsbParameter {
    /// On entry m is the maximum number of variable metric corrections allowed
    /// in the limited memory matrix.
    m: usize,

    /// The tolerances in the stopping criteria for function value.
    ///
    /// On entry factr >= 0 is specified by the user. The iteration will stop
    /// when
    ///
    ///   (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
    ///
    /// where epsmch is the machine precision, which is automatically generated
    /// by the code.
    factr: f64,

    /// The tolerances in the stopping criteria for gradient.
    ///
    /// On entry pgtol >= 0 is specified by the user. The iteration will stop
    /// when
    ///
    ///   max{|proj g_i | i = 1, ..., n} <= pgtol
    ///
    /// where pg_i is the ith component of the projected gradient.
    pgtol: f64,

    // iprint controls the frequency and type of output generated:
    //
    //    iprint<0    no output is generated;
    //    iprint=0    print only one line at the last iteration;
    //    0<iprint<99 print also f and |proj g| every iprint iterations;
    //    iprint=99   print details of every iteration except n-vectors;
    //    iprint=100  print also the changes of active set and final x;
    //    iprint>100  print details of every iteration including x and g;
    //
    // When iprint > 0, the file iterate.dat will be created to summarize the
    // iteration.
    iprint: i64,
}

impl Default for LbfgsbParameter {
    fn default() -> Self {
        Self {
            m: 5,
            factr: 1E7,
            pgtol: 1E-5,
            iprint: 1,
        }
    }
}
// param:1 ends here

// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*minimize][minimize:1]]
/// # Parameters
///
fn minimize_lbfgsb<E>(
    x: &mut [f64],
    l: &[f64],
    u: &[f64],
    g: &mut [f64],
    nbd: &[i64],
    param: &LbfgsbParameter,
    mut eval_fn: E,
) where
    E: FnMut(&[f64], &mut [f64]) -> f64,
{
    let n = x.len();
    let m = param.m;

    // 1. Allocate internal workspace arrays

    // wa is a double precision working array of length
    //   (2mmax + 5)nmax + 12mmax^2 + 12mmax.
    let mut wa = vec![0.0; 2 * m * n + 5 * n + 11 * m * m + 8 * m];
    // iwa is an integer working array of length 3nmax.
    let mut iwa = vec![0; 3 * n];

    // csave is a working string of characters of length 60.
    // static char csave[60];
    let mut csave = [0; 60];

    // static double dsave[29];
    //  dsave is a double precision working array of dimension 29.
    // On exit with 'task' = NEW_X, the following information is
    //                                                       available:
    //   dsave(1) = current 'theta' in the BFGS matrix;
    //   dsave(2) = f(x) in the previous iteration;
    //   dsave(3) = factr*epsmch;
    //   dsave(4) = 2-norm of the line search direction vector;
    //   dsave(5) = the machine precision epsmch generated by the code;
    //   dsave(7) = the accumulated time spent on searching for
    //                                                   Cauchy points;
    //   dsave(8) = the accumulated time spent on
    //                                           subspace minimization;
    //   dsave(9) = the accumulated time spent on line search;
    //   dsave(11) = the slope of the line search function at
    //                            the current point of line search;
    //   dsave(12) = the maximum relative step length imposed in
    //                                                     line search;
    //   dsave(13) = the infinity norm of the projected gradient;
    //   dsave(14) = the relative step length in the line search;
    //   dsave(15) = the slope of the line search function at
    //                           the starting point of the line search;
    //   dsave(16) = the square of the 2-norm of the line search
    //                                                direction vector.
    let mut dsave = [0f64; 29];

    // static integer isave[44];
    // isave is an integer working array of dimension 44.
    //   On exit with 'task' = NEW_X, the following information is
    //                                                         available:
    //     isave(22) = the total number of intervals explored in the
    //                     search of Cauchy points;
    //     isave(26) = the total number of skipped BFGS updates before
    //                     the current iteration;
    //     isave(30) = the number of current iteration;
    //     isave(31) = the total number of BFGS updates prior the current
    //                     iteration;
    //     isave(33) = the number of intervals explored in the search of
    //                     Cauchy point in the current iteration;
    //     isave(34) = the total number of function and gradient
    //                     evaluations;
    //     isave(36) = the number of function value or gradient
    //                              evaluations in the current iteration;
    //     if isave(37) = 0  then the subspace argmin is within the box;
    //     if isave(37) = 1  then the subspace argmin is beyond the box;
    //     isave(38) = the number of free variables in the current
    //                     iteration;
    //     isave(39) = the number of active constraints in the current
    //                     iteration;
    //     n + 1 - isave(40) = the number of variables leaving the set of
    //                       active constraints in the current iteration;
    //     isave(41) = the number of variables entering the set of active
    //                     constraints in the current iteration.
    let mut isave = [0; 44];

    // static logical lsave[4];
    // lsave is a logical working array of dimension 4. On exit with 'task' =
    // NEW_X, the following information is available:
    //
    //   If lsave(1) = .true. then the initial X has been replaced by its
    //   projection in the feasible set;
    //
    //   If lsave(2) = .true.  then  the problem is constrained;
    //
    //   If lsave(3) = .true. then each variable has upper and lower bounds;
    let mut lsave = [0; 4];

    // We start the iteration by initializing task.
    // *task = (integer)START;
    let mut task: i64 = START.into();

    // 2. the beginning of the loop
    // for the value of the function at x
    let mut f = 0.0;
    loop {
        // the call to the L-BFGS-B code
        unsafe {
            setulb(
                &(n as i64),        //x
                &(m as i64),        //x
                x.as_mut_ptr(),     //x
                l.as_ptr(),         //x
                u.as_ptr(),         //x
                nbd.as_ptr(),       //x
                &mut f,             //x
                g.as_mut_ptr(),     //x
                &param.factr,       //x
                &param.pgtol,       //x
                wa.as_mut_ptr(),    //x
                iwa.as_mut_ptr(),   //x
                &mut task,          //x
                &param.iprint,      //x
                csave.as_mut_ptr(), //x
                lsave.as_mut_ptr(), //x
                isave.as_mut_ptr(), //x
                dsave.as_mut_ptr(), //x
            );

            if is_fg(task) {
                // the minimization routine has returned to request the
                // function f and gradient g values at the current x.
                // Compute function value f for the sample problem.
                f = eval_fn(x, g);
            // go back to the minimization routine.
            } else if task == NEW_X as i64 {
                // the minimization routine has returned with a new iterate, and we have
                // opted to continue the iteration.
            } else {
                // If task is neither FG nor NEW_X we terminate execution.
                break;
            }
        }
    }
}
// minimize:1 ends here

// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*test][test:1]]
#[test]
fn test_lbfgsb() {
    const N: usize = 25;
    // for the value of the gradient at x.
    let mut g = [0f64; N];
    let mut x = [0f64; N];

    // nbd represents the type of bounds imposed on the
    // variables, and must be specified as follows:
    //     nbd(i)=0 if x(i) is unbounded,
    //            1 if x(i) has only a lower bound,
    //            2 if x(i) has both lower and upper bounds, and
    //            3 if x(i) has only an upper bound.
    let mut nbd = [0; N];

    // We now provide nbd which defines the bounds on the variables:
    // - l   specifies the lower bounds,
    // - u   specifies the upper bounds.
    //
    // First set bounds on the odd-numbered variables.
    //
    // nbd is an integer array of dimension n.
    //   On entry nbd represents the type of bounds imposed on the
    //     variables, and must be specified as follows:
    //     nbd(i)=0 if x(i) is unbounded,
    //            1 if x(i) has only a lower bound,
    //            2 if x(i) has both lower and upper bounds, and
    //            3 if x(i) has only an upper bound.
    //   On exit nbd is unchanged.
    let mut l = [0f64; N];
    let mut u = [0f64; N];
    for i in (1..=N).step_by(2) {
        nbd[i - 1] = 2;
        l[i - 1] = 1.;
        u[i - 1] = 100.;
    }
    // Next set bounds on the even-numbered variables.
    for i in (2..=N).step_by(2) {
        nbd[i - 1] = 2;
        l[i - 1] = -100.;
        u[i - 1] = 100.;
    }

    // We now define the starting point.
    for i in (1..=N) {
        x[i - 1] = 3.;
    }
    println!("     Solving sample problem (Rosenbrock test fcn).");
    println!("      (f = 0.0 at the optimal solution.)");

    // ------- the beginning of the loop ----------
    let param = LbfgsbParameter::default();
    minimize_lbfgsb(&mut x, &l, &u, &mut g, &nbd, &param, evaluate);
}
// test:1 ends here
