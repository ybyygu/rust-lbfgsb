// [[file:../lbfgsb.note::*imports][imports:1]]
use crate::*;
use bindings::{FG, FG_END, NEW_X, START};

use anyhow::Result;
// imports:1 ends here

// [[file:../lbfgsb.note::*util][util:1]]
// #define IS_FG(x) ( ((x)>=FG) ?  ( ((x)<=FG_END) ? 1 : 0 ) : 0 )
fn is_fg(task: i64) -> bool {
    let task = task as u32;
    task >= FG && task <= FG_END
}
// util:1 ends here

// [[file:../lbfgsb.note::*param][param:1]]
/// L-BFGS-B algorithm parameters
pub struct LbfgsbParameter {
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
            iprint: -1,
        }
    }
}
// param:1 ends here

// [[file:../lbfgsb.note::*problem][problem:1]]
pub struct LbfgsbProblem<E>
where
    E: FnMut(&[f64], &mut [f64]) -> Result<f64>,
{
    x: Vec<f64>,
    g: Vec<f64>,
    f: f64,
    l: Vec<f64>,
    u: Vec<f64>,
    nbd: Vec<i64>,
    eval_fn: E,
}

impl<E> LbfgsbProblem<E>
where
    E: FnMut(&[f64], &mut [f64]) -> Result<f64>,
{
    pub fn build(x: Vec<f64>, eval_fn: E) -> Self {
        let n = x.len();
        Self {
            x,
            g: vec![0.0; n],
            f: 0.0,
            l: vec![0.0; n],
            u: vec![0.0; n],
            nbd: vec![0; n],
            eval_fn,
        }
    }

    /// Set lower bounds and upper bounds for input variables
    pub fn set_bounds<B>(&mut self, bounds: B)
    where
        B: IntoIterator<Item = (Option<f64>, Option<f64>)>,
    {
        // nbd represents the type of bounds imposed on the variables, and must be
        // specified as follows:
        //
        //   nbd(i)=0 if x(i) is unbounded,
        //          1 if x(i) has only a lower bound,
        //          2 if x(i) has both lower and upper bounds, and
        //          3 if x(i) has only an upper bound.
        for (i, b) in bounds.into_iter().enumerate() {
            match b {
                // both lower and upper bonds
                (Some(l), Some(u)) => {
                    self.l[i] = l;
                    self.u[i] = u;
                    self.nbd[i] = 2;
                }
                // unbounded
                (None, None) => {
                    self.nbd[i] = 0;
                }
                // has only a lower bound
                (Some(l), None) => {
                    self.l[i] = l;
                    self.nbd[i] = 1;
                }
                // has only a upper bound
                (None, Some(u)) => {
                    self.u[i] = u;
                    self.nbd[i] = 3;
                }
            }
        }
    }
}
// problem:1 ends here

// [[file:../lbfgsb.note::*state][state:1]]
pub struct LbfgsbState<E>
where
    E: FnMut(&[f64], &mut [f64]) -> Result<f64>,
{
    problem: LbfgsbProblem<E>,

    param: LbfgsbParameter,

    /// wa is a double precision working array of length:
    ///   (2mmax + 5)nmax + 12mmax^2 + 12mmax.
    wa: Vec<f64>,

    // iwa is an integer working array of length 3nmax.
    iwa: Vec<i64>,

    // csave is a working string of characters of length 60.
    // static char csave[60];
    csave: [i64; 60],
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
    dsave: [f64; 29],

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
    isave: [i64; 44],
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
    lsave: [i64; 4],

    // Note in original fortran version:
    //
    // task is a working string of characters of length 60 indicating
    // the current job when entering and leaving this subroutine.
    //
    // Note in L-BFGS-B-C
    //
    // Modified L-BFGS-B to use integers instead of strings, for testing the
    // "task"
    task: i64,
}

impl<E> LbfgsbState<E>
where
    E: FnMut(&[f64], &mut [f64]) -> Result<f64>,
{
    pub(crate) fn new(problem: LbfgsbProblem<E>, param: LbfgsbParameter) -> Self {
        let n = problem.x.len();
        let m = param.m;
        // wa is a double precision working array of length
        //   (2mmax + 5)nmax + 12mmax^2 + 12mmax.
        let wa = vec![0.0; 2 * m * n + 5 * n + 11 * m * m + 8 * m];

        // iwa is an integer working array of length 3nmax.
        let iwa = vec![0; 3 * n];

        Self {
            csave: [0; 60],
            dsave: [0.0; 29],
            isave: [0; 44],
            lsave: [0; 4],
            task: START.into(),
            problem,
            param,
            wa,
            iwa,
        }
    }

    pub(crate) fn minimize(&mut self) -> Result<()> {
        let f = &mut self.problem.f;
        let x = &mut self.problem.x;
        let g = &mut self.problem.g;
        let l = &self.problem.l;
        let u = &self.problem.u;
        let nbd = &self.problem.nbd;

        let param = &self.param;
        let n = x.len();
        let m = param.m;
        loop {
            unsafe {
                crate::setulb(
                    &(n as i64),             //x
                    &(m as i64),             //x
                    x.as_mut_ptr(),          //x
                    l.as_ptr(),              //x
                    u.as_ptr(),              //x
                    nbd.as_ptr(),            //x
                    f,                       //x
                    g.as_mut_ptr(),          //x
                    &param.factr,            //x
                    &param.pgtol,            //x
                    self.wa.as_mut_ptr(),    //x
                    self.iwa.as_mut_ptr(),   //x
                    &mut self.task,          //x
                    &param.iprint,           //x
                    self.csave.as_mut_ptr(), //x
                    self.lsave.as_mut_ptr(), //x
                    self.isave.as_mut_ptr(), //x
                    self.dsave.as_mut_ptr(), //x
                );
            }
            if is_fg(self.task) {
                // the minimization routine has returned to request the
                // function f and gradient g values at the current x.
                // Compute function value f for the sample problem.
                *f = (self.problem.eval_fn)(x, g)?;
            // go back to the minimization routine.
            } else if self.task == NEW_X as i64 {
                // the minimization routine has returned with a new iterate, and we have
                // opted to continue the iteration.
            } else {
                // If task is neither FG nor NEW_X we terminate execution.
                break;
            }
        }

        Ok(())
    }

    /// Final function value f(x)
    pub fn fx(&self) -> f64 {
        self.problem.f
    }

    /// Final evaluated gradient g(x)
    pub fn gx(&self) -> &[f64] {
        &self.problem.g
    }

    /// Final optimized `x`
    pub fn x(&self) -> &[f64] {
        &self.problem.x
    }
}
// state:1 ends here

// [[file:../lbfgsb.note::*pub][pub:1]]
/// Minimize a scalar function of one or more variables using the L-BFGS-B
/// algorithm.
///
/// # Parameters
///
/// - bounds: a slice of tuple setting lower and upper bounds.
/// - eval_fn: a closure evaluating f(x) and g(x). Returning Err value will cancel minimization.
///
/// # Return
///
/// - Returns final state containing x, f(x), g(x).
pub fn lbfgsb<E>(x: Vec<f64>, bounds: &[(f64, f64)], eval_fn: E) -> Result<LbfgsbState<E>>
where
    E: FnMut(&[f64], &mut [f64]) -> Result<f64>,
{
    assert_eq!(x.len(), bounds.len());

    let param = LbfgsbParameter::default();
    let mut problem = LbfgsbProblem::build(x, eval_fn);
    let bounds = bounds.into_iter().copied().map(|(l, u)| (Some(l), Some(u)));
    problem.set_bounds(bounds);

    let mut state = LbfgsbState::new(problem, param);
    state.minimize()?;

    Ok(state)
}
// pub:1 ends here
