// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*imports][imports:1]]
use crate::*;
// imports:1 ends here

// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*problem][problem:1]]
pub struct LbfgsbProblem<E>
where
    E: FnMut(&[f64], &mut [f64]) -> f64,
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
    E: FnMut(&[f64], &mut [f64]) -> f64,
{
    fn build(x: Vec<f64>, eval_fn: E) -> Self {
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
}
// problem:1 ends here

// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*state][state:1]]
pub struct LbfgsbState<E>
where
    E: FnMut(&[f64], &mut [f64]) -> f64,
{
    problem: LbfgsbProblem<E>,

    param: LbfgsbParameter,

    /// wa is a double precision working array of length:
    ///   (2mmax + 5)nmax + 12mmax^2 + 12mmax.
    wa: Vec<f64>,

    // iwa is an integer working array of length 3nmax.
    iwa: Vec<i64>,

    csave: [i64; 60],
    dsave: [f64; 29],
    isave: [i64; 44],
    lsave: [i64; 4],

    task: i64,
}

impl<E> LbfgsbState<E>
where
    E: FnMut(&[f64], &mut [f64]) -> f64,
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
            task: crate::START.into(),
            problem,
            param,
            wa,
            iwa,
        }
    }

    pub fn propagate(&mut self) {
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
                *f = (self.problem.eval_fn)(x, g);
                break;
            // go back to the minimization routine.
            } else if self.task == NEW_X as i64 {
                // the minimization routine has returned with a new iterate, and we have
                // opted to continue the iteration.
            } else {
                // If task is neither FG nor NEW_X we terminate execution.
                break;
            }
        }
    }
}
// state:1 ends here
