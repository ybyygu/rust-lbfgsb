// [[file:../lbfgsb.note::86501f03][86501f03]]
#![allow(nonstandard_style)]

#[allow(clippy::all)]
mod bindings {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}
// pub use bindings::*;

use bindings::{integer, logical};

// update signature to remove mut access from const parameters
extern "C" {
    #[allow(clashing_extern_declarations)]
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
// 86501f03 ends here

// [[file:../lbfgsb.note::7c3692fa][7c3692fa]]
mod lbfgsb;
// 7c3692fa ends here

// [[file:../lbfgsb.note::*exports][exports:1]]
pub use crate::lbfgsb::lbfgsb;
pub use crate::lbfgsb::LbfgsbParameter;
pub use crate::lbfgsb::LbfgsbProblem;
pub use crate::lbfgsb::LbfgsbState;
// exports:1 ends here
