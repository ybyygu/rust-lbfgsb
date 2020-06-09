// [[file:~/Workspace/Programming/rust-libs/l-bfgs-b-c/lbfgsb.note::*build.rs][build.rs:1]]
use bindgen;
use cc;

use std::collections::HashSet;
use std::env;
use std::path::PathBuf;

#[derive(Debug)]
struct IgnoreMacros(HashSet<String>);

impl bindgen::callbacks::ParseCallbacks for IgnoreMacros {
    fn will_parse_macro(&self, name: &str) -> bindgen::callbacks::MacroParsingBehavior {
        if self.0.contains(name) {
            bindgen::callbacks::MacroParsingBehavior::Ignore
        } else {
            bindgen::callbacks::MacroParsingBehavior::Default
        }
    }
}

fn main() {
    println!("cargo:rustc-link-lib=m");
    cc::Build::new()
        .cpp(false)
        .include("lib/src")
        .file("lib/src/lbfgsb.c")
        .file("lib/src/linesearch.c")
        .file("lib/src/subalgorithms.c")
        .file("lib/src/print.c")
        .file("lib/src/linpack.c")
        .file("lib/src/miniCBLAS.c")
        .file("lib/src/timer.c")
        .compile("liblbfgsb.a");

    let ignored_macros = IgnoreMacros(
        vec![
            "FP_INFINITE".into(),
            "FP_NAN".into(),
            "FP_NORMAL".into(),
            "FP_SUBNORMAL".into(),
            "FP_ZERO".into(),
            "IPPORT_RESERVED".into(),
        ]
        .into_iter()
        .collect(),
    );

    let bindings = bindgen::Builder::default()
        .header("lib/src/lbfgsb.h")
        .parse_callbacks(Box::new(ignored_macros))
        .rustfmt_bindings(true)
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
// build.rs:1 ends here
