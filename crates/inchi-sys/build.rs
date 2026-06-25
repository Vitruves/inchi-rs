//! Build script for `inchi-sys`.
//!
//! Compiles the vendored IUPAC InChI 1.07 reference C library into a static
//! archive and (optionally) regenerates the Rust FFI bindings with `bindgen`.
//!
//! The set of translation units and the preprocessor defines mirror the
//! upstream `INCHI_API/libinchi/gcc/makefile` so the result is byte-for-byte
//! the same library the InChI project ships, just statically linked.

use std::path::{Path, PathBuf};

/// Object files compiled from `INCHI_BASE/src` (the core algorithm + the
/// classic API). Taken verbatim from `INCHI_LIB_OBJS` in the upstream makefile.
const BASE_SRC: &[&str] = &[
    "ichican2",
    "ichicano",
    "ichi_io",
    "ichierr",
    "ichicans",
    "ichiisot",
    "ichimak2",
    "ichimake",
    "ichimap1",
    "ichimap2",
    "ichimap4",
    "ichinorm",
    "ichiparm",
    "ichiprt1",
    "ichiprt2",
    "ichiprt3",
    "ichiqueu",
    "ichiring",
    "ichisort",
    "ichister",
    "ichitaut",
    "ichi_bns",
    "ichiread",
    "ichirvr1",
    "ichirvr2",
    "ichirvr3",
    "ichirvr4",
    "ichirvr5",
    "ichirvr6",
    "ichirvr7",
    "ikey_dll",
    "ikey_base26",
    "mol_fmt1",
    "mol_fmt2",
    "mol_fmt3",
    "mol_fmt4",
    "mol2atom",
    "readinch",
    "runichi",
    "runichi2",
    "runichi3",
    "runichi4",
    "sha2",
    "strutil",
    "util",
    "bcf_s",
];

/// Object files compiled from `INCHI_API/libinchi/src` (the DLL/API layer:
/// `GetINCHI`, `FreeINCHI`, `MakeINCHIFromMolfileText`, ...).
const LIBINCHI_SRC: &[&str] = &[
    "ichilnct",
    "inchi_dll",
    "inchi_dll_main",
    "inchi_dll_a",
    "inchi_dll_a2",
    "inchi_dll_b",
];

/// Object files compiled from `INCHI_API/libinchi/src/ixa` (the InChI
/// eXtensible API).
const IXA_SRC: &[&str] = &[
    "ixa_inchikey_builder",
    "ixa_read_mol",
    "ixa_status",
    "ixa_builder",
    "ixa_mol",
    "ixa_read_inchi",
];

fn main() {
    let vendor = PathBuf::from(env("CARGO_MANIFEST_DIR")).join("vendor/inchi");
    let base = vendor.join("INCHI_BASE/src");
    let libinchi = vendor.join("INCHI_API/libinchi/src");
    let ixa = libinchi.join("ixa");

    // Re-run only when inputs that affect the build change.
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed={}", base.display());

    compile_native(&base, &libinchi, &ixa);

    #[cfg(feature = "regenerate-bindings")]
    regenerate_bindings(&base);
}

fn compile_native(base: &Path, libinchi: &Path, ixa: &Path) {
    let mut build = cc::Build::new();
    build.include(base);

    for name in BASE_SRC {
        build.file(base.join(format!("{name}.c")));
    }
    for name in LIBINCHI_SRC {
        build.file(libinchi.join(format!("{name}.c")));
    }
    for name in IXA_SRC {
        build.file(ixa.join(format!("{name}.c")));
    }

    let target = env("TARGET");
    let msvc = target.contains("msvc");

    // Preprocessor defines, matching the upstream makefile build of the API.
    build.define("TARGET_API_LIB", None);
    if msvc {
        // Windows/MSVC support is structured but unverified; matches the
        // `C_SO_OPTIONS` Windows branch of the upstream makefile.
        build.define("BUILD_LINK_AS_DLL", None);
    } else {
        build.define("COMPILE_ANSI_ONLY", None);
    }

    // The InChI sources rely on type punning that strict aliasing would break.
    build.flag_if_supported("-fno-strict-aliasing");

    // The reference C is warning-heavy by design; do not let that fail or spam
    // the build. `cc` already omits `-Werror` unless asked.
    build.warnings(false);
    if !msvc {
        for flag in [
            "-Wno-unused-but-set-variable",
            "-Wno-unused-variable",
            "-Wno-unused-function",
            "-Wno-implicit-fallthrough",
            "-Wno-format",
            "-Wno-deprecated-declarations",
        ] {
            build.flag_if_supported(flag);
        }
    }

    build.compile("inchi");
    // `cc` emits the `cargo:rustc-link-lib`/`-search` directives for us.
}

#[cfg(feature = "regenerate-bindings")]
fn regenerate_bindings(base: &Path) {
    let header = base.join("inchi_api.h");
    println!("cargo:rerun-if-changed={}", header.display());

    let bindings = bindgen::Builder::default()
        .header(header.to_string_lossy())
        .clang_arg(format!("-I{}", base.display()))
        .clang_arg("-DTARGET_API_LIB")
        .clang_arg("-DCOMPILE_ANSI_ONLY")
        // Public classic-API surface.
        .allowlist_function("Get(Std)?INCHI(Ex)?")
        .allowlist_function("GetStructFrom(Std)?INCHI(Ex)?")
        .allowlist_function("Free(Std)?INCHI")
        .allowlist_function("FreeStructFrom(Std)?INCHI(Ex)?")
        .allowlist_function("Get(Std)?INCHIKeyFrom(Std)?INCHI")
        .allowlist_function("GetINCHIfromINCHI")
        .allowlist_function("MakeINCHIFromMolfileText")
        .allowlist_function("CheckINCHI(Key)?")
        .allowlist_function("GetStringLength")
        .allowlist_function("Free_(std_)?inchi_Input")
        .allowlist_function("Get_(std_)?inchi_Input_FromAuxInfo")
        .allowlist_type("inchi_.*")
        .allowlist_type("RetVal.*")
        .allowlist_type("INCHI_.*")
        .allowlist_var("INCHI_.*")
        .allowlist_var(
            "(MAXVAL|NO_ATOM|ATOM_EL_LEN|NUM_H_ISOTOPES|ISOTOPIC_SHIFT_FLAG|ISOTOPIC_SHIFT_MAX)",
        )
        // Plain integer constants: sound for values the C library returns
        // (an out-of-range return code can never be UB), and the constants are
        // accessible bare (e.g. `inchi_Ret_OKAY`).
        .default_enum_style(bindgen::EnumVariation::Consts)
        // The C enum names are already namespaced (`inchi_Ret_OKAY`,
        // `INCHI_BOND_TYPE_SINGLE`); don't double-prefix with the tag name.
        .prepend_enum_name(false)
        .derive_default(true)
        .derive_debug(true)
        .layout_tests(true)
        .generate_comments(false)
        .generate()
        .expect("failed to generate InChI bindings");

    let out = PathBuf::from(env("OUT_DIR")).join("bindings.rs");
    bindings
        .write_to_file(&out)
        .expect("failed to write generated bindings");
}

fn env(key: &str) -> String {
    std::env::var(key).unwrap_or_else(|_| panic!("environment variable `{key}` not set"))
}
