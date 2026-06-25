//! Low-level, `unsafe` FFI bindings to the vendored IUPAC InChI 1.07 reference
//! C library.
//!
//! This crate is the unsafe foundation underneath the safe [`inchi`] crate. It
//! exposes the raw C structs and `extern "C"` entry points (`GetINCHI`,
//! `GetINCHIKeyFromINCHI`, `MakeINCHIFromMolfileText`, `GetStructFromINCHI`,
//! the matching `Free*` deallocators, ...) exactly as declared in the upstream
//! `inchi_api.h`. The native C source is vendored and statically linked, so no
//! system InChI installation is required.
//!
//! Almost all users should depend on the high-level [`inchi`] crate instead.
//! Everything here is `unsafe` and mirrors C ownership rules: any `out`
//! parameter populated by a `Get*`/`Make*` call **must** be released with the
//! corresponding `Free*` call to avoid leaking memory allocated by the C side.
//!
//! [`inchi`]: https://docs.rs/inchi
//!
//! # Bindings
//!
//! By default the version-controlled, pre-generated bindings (produced by
//! `bindgen`) are used, so no `libclang` is needed to build. Enable the
//! `regenerate-bindings` feature to regenerate them from the vendored headers
//! at build time.
//!
//! # Safety & threading
//!
//! The InChI C library keeps some `static` internal state and is **not**
//! guaranteed to be thread-safe across concurrent calls. Callers of this crate
//! are responsible for synchronization; the [`inchi`] crate provides it.
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(clippy::all)]

#[cfg(feature = "regenerate-bindings")]
mod bindings {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

#[cfg(not(feature = "regenerate-bindings"))]
#[path = "bindings.rs"]
mod bindings;

pub use bindings::*;

#[cfg(test)]
mod smoke {
    //! Proves the native library links and produces a correct InChI for a
    //! trivial structure built directly into an `inchi_Input`.
    use super::*;
    use std::ffi::CStr;
    use std::os::raw::c_char;

    /// A V2000 molfile for methane (single carbon; InChI adds implicit H).
    const METHANE_MOLFILE: &str = "\n  inchi-sys\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    #[test]
    fn methane_via_molfile_text() {
        // Build NUL-terminated inputs for the C API.
        let moltext = std::ffi::CString::new(METHANE_MOLFILE).expect("no interior NUL");
        let options = std::ffi::CString::new("").expect("no interior NUL");

        let mut out: inchi_Output = unsafe { std::mem::zeroed() };
        let rc = unsafe {
            MakeINCHIFromMolfileText(moltext.as_ptr(), options.as_ptr() as *mut c_char, &mut out)
        };

        // 0 == okay, 1 == warning; both yield a valid InChI string.
        assert!(
            rc == inchi_Ret_OKAY as i32 || rc == inchi_Ret_WARNING as i32,
            "MakeINCHIFromMolfileText returned {rc}"
        );
        assert!(!out.szInChI.is_null(), "no InChI produced");

        let inchi = unsafe { CStr::from_ptr(out.szInChI) }
            .to_str()
            .expect("InChI is valid UTF-8")
            .to_owned();

        // Release everything the C side allocated before asserting.
        unsafe { FreeINCHI(&mut out) };

        assert_eq!(inchi, "InChI=1S/CH4/h1H4");
    }
}
