//! Safe, idiomatic Rust bindings to the IUPAC InChI reference library for
//! generating **InChI** and **InChIKey** chemical identifiers.
//!
//! The InChI C library is vendored and statically linked through the
//! [`inchi-sys`](inchi_sys) crate, so there is nothing to install and no
//! network access at build time. This crate wraps it in a panic-free,
//! [`Result`]-based API where every native allocation is released
//! automatically (RAII) and no `unsafe` is exposed at the public boundary.
//!
//! # Quickstart
//!
//! Generate an InChI and InChIKey from a Molfile. Pass `()` for default
//! options (or [`Options`] to customize):
//!
//! ```
//! use inchi::{from_molfile, inchikey};
//!
//! // A minimal V2000 Molfile for methane (one carbon; hydrogens are implicit).
//! let methane = "\n  example\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n\
//!     \x20   0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
//!
//! let result = from_molfile(methane, ())?;
//! assert_eq!(result.inchi(), "InChI=1S/CH4/h1H4");
//!
//! let key = inchikey(result.inchi())?;
//! assert_eq!(key, "VNWKTOKETHGBQD-UHFFFAOYSA-N");
//! # Ok::<(), inchi::InchiError>(())
//! ```
//!
//! Or build a structure programmatically:
//!
//! ```
//! use inchi::{Molecule, Atom, BondOrder};
//!
//! let mut ethanol = Molecule::new();
//! let c1 = ethanol.add_atom(Atom::new("C"));
//! let c2 = ethanol.add_atom(Atom::new("C"));
//! let o = ethanol.add_atom(Atom::new("O"));
//! ethanol.add_bond(c1, c2, BondOrder::Single)?;
//! ethanol.add_bond(c2, o, BondOrder::Single)?;
//!
//! assert_eq!(ethanol.to_inchi(())?.inchi(), "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3");
//! # Ok::<(), inchi::InchiError>(())
//! ```
//!
//! # API overview
//!
//! This crate wraps the complete InChI generation, parsing, and validation API.
//! For concepts and worked examples (standard vs. non-standard InChI, stereo
//! conventions, polymers), see the [`guide`].
//!
//! **Structure → InChI**
//! - [`from_molfile`] — Molfile/SDF text (`MakeINCHIFromMolfileText`).
//! - [`Molecule::to_inchi`] — a programmatically built structure (`GetINCHI`,
//!   or `GetINCHIEx` when [polymer units](Molecule::add_polymer_unit) are present).
//!
//! **InChI → structure**
//! - [`struct_from_inchi`] — atoms, bonds, and 0D stereo (`GetStructFromINCHI`).
//! - [`struct_from_inchi_ex`] — the above plus polymer data (`GetStructFromINCHIEx`).
//! - [`struct_from_aux_info`] — recover a structure from an `AuxInfo` block.
//!
//! **InChI → InChIKey** (feature `key`)
//! - [`inchikey`] / [`inchikey_with_hashes`] (`GetINCHIKeyFromINCHI`).
//!
//! **Conversion & validation**
//! - [`inchi_to_inchi`] — re-derive/normalize an InChI (`GetINCHIfromINCHI`).
//! - [`check_inchi`] — validate an InChI string (`CheckINCHI`).
//! - [`check_inchikey`] — validate an InChIKey string (`CheckINCHIKey`).
//!
//! Generation is tuned through [`Options`], including the experimental
//! [`Polymers`] extension. Native allocations are released automatically.
//!
//! # Feature flags
//!
//! | Feature               | Default | Description                                                      |
//! | --------------------- | ------- | ---------------------------------------------------------------- |
//! | `key`                 | yes     | Enables [`inchikey`]/[`inchikey_with_hashes`] (`GetINCHIKeyFromINCHI`). |
//! | `regenerate-bindings` | no      | Regenerates the FFI bindings with `bindgen` (needs `libclang`).  |
//!
//! # Thread safety
//!
//! The underlying InChI C library keeps internal `static` state and is not
//! guaranteed to be re-entrant, so every call into it is serialized behind a
//! global lock. All functions in this crate are therefore safe to call from
//! multiple threads, but calls do not run concurrently with one another. The
//! public types are [`Send`] and [`Sync`].
//!
//! # Minimum supported Rust version
//!
//! This crate supports Rust **1.74** and later. Raising the MSRV is treated as
//! a semver-breaking change.

#![deny(missing_docs)]
#![deny(clippy::unwrap_used)]
#![deny(clippy::expect_used)]
#![deny(clippy::panic)]
#![deny(clippy::indexing_slicing)]
#![cfg_attr(docsrs, feature(doc_cfg))]

pub mod guide;

mod convert;
mod error;
mod molecule;
mod options;
mod output;
mod polymer;
mod raw;
mod structure;
mod validate;

pub use convert::inchi_to_inchi;
pub use error::{InchiError, Result, Status};
pub use molecule::{Atom, BondOrder, ImplicitH, Molecule, Parity, Radical, Stereo};
pub use options::{Options, Polymers, StereoMode};
pub use output::InchiOutput;
pub use polymer::{PolymerConnection, PolymerSubtype, PolymerUnit, PolymerUnitKind};
pub use structure::{
    struct_from_aux_info, struct_from_inchi, struct_from_inchi_ex, ExtendedStructure, Structure,
    StructureAtom, StructureBond,
};
pub use validate::{check_inchi, check_inchikey, InchiKeyValidity, InchiValidity};

use std::os::raw::c_char;

/// Generates an InChI directly from a Molfile (V2000/V3000 SDF text).
///
/// This is the most convenient entry point: hand it the text of a `.mol`/`.sdf`
/// record and it returns the InChI, auxiliary information, and any warnings.
///
/// # Errors
///
/// Returns [`InchiError::Failed`] if the library could not produce an InChI,
/// [`InchiError::EmptyResult`] if it reported success but produced nothing, and
/// [`InchiError::InteriorNul`] if `molfile` contains a NUL byte.
///
/// ```
/// use inchi::from_molfile;
///
/// let water = "\n  ex\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n\
///     \x20   0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
/// let out = from_molfile(water, ())?;
/// assert_eq!(out.inchi(), "InChI=1S/H2O/h1H2");
/// # Ok::<(), inchi::InchiError>(())
/// ```
///
/// The second argument accepts `()` for defaults or an [`Options`] for control
/// (stereo, fixed-H, polymers, …).
pub fn from_molfile(molfile: impl AsRef<str>, options: impl Into<Options>) -> Result<InchiOutput> {
    let options = options.into();
    let moltext = raw::to_cstring(molfile.as_ref())?;
    let opts = raw::to_cstring(&options.to_arg_string())?;

    let _guard = raw::lock();
    let mut out = raw::OutputGuard::new();
    // SAFETY: both C-string arguments are valid and NUL-terminated for the
    // duration of the call; `out` is freed by its `Drop`. Serialized by `_guard`.
    let rc = unsafe {
        inchi_sys::MakeINCHIFromMolfileText(
            moltext.as_ptr(),
            opts.as_ptr() as *mut c_char,
            out.as_mut_ptr(),
        )
    };
    drop(moltext);
    drop(opts);
    build_output(rc, &out)
}

/// Computes the 27-character InChIKey for an InChI string.
///
/// # Errors
///
/// Returns [`InchiError::Failed`] if the input is not a valid InChI, and
/// [`InchiError::InteriorNul`] if it contains a NUL byte.
///
/// ```
/// use inchi::inchikey;
///
/// assert_eq!(inchikey("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")?, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N");
/// # Ok::<(), inchi::InchiError>(())
/// ```
#[cfg(feature = "key")]
#[cfg_attr(docsrs, doc(cfg(feature = "key")))]
pub fn inchikey(inchi: impl AsRef<str>) -> Result<String> {
    Ok(inchikey_with_hashes(inchi, false, false)?.key)
}

/// An InChIKey together with the optional 256-bit hash extensions.
///
/// The extension blocks are only populated when requested via
/// [`inchikey_with_hashes`]; otherwise they are [`None`]. Each, when present, is
/// a string of up to 64 hexadecimal characters.
#[cfg(feature = "key")]
#[cfg_attr(docsrs, doc(cfg(feature = "key")))]
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct InchiKey {
    /// The 27-character InChIKey.
    pub key: String,
    /// The first hash-extension block, if requested.
    pub extra1: Option<String>,
    /// The second hash-extension block, if requested.
    pub extra2: Option<String>,
}

/// Computes the InChIKey and, optionally, its two 256-bit hash extensions.
///
/// The hash extensions are an opt-in feature of the InChIKey algorithm that
/// reduce the (already tiny) collision probability of the 27-character key. Pass
/// `true` for `extra1`/`extra2` to compute the first/second extension block;
/// each requested block appears in the returned [`InchiKey`].
///
/// [`inchikey`] is the convenience form that returns just the key string.
///
/// # Errors
///
/// Returns [`InchiError::Failed`] if the input is not a valid InChI, and
/// [`InchiError::InteriorNul`] if it contains a NUL byte.
///
/// ```
/// use inchi::inchikey_with_hashes;
///
/// let k = inchikey_with_hashes("InChI=1S/CH4/h1H4", true, false)?;
/// assert_eq!(k.key, "VNWKTOKETHGBQD-UHFFFAOYSA-N");
/// assert!(k.extra1.is_some());
/// assert!(k.extra2.is_none());
/// # Ok::<(), inchi::InchiError>(())
/// ```
#[cfg(feature = "key")]
#[cfg_attr(docsrs, doc(cfg(feature = "key")))]
pub fn inchikey_with_hashes(
    inchi: impl AsRef<str>,
    extra1: bool,
    extra2: bool,
) -> Result<InchiKey> {
    use std::ffi::CStr;
    use std::os::raw::c_int;
    use std::ptr;

    let src = raw::to_cstring(inchi.as_ref())?;
    // The key buffer must be at least 28 bytes (27 characters + trailing NUL);
    // each extension buffer holds up to 64 characters + trailing NUL.
    let mut key = [0u8; 28];
    let mut buf1 = [0u8; 65];
    let mut buf2 = [0u8; 65];

    let ptr1 = if extra1 { buf1.as_mut_ptr() as *mut c_char } else { ptr::null_mut() };
    let ptr2 = if extra2 { buf2.as_mut_ptr() as *mut c_char } else { ptr::null_mut() };

    let _guard = raw::lock();
    // SAFETY: `src` is a valid NUL-terminated string; `key` is a 28-byte
    // writable buffer as required by the API; each extension pointer is either
    // null (extension disabled) or a 65-byte writable buffer.
    let rc = unsafe {
        inchi_sys::GetINCHIKeyFromINCHI(
            src.as_ptr(),
            c_int::from(extra1),
            c_int::from(extra2),
            key.as_mut_ptr() as *mut c_char,
            ptr1,
            ptr2,
        )
    };
    drop(src);

    // 0 == valid standard InChI, -1 == valid non-standard InChI; both succeed.
    if rc != inchi_sys::INCHIKEY_VALID_STANDARD && rc != inchi_sys::INCHIKEY_VALID_NON_STANDARD {
        let message = match rc {
            x if x == inchi_sys::INCHIKEY_INVALID_LENGTH => "InChI string has an invalid length",
            x if x == inchi_sys::INCHIKEY_INVALID_LAYOUT => "InChI string has an invalid layout",
            x if x == inchi_sys::INCHIKEY_INVALID_VERSION => "unsupported InChI version",
            _ => "could not compute InChIKey",
        };
        return Err(InchiError::Failed {
            status: Status::Error,
            message: message.to_string(),
        });
    }

    // SAFETY: on success the buffer holds a NUL-terminated ASCII InChIKey.
    let key = CStr::from_bytes_until_nul(&key).map_err(|_| InchiError::EmptyResult)?;
    Ok(InchiKey {
        key: key.to_string_lossy().into_owned(),
        extra1: extra1.then(|| read_hash_buf(&buf1)),
        extra2: extra2.then(|| read_hash_buf(&buf2)),
    })
}

/// Reads a NUL-terminated hash-extension buffer into an owned `String`.
#[cfg(feature = "key")]
fn read_hash_buf(buf: &[u8]) -> String {
    use std::ffi::CStr;
    CStr::from_bytes_until_nul(buf)
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default()
}

/// Shared post-processing for the `Get*`/`Make*` entry points: classifies the
/// return code and reads the populated [`OutputGuard`](raw::OutputGuard).
pub(crate) fn build_output(rc: i32, out: &raw::OutputGuard) -> Result<InchiOutput> {
    let status = Status::from_code(rc);
    if !status.is_success() {
        let mut message = out.message();
        if message.is_empty() {
            message = out.log();
        }
        return Err(InchiError::Failed { status, message });
    }

    let inchi = out.inchi();
    if inchi.is_empty() {
        return Err(InchiError::EmptyResult);
    }

    Ok(InchiOutput {
        inchi,
        aux_info: out.aux_info(),
        message: out.message(),
        log: out.log(),
        status,
    })
}
