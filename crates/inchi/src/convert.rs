//! InChI-to-InChI conversion via the native `GetINCHIfromINCHI` entry point.

use crate::error::Result;
use crate::options::Options;
use crate::output::InchiOutput;

/// Re-derives an InChI from an existing InChI string.
///
/// This is the programmatic form of the command-line `-InChI2InChI` option: it
/// parses the input InChI and regenerates it, optionally under different
/// [`Options`]. Common uses are normalizing a third-party InChI, or converting a
/// standard InChI into a non-standard one (e.g. with `FixedH` or `SUU`).
///
/// It is also the conversion the library performs internally for the strict
/// mode of [`check_inchi`](crate::check_inchi).
///
/// # Errors
///
/// Returns [`InchiError::Failed`](crate::InchiError::Failed) if the input is
/// not a usable InChI, [`InchiError::EmptyResult`](crate::InchiError::EmptyResult)
/// if conversion yields nothing, and
/// [`InchiError::InteriorNul`](crate::InchiError::InteriorNul) if either string
/// contains a NUL byte.
///
/// ```
/// use inchi::{inchi_to_inchi};
///
/// // Re-deriving a standard InChI with default options reproduces it.
/// let out = inchi_to_inchi("InChI=1S/CH4/h1H4", ())?;
/// assert_eq!(out.inchi(), "InChI=1S/CH4/h1H4");
/// # Ok::<(), inchi::InchiError>(())
/// ```
pub fn inchi_to_inchi(inchi: impl AsRef<str>, options: impl Into<Options>) -> Result<InchiOutput> {
    let options = options.into();
    let src = crate::raw::to_cstring(inchi.as_ref())?;
    let opts = crate::raw::to_cstring(&options.to_arg_string())?;

    let mut input = inchi_sys::inchi_InputINCHI {
        szInChI: src.as_ptr() as *mut std::os::raw::c_char,
        szOptions: opts.as_ptr() as *mut std::os::raw::c_char,
    };

    let _guard = crate::raw::lock();
    let mut out = crate::raw::OutputGuard::new();
    // SAFETY: `input` borrows `src`/`opts`, both alive until after the call;
    // `GetINCHIfromINCHI` does not take ownership, and `out` frees its result.
    let rc = unsafe { inchi_sys::GetINCHIfromINCHI(&mut input, out.as_mut_ptr()) };
    drop(src);
    drop(opts);

    // The InChI2InChI entry point reuses the standard `GetINCHI` return codes,
    // so the shared classifier handles success, warnings, and failures alike.
    crate::build_output(rc, &out)
}
