//! Validation of existing InChI and InChIKey strings via the native
//! `CheckINCHI` / `CheckINCHIKey` entry points.

use crate::error::Result;
use std::os::raw::c_int;

/// The verdict returned by [`check_inchi`] for an InChI string.
///
/// The three `*Valid*` variants report success; [`InchiValidity::is_valid`]
/// collapses them to a boolean when the distinction does not matter.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum InchiValidity {
    /// A valid *standard* InChI (`INCHI_VALID_STANDARD`).
    Standard,
    /// A valid *non-standard* InChI (`INCHI_VALID_NON_STANDARD`).
    NonStandard,
    /// A valid InChI produced by a beta/experimental feature
    /// (`INCHI_VALID_BETA`).
    Beta,
    /// The `InChI=` prefix is missing or malformed (`INCHI_INVALID_PREFIX`).
    InvalidPrefix,
    /// The version is unsupported (`INCHI_INVALID_VERSION`).
    InvalidVersion,
    /// The layout is malformed (`INCHI_INVALID_LAYOUT`).
    InvalidLayout,
    /// Strict checking re-derived a different InChI, so the input does not
    /// round-trip (`INCHI_FAIL_I2I`; only possible when `strict` is set).
    FailedRoundtrip,
    /// A return code outside the documented set, preserved verbatim.
    Other(i32),
}

impl InchiValidity {
    /// Returns `true` for the three valid verdicts.
    ///
    /// ```
    /// use inchi::{check_inchi, InchiValidity};
    ///
    /// assert!(check_inchi("InChI=1S/CH4/h1H4", false)?.is_valid());
    /// assert!(!check_inchi("not an inchi", false)?.is_valid());
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn is_valid(self) -> bool {
        matches!(self, InchiValidity::Standard | InchiValidity::NonStandard | InchiValidity::Beta)
    }
}

/// The verdict returned by [`check_inchikey`] for an InChIKey string.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum InchiKeyValidity {
    /// A valid key derived from a *standard* InChI (`INCHIKEY_VALID_STANDARD`).
    Standard,
    /// A valid key derived from a *non-standard* InChI
    /// (`INCHIKEY_VALID_NON_STANDARD`).
    NonStandard,
    /// The key has the wrong length (`INCHIKEY_INVALID_LENGTH`).
    InvalidLength,
    /// The key layout is malformed (`INCHIKEY_INVALID_LAYOUT`).
    InvalidLayout,
    /// The key encodes an unsupported version (`INCHIKEY_INVALID_VERSION`).
    InvalidVersion,
    /// A return code outside the documented set, preserved verbatim.
    Other(i32),
}

impl InchiKeyValidity {
    /// Returns `true` for the two valid verdicts.
    ///
    /// ```
    /// use inchi::check_inchikey;
    ///
    /// assert!(check_inchikey("VNWKTOKETHGBQD-UHFFFAOYSA-N")?.is_valid());
    /// assert!(!check_inchikey("definitely-not-a-key")?.is_valid());
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn is_valid(self) -> bool {
        matches!(self, InchiKeyValidity::Standard | InchiKeyValidity::NonStandard)
    }
}

/// Checks whether a string is a valid InChI.
///
/// This is a faithful binding to the reference library's `CheckINCHI`, including
/// its two distinct modes:
///
/// * **Lenient** (`strict = false`) — the everyday validity check. It inspects
///   the prefix, version, and overall layout and classifies the string as
///   [`Standard`](InchiValidity::Standard),
///   [`NonStandard`](InchiValidity::NonStandard),
///   [`Beta`](InchiValidity::Beta), or one of the `Invalid*` verdicts. Use this
///   unless you specifically need the round-trip check below.
///
/// * **Strict** (`strict = true`) — additionally re-derives the InChI through a
///   canonical `FixedH RecMet SUU SLUUD` conversion and requires the result to
///   match the input *byte for byte*. Because that conversion always emits a
///   **non-standard** InChI (`InChI=1/…`), a **standard** InChI (`InChI=1S/…`)
///   can never match and is reported as
///   [`FailedRoundtrip`](InchiValidity::FailedRoundtrip). Strict mode is
///   therefore meaningful only for non-standard InChIs produced with that
///   option set; for standard InChIs, lenient mode is the correct validator.
///   This mirrors the IUPAC reference implementation exactly.
///
/// An *invalid* InChI is reported through [`InchiValidity`], not as an [`Err`];
/// only an interior NUL byte in the input produces an error.
///
/// # Errors
///
/// Returns [`InchiError::InteriorNul`](crate::InchiError::InteriorNul) if
/// `inchi` contains a NUL byte.
///
/// ```
/// use inchi::{check_inchi, InchiValidity};
///
/// // Lenient mode is the general-purpose validity check.
/// assert_eq!(check_inchi("InChI=1S/H2O/h1H2", false)?, InchiValidity::Standard);
/// assert_eq!(check_inchi("InChI=1/H2O/h1H2", false)?, InchiValidity::NonStandard);
///
/// // Strict mode validates a round-tripping non-standard InChI...
/// assert!(check_inchi("InChI=1/H2O/h1H2", true)?.is_valid());
/// // ...but reports a standard InChI as FailedRoundtrip (by reference design).
/// assert_eq!(check_inchi("InChI=1S/H2O/h1H2", true)?, InchiValidity::FailedRoundtrip);
/// # Ok::<(), inchi::InchiError>(())
/// ```
pub fn check_inchi(inchi: impl AsRef<str>, strict: bool) -> Result<InchiValidity> {
    let src = crate::raw::to_cstring(inchi.as_ref())?;

    let _guard = crate::raw::lock();
    // SAFETY: `src` is a valid NUL-terminated string for the duration of the
    // call; `CheckINCHI` only reads it. Serialized by `_guard`.
    let rc = unsafe { inchi_sys::CheckINCHI(src.as_ptr(), c_int::from(strict)) };
    drop(src);

    Ok(match rc {
        x if x == inchi_sys::INCHI_VALID_STANDARD as i32 => InchiValidity::Standard,
        x if x == inchi_sys::INCHI_VALID_NON_STANDARD as i32 => InchiValidity::NonStandard,
        x if x == inchi_sys::INCHI_VALID_BETA as i32 => InchiValidity::Beta,
        x if x == inchi_sys::INCHI_INVALID_PREFIX as i32 => InchiValidity::InvalidPrefix,
        x if x == inchi_sys::INCHI_INVALID_VERSION as i32 => InchiValidity::InvalidVersion,
        x if x == inchi_sys::INCHI_INVALID_LAYOUT as i32 => InchiValidity::InvalidLayout,
        x if x == inchi_sys::INCHI_FAIL_I2I as i32 => InchiValidity::FailedRoundtrip,
        other => InchiValidity::Other(other),
    })
}

/// Checks whether a string is a syntactically valid InChIKey.
///
/// This validates the 27-character layout and its check character; it does
/// **not** prove that a key corresponds to any real structure (the hash is not
/// invertible). As with [`check_inchi`], an invalid key is reported through the
/// returned [`InchiKeyValidity`] rather than as an error.
///
/// # Errors
///
/// Returns [`InchiError::InteriorNul`](crate::InchiError::InteriorNul) if `key`
/// contains a NUL byte.
///
/// ```
/// use inchi::{check_inchikey, InchiKeyValidity};
///
/// assert_eq!(check_inchikey("VNWKTOKETHGBQD-UHFFFAOYSA-N")?, InchiKeyValidity::Standard);
/// assert_eq!(check_inchikey("TOOSHORT")?, InchiKeyValidity::InvalidLength);
/// # Ok::<(), inchi::InchiError>(())
/// ```
pub fn check_inchikey(key: impl AsRef<str>) -> Result<InchiKeyValidity> {
    let src = crate::raw::to_cstring(key.as_ref())?;

    let _guard = crate::raw::lock();
    // SAFETY: `src` is a valid NUL-terminated string the call only reads.
    let rc = unsafe { inchi_sys::CheckINCHIKey(src.as_ptr()) };
    drop(src);

    Ok(match rc {
        x if x == inchi_sys::INCHIKEY_VALID_STANDARD => InchiKeyValidity::Standard,
        x if x == inchi_sys::INCHIKEY_VALID_NON_STANDARD => InchiKeyValidity::NonStandard,
        x if x == inchi_sys::INCHIKEY_INVALID_LENGTH => InchiKeyValidity::InvalidLength,
        x if x == inchi_sys::INCHIKEY_INVALID_LAYOUT => InchiKeyValidity::InvalidLayout,
        x if x == inchi_sys::INCHIKEY_INVALID_VERSION => InchiKeyValidity::InvalidVersion,
        other => InchiKeyValidity::Other(other),
    })
}
