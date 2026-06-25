//! Error and status types for the [`crate`] API.

use std::fmt;

/// A specialized [`Result`](std::result::Result) for InChI operations.
///
/// ```
/// use inchi::{Result, inchikey};
///
/// fn key_for(inchi: &str) -> Result<String> {
///     inchikey(inchi)
/// }
/// # let _ = key_for("InChI=1S/CH4/h1H4");
/// ```
pub type Result<T> = std::result::Result<T, InchiError>;

/// The status class returned by the InChI C library.
///
/// These mirror the upstream `RetValGetINCHI` return codes. [`Status::Okay`]
/// and [`Status::Warning`] are *successful* outcomes (an identifier was
/// produced) and are therefore never wrapped in an [`InchiError`]; the
/// remaining variants always denote failure.
///
/// ```
/// use inchi::Status;
///
/// assert!(Status::Okay.is_success());
/// assert!(Status::Warning.is_success());
/// assert!(!Status::Error.is_success());
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Status {
    /// Success with no diagnostics (`inchi_Ret_OKAY`).
    Okay,
    /// Success, but the library issued a warning (`inchi_Ret_WARNING`). An
    /// identifier was still produced.
    Warning,
    /// Failure: no identifier was produced (`inchi_Ret_ERROR`).
    Error,
    /// Severe failure, typically a memory-allocation error (`inchi_Ret_FATAL`).
    Fatal,
    /// An unknown program error occurred (`inchi_Ret_UNKNOWN`).
    Unknown,
    /// A previous call had not yet returned (`inchi_Ret_BUSY`). The safe API
    /// serializes calls, so this should not normally be observed.
    Busy,
    /// No structural data was provided (`inchi_Ret_EOF`).
    Eof,
    /// A return code not covered by the documented set, preserved verbatim.
    Other(i32),
}

impl Status {
    /// Maps a raw C return code to a [`Status`].
    ///
    /// ```
    /// use inchi::Status;
    ///
    /// assert_eq!(Status::from_code(0), Status::Okay);
    /// assert_eq!(Status::from_code(2), Status::Error);
    /// assert_eq!(Status::from_code(42), Status::Other(42));
    /// ```
    #[must_use]
    pub fn from_code(code: i32) -> Self {
        match code {
            0 => Status::Okay,
            1 => Status::Warning,
            2 => Status::Error,
            3 => Status::Fatal,
            4 => Status::Unknown,
            5 => Status::Busy,
            -1 => Status::Eof,
            other => Status::Other(other),
        }
    }

    /// Returns `true` if an identifier was produced (`Okay` or `Warning`).
    ///
    /// ```
    /// use inchi::Status;
    /// assert!(Status::Okay.is_success());
    /// assert!(!Status::Fatal.is_success());
    /// ```
    #[must_use]
    pub fn is_success(self) -> bool {
        matches!(self, Status::Okay | Status::Warning)
    }
}

impl fmt::Display for Status {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Status::Okay => "okay",
            Status::Warning => "warning",
            Status::Error => "error",
            Status::Fatal => "fatal error",
            Status::Unknown => "unknown error",
            Status::Busy => "busy",
            Status::Eof => "no input",
            Status::Other(code) => return write!(f, "unrecognized status {code}"),
        };
        f.write_str(s)
    }
}

/// The error type returned by all fallible operations in this crate.
///
/// ```
/// use inchi::{from_molfile, InchiError};
///
/// // A molfile containing an interior NUL byte cannot be passed to C.
/// let err = from_molfile("a\0b", ()).unwrap_err();
/// assert!(matches!(err, InchiError::InteriorNul { .. }));
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum InchiError {
    /// The InChI library reported a failure. No identifier was produced.
    Failed {
        /// The failure class reported by the library.
        status: Status,
        /// The human-readable message from the library (may be empty).
        message: String,
    },
    /// The library returned success but produced an empty/absent identifier.
    EmptyResult,
    /// An input string contained an interior NUL byte at the given index and
    /// therefore could not be converted to a C string.
    InteriorNul {
        /// Byte offset of the first NUL in the input.
        position: usize,
    },
    /// The structure passed to a programmatic builder was invalid (e.g. a bond
    /// referencing a nonexistent atom, or an element symbol that is too long).
    InvalidStructure {
        /// A description of what was wrong.
        reason: String,
    },
}

impl InchiError {
    /// Returns the [`Status`] if this error originated from a library failure.
    ///
    /// ```
    /// use inchi::{InchiError, Status};
    ///
    /// let err = InchiError::Failed { status: Status::Error, message: "bad".into() };
    /// assert_eq!(err.status(), Some(Status::Error));
    /// assert_eq!(InchiError::EmptyResult.status(), None);
    /// ```
    #[must_use]
    pub fn status(&self) -> Option<Status> {
        match self {
            InchiError::Failed { status, .. } => Some(*status),
            _ => None,
        }
    }
}

impl fmt::Display for InchiError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InchiError::Failed { status, message } if message.is_empty() => {
                write!(f, "InChI generation failed: {status}")
            }
            InchiError::Failed { status, message } => {
                write!(f, "InChI generation failed: {status}: {message}")
            }
            InchiError::EmptyResult => {
                f.write_str("the InChI library produced an empty result")
            }
            InchiError::InteriorNul { position } => {
                write!(f, "input contains an interior NUL byte at position {position}")
            }
            InchiError::InvalidStructure { reason } => {
                write!(f, "invalid structure: {reason}")
            }
        }
    }
}

impl std::error::Error for InchiError {}
