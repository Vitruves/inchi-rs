//! The [`InchiOutput`] result type produced by InChI generation.

use crate::error::Status;

/// The successful result of generating an InChI.
///
/// An `InchiOutput` is returned when generation succeeds — including the case
/// where the library issued a warning ([`Status::Warning`]) but still produced
/// an identifier. Inspect [`InchiOutput::status`] and [`InchiOutput::message`]
/// to distinguish a clean result from a warned one.
///
/// ```
/// use inchi::{from_molfile, Status};
///
/// let methane = "\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
/// let out = from_molfile(methane, ())?;
/// assert_eq!(out.inchi(), "InChI=1S/CH4/h1H4");
/// assert_eq!(out.status(), Status::Okay);
/// # Ok::<(), inchi::InchiError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InchiOutput {
    pub(crate) inchi: String,
    pub(crate) aux_info: String,
    pub(crate) message: String,
    pub(crate) log: String,
    pub(crate) status: Status,
}

impl InchiOutput {
    /// The generated InChI identifier, e.g. `InChI=1S/CH4/h1H4`.
    ///
    /// ```
    /// # use inchi::{from_molfile};
    /// # let m = "\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
    /// let out = from_molfile(m, ())?;
    /// assert_eq!(out.inchi(), "InChI=1S/H2O/h1H2");
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn inchi(&self) -> &str {
        &self.inchi
    }

    /// The auxiliary information (`AuxInfo=...`) string, or empty if it was
    /// disabled via [`Options::aux_info(false)`](crate::Options::aux_info).
    ///
    /// ```
    /// # use inchi::{from_molfile};
    /// # let m = "\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
    /// let out = from_molfile(m, ())?;
    /// assert!(out.aux_info().starts_with("AuxInfo="));
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn aux_info(&self) -> &str {
        &self.aux_info
    }

    /// Any warning or informational message from the library (empty if none).
    ///
    /// ```
    /// # use inchi::{from_molfile};
    /// # let m = "\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
    /// let out = from_molfile(m, ())?;
    /// assert!(out.message().is_empty());
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn message(&self) -> &str {
        &self.message
    }

    /// The library's human-readable log of recognized options and diagnostics.
    ///
    /// ```
    /// # use inchi::{from_molfile};
    /// # let m = "\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
    /// let out = from_molfile(m, ())?;
    /// let _ = out.log();
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn log(&self) -> &str {
        &self.log
    }

    /// The status class of this (successful) result: [`Status::Okay`] or
    /// [`Status::Warning`].
    ///
    /// ```
    /// # use inchi::{from_molfile, Status};
    /// # let m = "\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
    /// let out = from_molfile(m, ())?;
    /// assert!(out.status().is_success());
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn status(&self) -> Status {
        self.status
    }

    /// Returns `true` if generation completed without any warning.
    ///
    /// ```
    /// # use inchi::{from_molfile};
    /// # let m = "\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
    /// let out = from_molfile(m, ())?;
    /// assert!(out.is_clean());
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn is_clean(&self) -> bool {
        self.status == Status::Okay
    }

    /// Consumes the result, returning the owned InChI string.
    ///
    /// ```
    /// # use inchi::{from_molfile};
    /// # let m = "\n  t\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
    /// let inchi: String = from_molfile(m, ())?.into_inchi();
    /// assert_eq!(inchi, "InChI=1S/CH4/h1H4");
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn into_inchi(self) -> String {
        self.inchi
    }
}
