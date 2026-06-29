//! Queryable molecules via the InChI eXtensible API (IXA).
//!
//! The classic API ([`from_molfile`](crate::from_molfile),
//! [`struct_from_inchi`](crate::struct_from_inchi)) is one-shot: structure in,
//! identifier or fixed-array structure out. [`IxaMolecule`] instead keeps a
//! *live* molecule object you can interrogate atom-by-atom and bond-by-bond,
//! then turn into an InChI or InChIKey on demand.
//!
//! Two entry points populate a molecule:
//! - [`IxaMolecule::from_molfile`] parses a Molfile (V2000/V3000) text record.
//! - [`IxaMolecule::from_inchi`] reconstructs atoms and bonds from an InChI.
//!
//! ```
//! use inchi::IxaMolecule;
//!
//! let mol = IxaMolecule::from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")?;
//! assert_eq!(mol.atom_count(), 3); // C, C, O (hydrogens are implicit)
//! assert_eq!(mol.to_inchi(())?.inchi(), "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3");
//! # Ok::<(), inchi::InchiError>(())
//! ```

use crate::error::{InchiError, Result, Status};
use crate::molecule::{BondOrder, Parity, Radical};
use crate::options::Options;
use crate::output::InchiOutput;
use crate::raw;
use std::os::raw::{c_int, c_long};

/// An atom identifier returned by an [`IxaMolecule`].
///
/// IDs are only meaningful for the molecule that produced them; using one with
/// a different molecule has unspecified results. `AtomId` is a cheap [`Copy`]
/// handle that borrows nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomId(inchi_sys::IXA_ATOMID);

/// A bond identifier returned by an [`IxaMolecule`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BondId(inchi_sys::IXA_BONDID);

/// A stereo-descriptor identifier returned by an [`IxaMolecule`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StereoId(inchi_sys::IXA_STEREOID);

// SAFETY: these newtypes wrap opaque IXA identifiers that are, at the C level,
// integers cast to pointers; they carry no pointee and are never dereferenced
// on the Rust side. Sending the value between threads is sound — it is only
// ever used (under the global lock) with the molecule that minted it.
unsafe impl Send for AtomId {}
unsafe impl Sync for AtomId {}
unsafe impl Send for BondId {}
unsafe impl Sync for BondId {}
unsafe impl Send for StereoId {}
unsafe impl Sync for StereoId {}

/// Owns the paired IXA status + molecule handles and frees them on drop.
struct IxaHandles {
    status: inchi_sys::IXA_STATUS_HANDLE,
    mol: inchi_sys::IXA_MOL_HANDLE,
}

impl Drop for IxaHandles {
    fn drop(&mut self) {
        let _guard = raw::lock();
        // SAFETY: both handles were created non-null in the constructor and are
        // not freed elsewhere; the molecule is destroyed before its status, as
        // the IXA API requires. Serialized by the global lock.
        unsafe {
            inchi_sys::IXA_MOL_Destroy(self.status, self.mol);
            inchi_sys::IXA_STATUS_Destroy(self.status);
        }
    }
}

// SAFETY: `IxaHandles` holds opaque C heap pointers valid for the lifetime of
// the value; every access goes through the global `INCHI_LOCK`, so the handles
// are never touched concurrently.
unsafe impl Send for IxaHandles {}
unsafe impl Sync for IxaHandles {}

/// A molecule loaded from a Molfile or InChI string, queryable atom-by-atom and
/// convertible back into an InChI or InChIKey.
///
/// See the [module documentation](self) for an overview.
pub struct IxaMolecule(IxaHandles);

impl IxaMolecule {
    /// Parses a Molfile (V2000/V3000) text record into a queryable molecule.
    ///
    /// # Errors
    ///
    /// Returns [`InchiError::Failed`] if the Molfile cannot be parsed and
    /// [`InchiError::InteriorNul`] if `molfile` contains a NUL byte.
    ///
    /// ```
    /// use inchi::IxaMolecule;
    ///
    /// let methane = include_str!("../tests/fixtures/methane.mol");
    /// let mol = IxaMolecule::from_molfile(methane)?;
    /// assert_eq!(mol.atom_count(), 1);
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    pub fn from_molfile(molfile: impl AsRef<str>) -> Result<Self> {
        let text = raw::to_cstring(molfile.as_ref())?;
        // SAFETY: `read` is invoked under the lock and `IXA_MOL_ReadMolfile`
        // takes a valid NUL-terminated pointer that outlives the call.
        Self::build(|status, mol| unsafe {
            inchi_sys::IXA_MOL_ReadMolfile(status, mol, text.as_ptr());
        })
    }

    /// Reconstructs a molecule from an InChI string.
    ///
    /// Hydrogens remain implicit and all atom coordinates are zero — an InChI
    /// stores neither explicit hydrogens nor geometry.
    ///
    /// # Errors
    ///
    /// Returns [`InchiError::Failed`] if the InChI cannot be parsed and
    /// [`InchiError::InteriorNul`] if `inchi` contains a NUL byte.
    ///
    /// ```
    /// use inchi::IxaMolecule;
    ///
    /// let mol = IxaMolecule::from_inchi("InChI=1S/CH4/h1H4")?;
    /// assert_eq!(mol.atom_count(), 1);
    /// assert!(IxaMolecule::from_inchi("not an inchi").is_err());
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    pub fn from_inchi(inchi: impl AsRef<str>) -> Result<Self> {
        let text = raw::to_cstring(inchi.as_ref())?;
        // SAFETY: as above; `IXA_MOL_ReadInChI` reads from a valid pointer.
        Self::build(|status, mol| unsafe {
            inchi_sys::IXA_MOL_ReadInChI(status, mol, text.as_ptr());
        })
    }

    /// Creates the status/molecule handles, runs `read`, and checks for errors.
    fn build(
        read: impl FnOnce(inchi_sys::IXA_STATUS_HANDLE, inchi_sys::IXA_MOL_HANDLE),
    ) -> Result<Self> {
        let _guard = raw::lock();
        // SAFETY: `IXA_STATUS_Create`/`IXA_MOL_Create` return owned handles (or
        // null on allocation failure, which we detect). On any early return the
        // already-created handles are freed before returning the error.
        unsafe {
            let status = inchi_sys::IXA_STATUS_Create();
            if status.is_null() {
                return Err(InchiError::Failed {
                    status: Status::Fatal,
                    message: "could not allocate IXA status object".to_string(),
                });
            }
            let mol = inchi_sys::IXA_MOL_Create(status);
            if mol.is_null() {
                inchi_sys::IXA_STATUS_Destroy(status);
                return Err(InchiError::Failed {
                    status: Status::Fatal,
                    message: "could not allocate IXA molecule object".to_string(),
                });
            }

            read(status, mol);

            if let Some(err) = take_error(status) {
                inchi_sys::IXA_MOL_Destroy(status, mol);
                inchi_sys::IXA_STATUS_Destroy(status);
                return Err(err);
            }

            // Some malformed inputs (e.g. text lacking the `InChI=` prefix) are
            // rejected by the reader without pushing an error, leaving an empty
            // molecule behind. A structure with no atoms is never useful, so we
            // treat it as a parse failure — matching the classic API.
            if inchi_sys::IXA_MOL_GetNumAtoms(status, mol) <= 0 {
                inchi_sys::IXA_MOL_Destroy(status, mol);
                inchi_sys::IXA_STATUS_Destroy(status);
                return Err(InchiError::Failed {
                    status: Status::Error,
                    message: "input did not yield any atoms".to_string(),
                });
            }

            Ok(IxaMolecule(IxaHandles { status, mol }))
        }
    }

    /// The number of (heavy) atoms in the molecule.
    #[must_use]
    pub fn atom_count(&self) -> usize {
        self.count(inchi_sys::IXA_MOL_GetNumAtoms)
    }

    /// The number of bonds in the molecule.
    #[must_use]
    pub fn bond_count(&self) -> usize {
        self.count(inchi_sys::IXA_MOL_GetNumBonds)
    }

    /// The number of stereo descriptors in the molecule.
    #[must_use]
    pub fn stereo_count(&self) -> usize {
        self.count(inchi_sys::IXA_MOL_GetNumStereos)
    }

    /// Whether the molecule carries a set chiral flag.
    ///
    /// ```
    /// # use inchi::IxaMolecule;
    /// let mol = IxaMolecule::from_inchi("InChI=1S/CH4/h1H4")?;
    /// assert!(!mol.is_chiral());
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn is_chiral(&self) -> bool {
        let _guard = raw::lock();
        // SAFETY: handles are valid; serialized by the lock.
        unsafe {
            inchi_sys::IXA_STATUS_Clear(self.0.status);
            inchi_sys::IXA_MOL_GetChiral(self.0.status, self.0.mol) == inchi_sys::IXA_TRUE
        }
    }

    /// Returns every atom identifier, in index order.
    ///
    /// ```
    /// # use inchi::IxaMolecule;
    /// let mol = IxaMolecule::from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")?;
    /// assert_eq!(mol.atoms().count(), 3);
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    pub fn atoms(&self) -> impl Iterator<Item = AtomId> {
        let n = self.atom_count();
        let _guard = raw::lock();
        let mut ids = Vec::with_capacity(n);
        // SAFETY: indices `0..n` are in range for `IXA_MOL_GetAtomId`; the lock
        // serializes the calls.
        unsafe {
            inchi_sys::IXA_STATUS_Clear(self.0.status);
            for i in 0..n {
                let Ok(idx) = c_int::try_from(i) else { break };
                ids.push(AtomId(inchi_sys::IXA_MOL_GetAtomId(
                    self.0.status,
                    self.0.mol,
                    idx,
                )));
            }
        }
        ids.into_iter()
    }

    /// Returns every bond identifier, in index order.
    pub fn bonds(&self) -> impl Iterator<Item = BondId> {
        let n = self.bond_count();
        let _guard = raw::lock();
        let mut ids = Vec::with_capacity(n);
        // SAFETY: indices `0..n` are in range for `IXA_MOL_GetBondId`.
        unsafe {
            inchi_sys::IXA_STATUS_Clear(self.0.status);
            for i in 0..n {
                let Ok(idx) = c_int::try_from(i) else { break };
                ids.push(BondId(inchi_sys::IXA_MOL_GetBondId(
                    self.0.status,
                    self.0.mol,
                    idx,
                )));
            }
        }
        ids.into_iter()
    }

    /// Returns every stereo-descriptor identifier, in index order.
    pub fn stereos(&self) -> impl Iterator<Item = StereoId> {
        let n = self.stereo_count();
        let _guard = raw::lock();
        let mut ids = Vec::with_capacity(n);
        // SAFETY: indices `0..n` are in range for `IXA_MOL_GetStereoId`.
        unsafe {
            inchi_sys::IXA_STATUS_Clear(self.0.status);
            for i in 0..n {
                let Ok(idx) = c_int::try_from(i) else { break };
                ids.push(StereoId(inchi_sys::IXA_MOL_GetStereoId(
                    self.0.status,
                    self.0.mol,
                    idx,
                )));
            }
        }
        ids.into_iter()
    }

    /// The element symbol of an atom (e.g. `"C"`, `"Fe"`).
    ///
    /// ```
    /// # use inchi::IxaMolecule;
    /// let mol = IxaMolecule::from_inchi("InChI=1S/H2O/h1H2")?;
    /// let oxygen = mol.atoms().next().unwrap();
    /// assert_eq!(mol.atom_element(oxygen)?, "O");
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    pub fn atom_element(&self, atom: AtomId) -> Result<String> {
        self.checked(|status, mol| {
            // SAFETY: `atom` is an IXA id; the returned pointer is read before
            // the lock is released and any error is reported via the status.
            unsafe { raw::cstr_to_string(inchi_sys::IXA_MOL_GetAtomElement(status, mol, atom.0)) }
        })
    }

    /// The atomic number of an atom.
    pub fn atom_atomic_number(&self, atom: AtomId) -> Result<u32> {
        self.checked(|status, mol| unsafe {
            inchi_sys::IXA_MOL_GetAtomAtomicNumber(status, mol, atom.0)
        })
        .map(|n| u32::try_from(n).unwrap_or(0))
    }

    /// The mass number of an atom (`0` for the natural isotopic composition).
    pub fn atom_mass(&self, atom: AtomId) -> Result<i32> {
        self.checked(|status, mol| unsafe { inchi_sys::IXA_MOL_GetAtomMass(status, mol, atom.0) })
    }

    /// The formal charge of an atom.
    pub fn atom_charge(&self, atom: AtomId) -> Result<i32> {
        self.checked(|status, mol| unsafe { inchi_sys::IXA_MOL_GetAtomCharge(status, mol, atom.0) })
    }

    /// The radical state of an atom.
    pub fn atom_radical(&self, atom: AtomId) -> Result<Radical> {
        self.checked(|status, mol| unsafe {
            inchi_sys::IXA_MOL_GetAtomRadical(status, mol, atom.0)
        })
        .map(decode_radical)
    }

    /// The hydrogen counts attached to an atom, indexed by isotope mass number:
    /// `[implicit, ¹H, ²H, ³H]`. Index `0` is the count of ordinary
    /// (non-isotopic) implicit hydrogens.
    ///
    /// ```
    /// # use inchi::IxaMolecule;
    /// let mol = IxaMolecule::from_inchi("InChI=1S/CH4/h1H4")?;
    /// let carbon = mol.atoms().next().unwrap();
    /// assert_eq!(mol.atom_hydrogens(carbon)?[0], 4);
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    pub fn atom_hydrogens(&self, atom: AtomId) -> Result<[i32; 4]> {
        self.checked(|status, mol| unsafe {
            [
                inchi_sys::IXA_MOL_GetAtomHydrogens(status, mol, atom.0, 0),
                inchi_sys::IXA_MOL_GetAtomHydrogens(status, mol, atom.0, 1),
                inchi_sys::IXA_MOL_GetAtomHydrogens(status, mol, atom.0, 2),
                inchi_sys::IXA_MOL_GetAtomHydrogens(status, mol, atom.0, 3),
            ]
        })
    }

    /// The Cartesian coordinates of an atom (all zero for InChI-derived molecules).
    pub fn atom_position(&self, atom: AtomId) -> Result<(f64, f64, f64)> {
        self.checked(|status, mol| unsafe {
            (
                inchi_sys::IXA_MOL_GetAtomX(status, mol, atom.0),
                inchi_sys::IXA_MOL_GetAtomY(status, mol, atom.0),
                inchi_sys::IXA_MOL_GetAtomZ(status, mol, atom.0),
            )
        })
    }

    /// The number of bonds attached to an atom.
    pub fn atom_bond_count(&self, atom: AtomId) -> Result<usize> {
        self.checked(|status, mol| unsafe {
            inchi_sys::IXA_MOL_GetAtomNumBonds(status, mol, atom.0)
        })
        .map(|n| usize::try_from(n).unwrap_or(0))
    }

    /// The bond connecting `atom` to its `n`th neighbor.
    pub fn atom_bond(&self, atom: AtomId, n: usize) -> Result<BondId> {
        let idx = c_int::try_from(n).map_err(|_| InchiError::InvalidStructure {
            reason: format!("bond index {n} out of range"),
        })?;
        self.checked(|status, mol| unsafe {
            BondId(inchi_sys::IXA_MOL_GetAtomBond(status, mol, atom.0, idx))
        })
    }

    /// The two atoms at the ends of a bond.
    pub fn bond_atoms(&self, bond: BondId) -> Result<(AtomId, AtomId)> {
        self.checked(|status, mol| unsafe {
            (
                AtomId(inchi_sys::IXA_MOL_GetBondAtom1(status, mol, bond.0)),
                AtomId(inchi_sys::IXA_MOL_GetBondAtom2(status, mol, bond.0)),
            )
        })
    }

    /// The order of a bond.
    pub fn bond_order(&self, bond: BondId) -> Result<BondOrder> {
        self.checked(|status, mol| unsafe { inchi_sys::IXA_MOL_GetBondType(status, mol, bond.0) })
            .map(decode_bond_type)
    }

    /// The parity of a stereo descriptor, or `None` if it carries no parity.
    pub fn stereo_parity(&self, stereo: StereoId) -> Result<Option<Parity>> {
        self.checked(|status, mol| unsafe {
            inchi_sys::IXA_MOL_GetStereoParity(status, mol, stereo.0)
        })
        .map(decode_parity)
    }

    /// The central atom of a stereo descriptor, or `None` for a descriptor
    /// (such as a stereogenic double bond) centered on a bond rather than an atom.
    pub fn stereo_central_atom(&self, stereo: StereoId) -> Result<Option<AtomId>> {
        self.checked(|status, mol| unsafe {
            let id = inchi_sys::IXA_MOL_GetStereoCentralAtom(status, mol, stereo.0);
            if id.is_null() {
                None
            } else {
                Some(AtomId(id))
            }
        })
    }

    /// Generates an InChI from the current molecule state.
    ///
    /// Options are honored through the IXA InChI builder. The `raw()` option
    /// escape hatch is *not* applied here (it has no IXA builder equivalent);
    /// use [`Molecule::to_inchi`](crate::Molecule::to_inchi) or
    /// [`from_molfile`](crate::from_molfile) if you need raw option tokens.
    ///
    /// ```
    /// # use inchi::IxaMolecule;
    /// let mol = IxaMolecule::from_inchi("InChI=1S/H2O/h1H2")?;
    /// assert_eq!(mol.to_inchi(())?.inchi(), "InChI=1S/H2O/h1H2");
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    pub fn to_inchi(&self, options: impl Into<Options>) -> Result<InchiOutput> {
        let options = options.into();
        let tokens = options.to_arg_string();
        let polymers = tokens.split_whitespace().any(|t| {
            let bare = t.trim_start_matches(['-', '/']);
            bare == "Polymers" || bare.starts_with("Polymers")
        });

        let _guard = raw::lock();
        // SAFETY: the builder is created, used, and destroyed within this locked
        // block; every pointer read from it is copied into an owned `String`
        // before `IXA_INCHIBUILDER_Destroy` invalidates it.
        unsafe {
            inchi_sys::IXA_STATUS_Clear(self.0.status);
            let builder = inchi_sys::IXA_INCHIBUILDER_Create(self.0.status);
            if builder.is_null() {
                return Err(InchiError::Failed {
                    status: Status::Fatal,
                    message: "could not allocate IXA InChI builder".to_string(),
                });
            }
            inchi_sys::IXA_INCHIBUILDER_SetMolecule(self.0.status, builder, self.0.mol);
            apply_options(self.0.status, builder, self.0.mol, &tokens);

            let inchi_ptr = if polymers {
                inchi_sys::IXA_INCHIBUILDER_GetInChIEx(self.0.status, builder)
            } else {
                inchi_sys::IXA_INCHIBUILDER_GetInChI(self.0.status, builder)
            };
            let has_error = inchi_sys::IXA_STATUS_HasError(self.0.status) == inchi_sys::IXA_TRUE;
            let has_warning =
                inchi_sys::IXA_STATUS_HasWarning(self.0.status) == inchi_sys::IXA_TRUE;
            let inchi = raw::cstr_to_string(inchi_ptr);
            let aux_info = raw::cstr_to_string(inchi_sys::IXA_INCHIBUILDER_GetAuxInfo(
                self.0.status,
                builder,
            ));
            let log =
                raw::cstr_to_string(inchi_sys::IXA_INCHIBUILDER_GetLog(self.0.status, builder));
            let message = collect_messages(self.0.status);

            inchi_sys::IXA_INCHIBUILDER_Destroy(self.0.status, builder);

            if has_error {
                return Err(InchiError::Failed {
                    status: Status::Error,
                    message: if message.is_empty() {
                        "IXA InChI generation failed".to_string()
                    } else {
                        message
                    },
                });
            }
            if inchi.is_empty() {
                return Err(InchiError::EmptyResult);
            }

            Ok(InchiOutput {
                inchi,
                aux_info,
                message,
                log,
                status: if has_warning {
                    Status::Warning
                } else {
                    Status::Okay
                },
            })
        }
    }

    /// Generates an InChIKey from the current molecule state.
    ///
    /// ```
    /// # use inchi::IxaMolecule;
    /// let mol = IxaMolecule::from_inchi("InChI=1S/CH4/h1H4")?;
    /// assert_eq!(mol.to_inchikey(())?, "VNWKTOKETHGBQD-UHFFFAOYSA-N");
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[cfg(feature = "key")]
    #[cfg_attr(docsrs, doc(cfg(feature = "key")))]
    pub fn to_inchikey(&self, options: impl Into<Options>) -> Result<String> {
        let inchi = self.to_inchi(options)?.into_inchi();
        let src = raw::to_cstring(&inchi)?;

        let _guard = raw::lock();
        // SAFETY: the key builder lifecycle is fully contained here; the key
        // pointer is copied into a `String` before the builder is destroyed.
        unsafe {
            inchi_sys::IXA_STATUS_Clear(self.0.status);
            let kb = inchi_sys::IXA_INCHIKEYBUILDER_Create(self.0.status);
            if kb.is_null() {
                return Err(InchiError::Failed {
                    status: Status::Fatal,
                    message: "could not allocate IXA InChIKey builder".to_string(),
                });
            }
            inchi_sys::IXA_INCHIKEYBUILDER_SetInChI(self.0.status, kb, src.as_ptr());
            let key = raw::cstr_to_string(inchi_sys::IXA_INCHIKEYBUILDER_GetInChIKey(
                self.0.status,
                kb,
            ));
            let message = collect_messages(self.0.status);
            let has_error = inchi_sys::IXA_STATUS_HasError(self.0.status) == inchi_sys::IXA_TRUE;
            inchi_sys::IXA_INCHIKEYBUILDER_Destroy(self.0.status, kb);

            if has_error || key.is_empty() {
                return Err(InchiError::Failed {
                    status: Status::Error,
                    message: if message.is_empty() {
                        "IXA InChIKey generation failed".to_string()
                    } else {
                        message
                    },
                });
            }
            Ok(key)
        }
    }

    /// Runs a count-returning IXA call, mapping a negative/failed result to `0`.
    fn count(
        &self,
        f: unsafe extern "C" fn(inchi_sys::IXA_STATUS_HANDLE, inchi_sys::IXA_MOL_HANDLE) -> c_int,
    ) -> usize {
        let _guard = raw::lock();
        // SAFETY: handles are valid; serialized by the lock.
        let n = unsafe {
            inchi_sys::IXA_STATUS_Clear(self.0.status);
            f(self.0.status, self.0.mol)
        };
        usize::try_from(n).unwrap_or(0)
    }

    /// Runs `f` under the lock with a freshly cleared status, then converts any
    /// reported IXA error into an [`InchiError`].
    fn checked<T>(
        &self,
        f: impl FnOnce(inchi_sys::IXA_STATUS_HANDLE, inchi_sys::IXA_MOL_HANDLE) -> T,
    ) -> Result<T> {
        let _guard = raw::lock();
        // SAFETY: handles are valid for the molecule's lifetime; the lock
        // serializes access. `f` performs the actual FFI calls.
        unsafe { inchi_sys::IXA_STATUS_Clear(self.0.status) };
        let out = f(self.0.status, self.0.mol);
        // SAFETY: the status handle is valid and just-populated by `f`.
        match unsafe { take_error(self.0.status) } {
            Some(err) => Err(err),
            None => Ok(out),
        }
    }
}

/// Applies an option string (as rendered by [`Options::to_arg_string`]) to an
/// IXA InChI builder. Tokens with no IXA builder equivalent are ignored.
///
/// # Safety
///
/// `status`, `builder`, and `mol` must be valid handles, and the caller must
/// hold the global lock for the duration of the call.
unsafe fn apply_options(
    status: inchi_sys::IXA_STATUS_HANDLE,
    builder: inchi_sys::IXA_INCHIBUILDER_HANDLE,
    mol: inchi_sys::IXA_MOL_HANDLE,
    tokens: &str,
) {
    use inchi_sys::*;
    let set = |opt: IXA_INCHIBUILDER_OPTION| {
        IXA_INCHIBUILDER_SetOption(status, builder, opt, IXA_TRUE);
    };
    let set_stereo = |opt: IXA_INCHIBUILDER_STEREOOPTION| {
        IXA_INCHIBUILDER_SetOption_Stereo(status, builder, opt);
    };

    for token in tokens.split_whitespace() {
        let bare = token.trim_start_matches(['-', '/']);
        match bare {
            "DoNotAddH" => set(IXA_INCHIBUILDER_OPTION_DoNotAddH),
            "AuxNone" => set(IXA_INCHIBUILDER_OPTION_AuxNone),
            "SNon" => set_stereo(IXA_INCHIBUILDER_STEREOOPTION_SNon),
            "SRel" => set_stereo(IXA_INCHIBUILDER_STEREOOPTION_SRel),
            "SRac" => set_stereo(IXA_INCHIBUILDER_STEREOOPTION_SRac),
            "SUCF" => set_stereo(IXA_INCHIBUILDER_STEREOOPTION_SUCF),
            "FixedH" => set(IXA_INCHIBUILDER_OPTION_FixedH),
            "RecMet" => set(IXA_INCHIBUILDER_OPTION_RecMet),
            "KET" => set(IXA_INCHIBUILDER_OPTION_KET),
            "15T" => set(IXA_INCHIBUILDER_OPTION_15T),
            "SaveOpt" => set(IXA_INCHIBUILDER_OPTION_SaveOpt),
            "Polymers" => set(IXA_INCHIBUILDER_OPTION_Polymers),
            "Polymers105" => set(IXA_INCHIBUILDER_OPTION_Polymers105),
            "NoFrameShift" => set(IXA_INCHIBUILDER_OPTION_NoFrameShift),
            "FoldSRU" => set(IXA_INCHIBUILDER_OPTION_FoldCRU),
            "NoEdits" => set(IXA_INCHIBUILDER_OPTION_NoEdits),
            "NPZz" => set(IXA_INCHIBUILDER_OPTION_NPZZ),
            "ChiralFlagON" => IXA_MOL_SetChiral(status, mol, IXA_TRUE),
            "ChiralFlagOFF" => IXA_MOL_SetChiral(status, mol, IXA_FALSE),
            other => {
                if let Some(ms) = other.strip_prefix("WM") {
                    if let Ok(value) = ms.parse::<c_long>() {
                        IXA_INCHIBUILDER_SetOption_Timeout_MilliSeconds(status, builder, value);
                    }
                }
                // Any other token has no IXA builder equivalent; ignore it.
            }
        }
    }
}

/// Collects every status message into a single `; `-joined string.
///
/// # Safety
///
/// `status` must be a valid handle and the caller must hold the global lock.
unsafe fn collect_messages(status: inchi_sys::IXA_STATUS_HANDLE) -> String {
    let count = inchi_sys::IXA_STATUS_GetCount(status);
    let mut parts = Vec::new();
    let mut i = 0;
    while i < count {
        let msg = raw::cstr_to_string(inchi_sys::IXA_STATUS_GetMessage(status, i));
        if !msg.is_empty() {
            parts.push(msg);
        }
        i += 1;
    }
    parts.join("; ")
}

/// Reads the error state of `status`, returning a populated [`InchiError`] if an
/// error (not merely a warning) was recorded.
///
/// # Safety
///
/// `status` must be a valid handle and the caller must hold the global lock.
unsafe fn take_error(status: inchi_sys::IXA_STATUS_HANDLE) -> Option<InchiError> {
    if inchi_sys::IXA_STATUS_HasError(status) != inchi_sys::IXA_TRUE {
        return None;
    }
    let message = collect_messages(status);
    Some(InchiError::Failed {
        status: Status::Error,
        message: if message.is_empty() {
            "IXA operation failed".to_string()
        } else {
            message
        },
    })
}

fn decode_radical(code: inchi_sys::IXA_ATOM_RADICAL) -> Radical {
    match code {
        inchi_sys::IXA_ATOM_RADICAL_SINGLET => Radical::Singlet,
        inchi_sys::IXA_ATOM_RADICAL_DOUBLET => Radical::Doublet,
        inchi_sys::IXA_ATOM_RADICAL_TRIPLET => Radical::Triplet,
        _ => Radical::None,
    }
}

fn decode_bond_type(code: inchi_sys::IXA_BOND_TYPE) -> BondOrder {
    match code {
        inchi_sys::IXA_BOND_TYPE_DOUBLE => BondOrder::Double,
        inchi_sys::IXA_BOND_TYPE_TRIPLE => BondOrder::Triple,
        inchi_sys::IXA_BOND_TYPE_AROMATIC => BondOrder::Alternating,
        _ => BondOrder::Single,
    }
}

fn decode_parity(code: inchi_sys::IXA_STEREO_PARITY) -> Option<Parity> {
    match code {
        inchi_sys::IXA_STEREO_PARITY_ODD => Some(Parity::Odd),
        inchi_sys::IXA_STEREO_PARITY_EVEN => Some(Parity::Even),
        inchi_sys::IXA_STEREO_PARITY_UNKNOWN => Some(Parity::Unknown),
        _ => None,
    }
}
