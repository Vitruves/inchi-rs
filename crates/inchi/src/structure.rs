//! Reverse direction: parse an InChI string back into an atom/bond/stereo
//! [`Structure`] via the native `GetStructFromINCHI` entry point.

use crate::error::{InchiError, Result, Status};
use crate::molecule::{BondOrder, Parity, Stereo};

/// One atom of a [`Structure`] recovered from an InChI.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct StructureAtom {
    /// The element symbol (e.g. `"C"`).
    pub element: String,
    /// Cartesian coordinates (zero for InChIs without coordinates).
    pub position: [f64; 3],
    /// The formal charge.
    pub charge: i8,
    /// The absolute isotopic mass, or `None` for the natural composition.
    pub isotope: Option<u16>,
    /// Total number of implicit (non-isotopic) hydrogens attached.
    pub implicit_hydrogens: u8,
}

/// One bond of a [`Structure`], referencing atoms by index.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct StructureBond {
    /// Index of the first atom.
    pub from: usize,
    /// Index of the second atom.
    pub to: usize,
    /// The bond order.
    pub order: BondOrder,
}

/// A molecular structure recovered from an InChI string.
///
/// ```
/// use inchi::struct_from_inchi;
///
/// let s = struct_from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")?;
/// assert_eq!(s.atoms().len(), 3); // two carbons and an oxygen
/// assert_eq!(s.bonds().len(), 2);
/// # Ok::<(), inchi::InchiError>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Structure {
    atoms: Vec<StructureAtom>,
    bonds: Vec<StructureBond>,
    stereo: Vec<Stereo>,
}

impl Structure {
    /// The atoms, in InChI canonical order.
    ///
    /// ```
    /// # use inchi::struct_from_inchi;
    /// let s = struct_from_inchi("InChI=1S/H2O/h1H2")?;
    /// assert_eq!(s.atoms()[0].element, "O");
    /// assert_eq!(s.atoms()[0].implicit_hydrogens, 2);
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn atoms(&self) -> &[StructureAtom] {
        &self.atoms
    }

    /// The bonds, each listed once.
    ///
    /// ```
    /// # use inchi::struct_from_inchi;
    /// let s = struct_from_inchi("InChI=1S/CH4/h1H4")?;
    /// assert!(s.bonds().is_empty()); // lone carbon, hydrogens are implicit
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn bonds(&self) -> &[StructureBond] {
        &self.bonds
    }

    /// The 0D stereo descriptors recovered from the InChI.
    ///
    /// ```
    /// # use inchi::struct_from_inchi;
    /// let s = struct_from_inchi("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1")?;
    /// assert_eq!(s.stereo().len(), 1);
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn stereo(&self) -> &[Stereo] {
        &self.stereo
    }
}

/// A [`Structure`] together with any polymer data recovered from an InChI by
/// [`struct_from_inchi_ex`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ExtendedStructure {
    /// The atom/bond/stereo structure (identical to [`struct_from_inchi`]).
    pub structure: Structure,
    /// The polymer structural repeating units, empty for non-polymer InChIs.
    pub polymer_units: Vec<crate::polymer::PolymerUnit>,
}

/// RAII guard ensuring `FreeStructFromINCHI` runs for the output struct.
struct StructGuard {
    raw: inchi_sys::inchi_OutputStruct,
}

impl StructGuard {
    fn new() -> Self {
        StructGuard {
            raw: unsafe { std::mem::zeroed() },
        }
    }
}

impl Drop for StructGuard {
    fn drop(&mut self) {
        unsafe { inchi_sys::FreeStructFromINCHI(&mut self.raw) }
    }
}

/// Parses an InChI string into its atom/bond/stereo [`Structure`].
///
/// This reverses [`from_molfile`](crate::from_molfile) / [`Molecule::to_inchi`](crate::Molecule::to_inchi)
/// (modulo information InChI does not retain, such as exact coordinates).
///
/// # Errors
///
/// Returns [`InchiError::Failed`] if the InChI is invalid or cannot be expanded,
/// and [`InchiError::InteriorNul`] if the input contains a NUL byte.
///
/// ```
/// use inchi::struct_from_inchi;
///
/// let s = struct_from_inchi("InChI=1S/CH4/h1H4")?;
/// assert_eq!(s.atoms().len(), 1);
/// assert_eq!(s.atoms()[0].element, "C");
/// assert_eq!(s.atoms()[0].implicit_hydrogens, 4);
/// # Ok::<(), inchi::InchiError>(())
/// ```
pub fn struct_from_inchi(inchi: impl AsRef<str>) -> Result<Structure> {
    let src = crate::raw::to_cstring(inchi.as_ref())?;
    // No options are required for the reverse direction.
    let empty = crate::raw::to_cstring("")?;

    let mut input = inchi_sys::inchi_InputINCHI {
        szInChI: src.as_ptr() as *mut std::os::raw::c_char,
        szOptions: empty.as_ptr() as *mut std::os::raw::c_char,
    };

    let _guard = crate::raw::lock();
    let mut out = StructGuard::new();
    let rc = unsafe { inchi_sys::GetStructFromINCHI(&mut input, &mut out.raw) };
    drop(src);
    drop(empty);

    let status = Status::from_code(rc);
    if !status.is_success() {
        let message = unsafe { crate::raw::cstr_to_string(out.raw.szMessage) };
        return Err(InchiError::Failed { status, message });
    }

    // SAFETY: on success the library populated `num_atoms` atoms at `atom` and
    // `num_stereo0D` stereo elements at `stereo0D` (either may be null/zero).
    unsafe { convert_raw(out.raw.atom, out.raw.num_atoms, out.raw.stereo0D, out.raw.num_stereo0D) }
}

/// RAII guard ensuring `FreeStructFromINCHIEx` runs for the extended output.
struct StructExGuard {
    raw: inchi_sys::inchi_OutputStructEx,
}

impl StructExGuard {
    fn new() -> Self {
        StructExGuard {
            raw: unsafe { std::mem::zeroed() },
        }
    }
}

impl Drop for StructExGuard {
    fn drop(&mut self) {
        unsafe { inchi_sys::FreeStructFromINCHIEx(&mut self.raw) }
    }
}

/// Parses an InChI into its structure *and* polymer data via the extended
/// `GetStructFromINCHIEx` entry point.
///
/// For an ordinary (non-polymer) InChI this returns the same atoms, bonds, and
/// stereo as [`struct_from_inchi`], with an empty
/// [`polymer_units`](ExtendedStructure::polymer_units). For a polymer InChI
/// (produced with [`Options::polymers`](crate::Options::polymers)) it
/// additionally recovers the structural repeating units.
///
/// # Errors
///
/// Returns [`InchiError::Failed`] if the InChI is invalid or cannot be expanded,
/// and [`InchiError::InteriorNul`] if the input contains a NUL byte.
///
/// ```
/// use inchi::struct_from_inchi_ex;
///
/// // An ordinary InChI yields no polymer units.
/// let ext = struct_from_inchi_ex("InChI=1S/CH4/h1H4")?;
/// assert_eq!(ext.structure.atoms().len(), 1);
/// assert!(ext.polymer_units.is_empty());
/// # Ok::<(), inchi::InchiError>(())
/// ```
pub fn struct_from_inchi_ex(inchi: impl AsRef<str>) -> Result<ExtendedStructure> {
    let src = crate::raw::to_cstring(inchi.as_ref())?;
    let empty = crate::raw::to_cstring("")?;

    let mut input = inchi_sys::inchi_InputINCHI {
        szInChI: src.as_ptr() as *mut std::os::raw::c_char,
        szOptions: empty.as_ptr() as *mut std::os::raw::c_char,
    };

    let _guard = crate::raw::lock();
    let mut out = StructExGuard::new();
    let rc = unsafe { inchi_sys::GetStructFromINCHIEx(&mut input, &mut out.raw) };
    drop(src);
    drop(empty);

    let status = Status::from_code(rc);
    if !status.is_success() {
        let message = unsafe { crate::raw::cstr_to_string(out.raw.szMessage) };
        return Err(InchiError::Failed { status, message });
    }

    // SAFETY: on success the atom/stereo arrays are populated as in the plain
    // struct output; the field layout matches `convert_raw`'s expectations.
    let structure = unsafe {
        convert_raw(out.raw.atom, out.raw.num_atoms, out.raw.stereo0D, out.raw.num_stereo0D)
    }?;

    // SAFETY: `polymer`, when non-null, points to a populated polymer block.
    let polymer_units = unsafe { read_polymer(out.raw.polymer) };

    Ok(ExtendedStructure { structure, polymer_units })
}

/// Reads the polymer block of an extended output into safe [`PolymerUnit`]s.
///
/// # Safety
///
/// `polymer` must be null or a valid pointer to a populated `inchi_Output_Polymer`.
unsafe fn read_polymer(
    polymer: *const inchi_sys::inchi_Input_Polymer,
) -> Vec<crate::polymer::PolymerUnit> {
    if polymer.is_null() {
        return Vec::new();
    }
    // SAFETY: caller guarantees `polymer` is valid when non-null.
    let p = unsafe { &*polymer };
    let n = usize::try_from(p.n).unwrap_or(0);
    if n == 0 || p.units.is_null() {
        return Vec::new();
    }
    // SAFETY: `units` is an array of `n` unit pointers.
    let units = unsafe { std::slice::from_raw_parts(p.units, n) };
    units
        .iter()
        .filter_map(|&u| unsafe { crate::polymer::read_unit(u) })
        .collect()
}

/// Reconstructs a [`Structure`] from the AuxInfo emitted alongside an InChI.
///
/// The auxiliary-information string ([`InchiOutput::aux_info`](crate::InchiOutput::aux_info))
/// records the original atom numbering, so it recovers a structure closer to
/// the input than [`struct_from_inchi`] does. The input may be either a full
/// InChI output block or just the `AuxInfo=` line.
///
/// When `add_hydrogens` is `true` the library is allowed to add implicit
/// hydrogens to complete normal valences (the usual choice); pass `false` to
/// keep only the hydrogens explicitly recorded.
///
/// # Errors
///
/// Returns [`InchiError::Failed`] if the AuxInfo cannot be parsed,
/// [`InchiError::EmptyResult`] if it yields no atoms, and
/// [`InchiError::InteriorNul`] if the input contains a NUL byte.
///
/// ```
/// use inchi::{from_molfile, struct_from_aux_info};
///
/// let methane = "\n  ex\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n\
///     \x20   0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
/// let out = from_molfile(methane, ())?;
/// let s = struct_from_aux_info(out.aux_info(), true)?;
/// assert_eq!(s.atoms()[0].element, "C");
/// # Ok::<(), inchi::InchiError>(())
/// ```
pub fn struct_from_aux_info(aux_info: impl AsRef<str>, add_hydrogens: bool) -> Result<Structure> {
    let aux = crate::raw::to_cstring(aux_info.as_ref())?;

    let _guard = crate::raw::lock();
    let mut guard = InputGuard::new();
    // `bDoNotAddH` is the inverse of `add_hydrogens`; `bDiffUnkUndfStereo = 0`
    // keeps the default of merging unknown and undefined stereo labels.
    let b_do_not_add_h = i32::from(!add_hydrogens);
    // SAFETY: `aux` is a valid NUL-terminated string for the call; the library
    // treats it as read-only despite the `*mut` signature. `guard.data.pInp`
    // points at the zeroed `inchi_Input` that `guard` owns and frees on drop.
    let rc = unsafe {
        inchi_sys::Get_inchi_Input_FromAuxInfo(
            aux.as_ptr() as *mut std::os::raw::c_char,
            b_do_not_add_h,
            0,
            &mut guard.data,
        )
    };
    drop(aux);

    let status = Status::from_code(rc);
    if !status.is_success() {
        let message = unsafe { read_err_msg(&guard.data.szErrMsg) };
        return Err(InchiError::Failed { status, message });
    }

    // SAFETY: on success `pInp` points to a populated `inchi_Input`.
    let inp = &guard.input;
    let s = unsafe { convert_raw(inp.atom, inp.num_atoms, inp.stereo0D, inp.num_stereo0D) }?;
    if s.atoms.is_empty() {
        return Err(InchiError::EmptyResult);
    }
    Ok(s)
}

/// RAII guard owning the `inchi_Input` that `Get_inchi_Input_FromAuxInfo`
/// populates, ensuring `Free_inchi_Input` runs for whatever it allocated.
struct InputGuard {
    input: Box<inchi_sys::inchi_Input>,
    data: inchi_sys::InchiInpData,
}

impl InputGuard {
    fn new() -> Self {
        // `input` is boxed so its address is stable while `data.pInp` borrows it.
        let mut input = Box::new(unsafe { std::mem::zeroed::<inchi_sys::inchi_Input>() });
        let mut data: inchi_sys::InchiInpData = unsafe { std::mem::zeroed() };
        data.pInp = input.as_mut();
        InputGuard { input, data }
    }
}

impl Drop for InputGuard {
    fn drop(&mut self) {
        // Frees the atom/stereo arrays the library allocated and zeroes the
        // changed members; safe to call even when nothing was populated.
        unsafe { inchi_sys::Free_inchi_Input(self.input.as_mut()) }
    }
}

/// Reads the fixed-size `szErrMsg` buffer into an owned `String`.
///
/// # Safety
///
/// `buf` must be a valid array holding a NUL-terminated C string.
unsafe fn read_err_msg(buf: &[std::os::raw::c_char]) -> String {
    crate::raw::cstr_to_string(buf.as_ptr())
}

/// Converts a populated atom/stereo array (shared by `inchi_Input` and
/// `inchi_OutputStruct`) into a safe [`Structure`].
///
/// # Safety
///
/// When `num_atoms > 0`, `atom` must point to that many valid [`inchi_Atom`](inchi_sys::inchi_Atom)
/// values; likewise `stereo0D` must hold `num_stereo0D` valid elements.
unsafe fn convert_raw(
    atom: *mut inchi_sys::inchi_Atom,
    num_atoms: inchi_sys::AT_NUM,
    stereo0d: *mut inchi_sys::inchi_Stereo0D,
    num_stereo0d: inchi_sys::AT_NUM,
) -> Result<Structure> {
    let num_atoms = usize::try_from(num_atoms).unwrap_or(0);
    if num_atoms == 0 || atom.is_null() {
        return Ok(Structure {
            atoms: Vec::new(),
            bonds: Vec::new(),
            stereo: Vec::new(),
        });
    }

    // SAFETY: the caller guarantees `num_atoms` valid atoms at `atom`.
    let c_atoms = unsafe { std::slice::from_raw_parts(atom, num_atoms) };

    let mut atoms = Vec::with_capacity(num_atoms);
    let mut bonds = Vec::new();
    for (i, ca) in c_atoms.iter().enumerate() {
        atoms.push(StructureAtom {
            element: read_elname(&ca.elname),
            position: [ca.x, ca.y, ca.z],
            charge: ca.charge,
            isotope: read_isotope(ca.isotopic_mass),
            implicit_hydrogens: read_implicit_h(&ca.num_iso_H),
        });

        // Each bond appears in both atoms' adjacency lists; keep it once.
        let degree = (ca.num_bonds.max(0) as usize).min(ca.neighbor.len());
        for slot in 0..degree {
            let (Some(&nbr), Some(&bt)) = (ca.neighbor.get(slot), ca.bond_type.get(slot)) else {
                continue;
            };
            let j = usize::try_from(nbr).unwrap_or(usize::MAX);
            if j < num_atoms && i < j {
                bonds.push(StructureBond {
                    from: i,
                    to: j,
                    order: decode_order(bt),
                });
            }
        }
    }

    // SAFETY: forwarded from the caller's guarantee on `stereo0d`/`num_stereo0d`.
    let stereo = unsafe { read_stereo(stereo0d, num_stereo0d, num_atoms) };

    Ok(Structure { atoms, bonds, stereo })
}

/// # Safety
///
/// When `num_stereo0d > 0`, `stereo0d` must point to that many valid elements.
unsafe fn read_stereo(
    stereo0d: *mut inchi_sys::inchi_Stereo0D,
    num_stereo0d: inchi_sys::AT_NUM,
    num_atoms: usize,
) -> Vec<Stereo> {
    let count = usize::try_from(num_stereo0d).unwrap_or(0);
    if count == 0 || stereo0d.is_null() {
        return Vec::new();
    }
    // SAFETY: the caller guarantees `count` valid stereo elements.
    let c_stereo = unsafe { std::slice::from_raw_parts(stereo0d, count) };

    let idx = |v: inchi_sys::AT_NUM| usize::try_from(v).ok().filter(|&i| i < num_atoms);
    let mut out = Vec::with_capacity(count);
    for cs in c_stereo {
        let Some(parity) = decode_parity(cs.parity) else {
            continue;
        };
        let Some(ends) = (|| {
            Some([
                idx(cs.neighbor[0])?,
                idx(cs.neighbor[1])?,
                idx(cs.neighbor[2])?,
                idx(cs.neighbor[3])?,
            ])
        })() else {
            continue;
        };
        let ty = cs.type_ as u32;
        if ty == inchi_sys::INCHI_StereoType_Tetrahedral {
            if let Some(center) = idx(cs.central_atom) {
                out.push(Stereo::Tetrahedral { center, neighbors: ends, parity });
            }
        } else if ty == inchi_sys::INCHI_StereoType_DoubleBond {
            out.push(Stereo::DoubleBond { ends, parity });
        } else if ty == inchi_sys::INCHI_StereoType_Allene {
            if let Some(center) = idx(cs.central_atom) {
                out.push(Stereo::Allene { center, ends, parity });
            }
        }
    }
    out
}

fn read_elname(raw: &[std::os::raw::c_char]) -> String {
    let mut s = String::new();
    for &c in raw {
        if c == 0 {
            break;
        }
        s.push(c as u8 as char);
    }
    s
}

fn read_isotope(mass: inchi_sys::AT_NUM) -> Option<u16> {
    if mass == 0 {
        None
    } else {
        u16::try_from(mass).ok()
    }
}

fn read_implicit_h(num_iso_h: &[i8]) -> u8 {
    // num_iso_H[0] is the count of non-isotopic implicit H (or -1 for "auto",
    // which never appears in library output). Sum all isotopes for the total.
    num_iso_h
        .iter()
        .map(|&n| if n < 0 { 0 } else { n as u16 })
        .sum::<u16>()
        .min(u8::MAX as u16) as u8
}

fn decode_order(code: i8) -> BondOrder {
    let c = code as u32;
    if c == inchi_sys::INCHI_BOND_TYPE_DOUBLE {
        BondOrder::Double
    } else if c == inchi_sys::INCHI_BOND_TYPE_TRIPLE {
        BondOrder::Triple
    } else if c == inchi_sys::INCHI_BOND_TYPE_ALTERN {
        BondOrder::Alternating
    } else {
        BondOrder::Single
    }
}

fn decode_parity(code: i8) -> Option<Parity> {
    // The low 3 bits hold the parity of the connected structure.
    match (code & 0x07) as u32 {
        inchi_sys::INCHI_PARITY_ODD => Some(Parity::Odd),
        inchi_sys::INCHI_PARITY_EVEN => Some(Parity::Even),
        inchi_sys::INCHI_PARITY_UNKNOWN => Some(Parity::Unknown),
        _ => None,
    }
}
