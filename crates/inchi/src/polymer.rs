//! Polymer (CTFile Sgroup) data model shared by the extended generation and
//! parsing entry points (`GetINCHIEx` / `GetStructFromINCHIEx`).
//!
//! Polymer support is an **experimental, non-standard** InChI extension: the
//! generated identifier carries a `B` ("beta") version flag and a `/z` polymer
//! layer. A polymer is described by one or more structural repeating units
//! (SRUs), each spanning a set of atoms and crossed by the bonds that connect it
//! to the rest of the structure.

/// The CTFile `STY` unit type of a [`PolymerUnit`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[non_exhaustive]
pub enum PolymerUnitKind {
    /// No type specified (`POLYMER_STY_NON`).
    #[default]
    None,
    /// Structure-based repeating unit or copolymer subunit (`SRU`).
    StructureBasedRepeat,
    /// Source-based polymer unit or copolymer subunit (`MON`).
    SourceBasedMonomer,
    /// Copolymer unit embedding more than one subunit (`COP`).
    Copolymer,
    /// Chemical modification of an SRU (`MOD`).
    Modification,
    /// Cross-linked version of an SRU (`CRO`).
    CrossLink,
    /// Source-based subunit with no homopolymerization (`MER`).
    Mer,
    /// A code outside the documented set, preserved verbatim.
    Other(i32),
}

impl PolymerUnitKind {
    fn from_code(code: i32) -> Self {
        match code {
            0 => PolymerUnitKind::None,
            1 => PolymerUnitKind::StructureBasedRepeat,
            2 => PolymerUnitKind::SourceBasedMonomer,
            3 => PolymerUnitKind::Copolymer,
            4 => PolymerUnitKind::Modification,
            5 => PolymerUnitKind::CrossLink,
            6 => PolymerUnitKind::Mer,
            other => PolymerUnitKind::Other(other),
        }
    }

    pub(crate) fn code(self) -> i32 {
        match self {
            PolymerUnitKind::None => 0,
            PolymerUnitKind::StructureBasedRepeat => 1,
            PolymerUnitKind::SourceBasedMonomer => 2,
            PolymerUnitKind::Copolymer => 3,
            PolymerUnitKind::Modification => 4,
            PolymerUnitKind::CrossLink => 5,
            PolymerUnitKind::Mer => 6,
            PolymerUnitKind::Other(c) => c,
        }
    }
}

/// The CTFile `SST` unit subtype of a [`PolymerUnit`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[non_exhaustive]
pub enum PolymerSubtype {
    /// No subtype specified (`POLYMER_SST_NON`).
    #[default]
    None,
    /// Alternating (`ALT`).
    Alternating,
    /// Random (`RAN`).
    Random,
    /// Block (`BLK`).
    Block,
    /// A code outside the documented set, preserved verbatim.
    Other(i32),
}

impl PolymerSubtype {
    fn from_code(code: i32) -> Self {
        match code {
            0 => PolymerSubtype::None,
            1 => PolymerSubtype::Alternating,
            2 => PolymerSubtype::Random,
            3 => PolymerSubtype::Block,
            other => PolymerSubtype::Other(other),
        }
    }

    pub(crate) fn code(self) -> i32 {
        match self {
            PolymerSubtype::None => 0,
            PolymerSubtype::Alternating => 1,
            PolymerSubtype::Random => 2,
            PolymerSubtype::Block => 3,
            PolymerSubtype::Other(c) => c,
        }
    }
}

/// The CTFile `SCN` connection scheme of a [`PolymerUnit`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[non_exhaustive]
pub enum PolymerConnection {
    /// Unspecified connection.
    #[default]
    Unspecified,
    /// Head-to-tail (`HT`).
    HeadToTail,
    /// Head-to-head (`HH`).
    HeadToHead,
    /// Either/unknown (`EU`).
    EitherUnknown,
    /// A code outside the documented set, preserved verbatim.
    Other(i32),
}

impl PolymerConnection {
    fn from_code(code: i32) -> Self {
        match code {
            0 => PolymerConnection::Unspecified,
            1 => PolymerConnection::HeadToTail,
            2 => PolymerConnection::HeadToHead,
            3 => PolymerConnection::EitherUnknown,
            other => PolymerConnection::Other(other),
        }
    }

    pub(crate) fn code(self) -> i32 {
        match self {
            PolymerConnection::Unspecified => 0,
            PolymerConnection::HeadToTail => 1,
            PolymerConnection::HeadToHead => 2,
            PolymerConnection::EitherUnknown => 3,
            PolymerConnection::Other(c) => c,
        }
    }
}

/// A single polymer structural repeating unit (CTFile Sgroup).
///
/// Atom and bond indices are **0-based** and refer to the surrounding
/// [`Molecule`](crate::Molecule) / [`Structure`](crate::Structure); the
/// conversion to and from the library's 1-based numbering is handled internally.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct PolymerUnit {
    /// The Sgroup number (`id`); informational, preserved for round-tripping.
    pub id: i32,
    /// The unique Sgroup identifier (`label`); informational.
    pub label: i32,
    /// The unit type (CTFile `STY`).
    pub kind: PolymerUnitKind,
    /// The unit subtype (CTFile `SST`).
    pub subtype: PolymerSubtype,
    /// The connection scheme (CTFile `SCN`).
    pub connection: PolymerConnection,
    /// The Sgroup subscript (CTFile `SMT`), e.g. `"n"`.
    pub subscript: String,
    /// Indices of the atoms inside the unit (CTFile `SAL`), 0-based.
    pub atoms: Vec<usize>,
    /// The unit's crossing bonds (CTFile `SBL`), each as a `[from, to]` pair of
    /// 0-based atom indices.
    pub crossing_bonds: Vec<[usize; 2]>,
}

impl PolymerUnit {
    /// Creates a structure-based repeating unit (`SRU`) over the given atoms.
    ///
    /// This is the common case; refine the result with the public fields for
    /// other unit types, subtypes, or connection schemes.
    ///
    /// ```
    /// use inchi::PolymerUnit;
    /// let unit = PolymerUnit::sru([1, 2], [[0, 1], [2, 3]]);
    /// assert_eq!(unit.atoms, vec![1, 2]);
    /// ```
    #[must_use]
    pub fn sru(
        atoms: impl IntoIterator<Item = usize>,
        crossing_bonds: impl IntoIterator<Item = [usize; 2]>,
    ) -> Self {
        PolymerUnit {
            id: 1,
            label: 0,
            kind: PolymerUnitKind::StructureBasedRepeat,
            subtype: PolymerSubtype::None,
            connection: PolymerConnection::HeadToTail,
            subscript: "n".to_string(),
            atoms: atoms.into_iter().collect(),
            crossing_bonds: crossing_bonds.into_iter().collect(),
        }
    }
}

/// Reads a C polymer unit (1-based indices) into a safe [`PolymerUnit`].
///
/// # Safety
///
/// `raw` must be a valid pointer to an `inchi_Input_PolymerUnit` whose `alist`
/// and `blist` arrays hold `na` and `2 * nb` valid entries respectively.
pub(crate) unsafe fn read_unit(raw: *const inchi_sys::inchi_Input_PolymerUnit) -> Option<PolymerUnit> {
    if raw.is_null() {
        return None;
    }
    // SAFETY: caller guarantees `raw` is a valid unit pointer.
    let u = unsafe { &*raw };

    let na = usize::try_from(u.na).unwrap_or(0);
    let nb = usize::try_from(u.nb).unwrap_or(0);

    let atoms = if u.alist.is_null() || na == 0 {
        Vec::new()
    } else {
        // SAFETY: `alist` holds `na` 1-based atom numbers.
        let slice = unsafe { std::slice::from_raw_parts(u.alist, na) };
        slice.iter().map(|&a| zero_based(a)).collect()
    };

    let crossing_bonds = if u.blist.is_null() || nb == 0 {
        Vec::new()
    } else {
        // SAFETY: `blist` holds `2 * nb` 1-based atom numbers (pairs).
        let slice = unsafe { std::slice::from_raw_parts(u.blist, nb * 2) };
        slice
            .chunks_exact(2)
            .filter_map(|p| match p {
                [a, b] => Some([zero_based(*a), zero_based(*b)]),
                _ => None,
            })
            .collect()
    };

    Some(PolymerUnit {
        id: u.id,
        label: u.label,
        kind: PolymerUnitKind::from_code(u.type_),
        subtype: PolymerSubtype::from_code(u.subtype),
        connection: PolymerConnection::from_code(u.conn),
        subscript: read_smt(&u.smt),
        atoms,
        crossing_bonds,
    })
}

/// Converts a 1-based C atom number to a 0-based index (clamping the unused
/// `0`/negative sentinels to index 0).
fn zero_based(n: std::os::raw::c_int) -> usize {
    (n.max(1) as usize) - 1
}

/// Reads the fixed-size `smt` subscript field into an owned `String`.
fn read_smt(raw: &[std::os::raw::c_char]) -> String {
    let mut s = String::new();
    for &c in raw {
        if c == 0 {
            break;
        }
        s.push(c as u8 as char);
    }
    s
}
