//! Programmatic structure input: build a [`Molecule`] from atoms, bonds, and
//! 0D stereo, then generate its InChI via the native `GetINCHI` entry point.

use crate::error::{InchiError, Result};
use crate::options::Options;
use crate::output::InchiOutput;

/// Maximum number of atoms accepted by the InChI library (`MAX_ATOMS`).
const MAX_ATOMS: usize = 1024;
/// Maximum bonds recordable per atom in the FFI struct (`MAXVAL`).
const MAX_BONDS_PER_ATOM: usize = inchi_sys::MAXVAL as usize;
/// Capacity of the element-symbol field, including the trailing NUL.
const ELNAME_CAP: usize = inchi_sys::ATOM_EL_LEN as usize;

/// The unpaired-electron (radical) state of an atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum Radical {
    /// No radical.
    #[default]
    None,
    /// Singlet (carbene-like).
    Singlet,
    /// Doublet (one unpaired electron).
    Doublet,
    /// Triplet (two unpaired electrons).
    Triplet,
}

impl Radical {
    fn code(self) -> i8 {
        let v = match self {
            Radical::None => inchi_sys::INCHI_RADICAL_NONE,
            Radical::Singlet => inchi_sys::INCHI_RADICAL_SINGLET,
            Radical::Doublet => inchi_sys::INCHI_RADICAL_DOUBLET,
            Radical::Triplet => inchi_sys::INCHI_RADICAL_TRIPLET,
        };
        v as i8
    }
}

/// The order of a covalent bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum BondOrder {
    /// A single bond.
    #[default]
    Single,
    /// A double bond.
    Double,
    /// A triple bond.
    Triple,
    /// An "alternating"/aromatic bond. The InChI documentation recommends
    /// avoiding this in favor of explicit single/double bonds.
    Alternating,
}

impl BondOrder {
    fn code(self) -> i8 {
        let v = match self {
            BondOrder::Single => inchi_sys::INCHI_BOND_TYPE_SINGLE,
            BondOrder::Double => inchi_sys::INCHI_BOND_TYPE_DOUBLE,
            BondOrder::Triple => inchi_sys::INCHI_BOND_TYPE_TRIPLE,
            BondOrder::Alternating => inchi_sys::INCHI_BOND_TYPE_ALTERN,
        };
        v as i8
    }
}

/// How many implicit hydrogens an atom carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum ImplicitH {
    /// Let the library add implicit hydrogens to satisfy normal valence
    /// (`num_iso_H[0] = -1`). This is the usual choice for a heavy-atom
    /// skeleton and mirrors molfile behavior.
    #[default]
    Auto,
    /// Attach exactly this many implicit (non-isotopic) hydrogens.
    Exactly(u8),
}

/// A 0D stereo parity (used when no coordinates disambiguate the geometry).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Parity {
    /// Odd parity (`'-'` / `'o'`).
    Odd,
    /// Even parity (`'+'` / `'e'`).
    Even,
    /// Known to be stereogenic but of unspecified configuration (`'u'`).
    Unknown,
}

impl Parity {
    fn code(self) -> i8 {
        let v = match self {
            Parity::Odd => inchi_sys::INCHI_PARITY_ODD,
            Parity::Even => inchi_sys::INCHI_PARITY_EVEN,
            Parity::Unknown => inchi_sys::INCHI_PARITY_UNKNOWN,
        };
        v as i8
    }
}

/// A single 0D stereo element, referencing atoms by their index in the
/// [`Molecule`].
///
/// The neighbor ordering follows the InChI convention exactly; getting it wrong
/// flips the parity, so validate against known identifiers. See the upstream
/// `inchi_api.h` for the precise diagrams.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Stereo {
    /// A tetrahedral stereocenter. `neighbors` lists the four substituents in
    /// the order whose handedness defines `parity`, seen from the first
    /// neighbor toward `center`.
    Tetrahedral {
        /// Index of the central atom.
        center: usize,
        /// The four neighbor atom indices, in convention order.
        neighbors: [usize; 4],
        /// The parity of the configuration.
        parity: Parity,
    },
    /// A stereogenic double bond `>A=B<` (or even-length cumulene). `ends` is
    /// `[X, A, B, Y]` where `A=B` is the double bond and `X`, `Y` are the
    /// reference substituents.
    DoubleBond {
        /// `[X, A, B, Y]` atom indices.
        ends: [usize; 4],
        /// The parity of the configuration.
        parity: Parity,
    },
    /// An allene / odd-length cumulene stereocenter. `ends` is `[X, A, B, Y]`
    /// and `center` is the central cumulene atom.
    Allene {
        /// Index of the central atom.
        center: usize,
        /// `[X, A, B, Y]` atom indices.
        ends: [usize; 4],
        /// The parity of the configuration.
        parity: Parity,
    },
}

/// A single atom in a [`Molecule`].
///
/// Construct with [`Atom::new`] and refine with the chainable setters.
///
/// ```
/// use inchi::{Atom, Radical, ImplicitH};
///
/// let carbon = Atom::new("C").position(0.0, 0.0, 0.0);
/// let chloride = Atom::new("Cl").charge(-1).implicit_hydrogens(ImplicitH::Exactly(0));
/// let _ = (carbon, chloride, Radical::None);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    element: String,
    x: f64,
    y: f64,
    z: f64,
    charge: i8,
    isotope: Option<u16>,
    radical: Radical,
    implicit_h: ImplicitH,
}

impl Atom {
    /// Creates an atom of the given element (e.g. `"C"`, `"Cl"`, `"Na"`).
    ///
    /// The symbol is validated when the molecule is converted; an empty,
    /// non-ASCII, or over-long symbol yields [`InchiError::InvalidStructure`].
    ///
    /// ```
    /// use inchi::Atom;
    /// let _ = Atom::new("O");
    /// ```
    #[must_use]
    pub fn new(element: impl Into<String>) -> Self {
        Atom {
            element: element.into(),
            x: 0.0,
            y: 0.0,
            z: 0.0,
            charge: 0,
            isotope: None,
            radical: Radical::None,
            implicit_h: ImplicitH::Auto,
        }
    }

    /// Sets the 3D coordinates of the atom (defaults to the origin).
    ///
    /// ```
    /// use inchi::Atom;
    /// let _ = Atom::new("C").position(1.0, 0.5, -0.25);
    /// ```
    #[must_use]
    pub fn position(mut self, x: f64, y: f64, z: f64) -> Self {
        self.x = x;
        self.y = y;
        self.z = z;
        self
    }

    /// Sets the formal charge (defaults to `0`).
    ///
    /// ```
    /// use inchi::Atom;
    /// let _ = Atom::new("N").charge(1);
    /// ```
    #[must_use]
    pub fn charge(mut self, charge: i8) -> Self {
        self.charge = charge;
        self
    }

    /// Sets the absolute isotopic mass (e.g. `13` for carbon-13). Omit for the
    /// natural isotopic composition.
    ///
    /// ```
    /// use inchi::Atom;
    /// let _ = Atom::new("C").isotope(13);
    /// ```
    #[must_use]
    pub fn isotope(mut self, mass: u16) -> Self {
        self.isotope = Some(mass);
        self
    }

    /// Sets the radical state (defaults to [`Radical::None`]).
    ///
    /// ```
    /// use inchi::{Atom, Radical};
    /// let _ = Atom::new("C").radical(Radical::Triplet);
    /// ```
    #[must_use]
    pub fn radical(mut self, radical: Radical) -> Self {
        self.radical = radical;
        self
    }

    /// Sets how implicit hydrogens are handled (defaults to [`ImplicitH::Auto`]).
    ///
    /// ```
    /// use inchi::{Atom, ImplicitH};
    /// let _ = Atom::new("C").implicit_hydrogens(ImplicitH::Exactly(3));
    /// ```
    #[must_use]
    pub fn implicit_hydrogens(mut self, h: ImplicitH) -> Self {
        self.implicit_h = h;
        self
    }
}

/// A molecular structure assembled programmatically from atoms, bonds, and 0D
/// stereo descriptors.
///
/// Atoms are referenced by the index returned from [`Molecule::add_atom`] (also
/// the order in which they are added, starting at `0`).
///
/// ```
/// use inchi::{Molecule, Atom, BondOrder};
///
/// // Ethanol: C-C-O (implicit hydrogens added automatically).
/// let mut mol = Molecule::new();
/// let c1 = mol.add_atom(Atom::new("C"));
/// let c2 = mol.add_atom(Atom::new("C"));
/// let o = mol.add_atom(Atom::new("O"));
/// mol.add_bond(c1, c2, BondOrder::Single)?;
/// mol.add_bond(c2, o, BondOrder::Single)?;
///
/// let out = mol.to_inchi(())?;
/// assert_eq!(out.inchi(), "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3");
/// # Ok::<(), inchi::InchiError>(())
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<(usize, usize, BondOrder)>,
    stereo: Vec<Stereo>,
    polymer_units: Vec<crate::polymer::PolymerUnit>,
}

impl Molecule {
    /// Creates an empty molecule.
    ///
    /// ```
    /// use inchi::Molecule;
    /// let mol = Molecule::new();
    /// assert_eq!(mol.atom_count(), 0);
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Molecule::default()
    }

    /// Adds an atom and returns its index.
    ///
    /// ```
    /// use inchi::{Molecule, Atom};
    /// let mut mol = Molecule::new();
    /// assert_eq!(mol.add_atom(Atom::new("C")), 0);
    /// assert_eq!(mol.add_atom(Atom::new("O")), 1);
    /// ```
    pub fn add_atom(&mut self, atom: Atom) -> usize {
        self.atoms.push(atom);
        self.atoms.len() - 1
    }

    /// Adds a bond between two existing atoms.
    ///
    /// Returns [`InchiError::InvalidStructure`] if either index is out of range
    /// or if `a == b`.
    ///
    /// ```
    /// use inchi::{Molecule, Atom, BondOrder};
    /// let mut mol = Molecule::new();
    /// let a = mol.add_atom(Atom::new("C"));
    /// let b = mol.add_atom(Atom::new("O"));
    /// mol.add_bond(a, b, BondOrder::Double)?;
    /// assert!(mol.add_bond(a, 99, BondOrder::Single).is_err());
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    pub fn add_bond(&mut self, a: usize, b: usize, order: BondOrder) -> Result<()> {
        if a == b {
            return Err(InchiError::InvalidStructure {
                reason: format!("bond connects atom {a} to itself"),
            });
        }
        let n = self.atoms.len();
        if a >= n || b >= n {
            return Err(InchiError::InvalidStructure {
                reason: format!("bond ({a}, {b}) references a nonexistent atom (have {n})"),
            });
        }
        self.bonds.push((a, b, order));
        Ok(())
    }

    /// Adds a 0D stereo descriptor. Indices are validated at conversion time.
    ///
    /// ```
    /// use inchi::{Molecule, Atom, Stereo, Parity};
    /// let mut mol = Molecule::new();
    /// for el in ["C", "C", "N", "O"] { mol.add_atom(Atom::new(el)); }
    /// mol.add_stereo(Stereo::Tetrahedral { center: 0, neighbors: [1, 2, 3, 0], parity: Parity::Odd });
    /// assert_eq!(mol.stereo_count(), 1);
    /// ```
    pub fn add_stereo(&mut self, stereo: Stereo) {
        self.stereo.push(stereo);
    }

    /// Adds a polymer structural repeating unit, switching InChI generation to
    /// the extended `GetINCHIEx` entry point.
    ///
    /// Polymers require the [`Options::polymers`](crate::Options::polymers)
    /// flag to be set and yield a non-standard, beta-flagged InChI. The unit's
    /// atom indices refer to atoms already added to this molecule; the two
    /// chain ends are normally capped with `"Zz"` star atoms.
    ///
    /// ```
    /// use inchi::{Molecule, Atom, BondOrder, Options, Polymers, PolymerUnit};
    /// // A polyethylene repeat unit: *-CH2-CH2-*
    /// let mut mol = Molecule::new();
    /// let s1 = mol.add_atom(Atom::new("Zz"));
    /// let c1 = mol.add_atom(Atom::new("C"));
    /// let c2 = mol.add_atom(Atom::new("C"));
    /// let s2 = mol.add_atom(Atom::new("Zz"));
    /// mol.add_bond(s1, c1, BondOrder::Single)?;
    /// mol.add_bond(c1, c2, BondOrder::Single)?;
    /// mol.add_bond(c2, s2, BondOrder::Single)?;
    /// mol.add_polymer_unit(PolymerUnit::sru([c1, c2], [[s1, c1], [c2, s2]]));
    /// let inchi = mol.to_inchi(Options::new().polymers(Polymers::On))?.into_inchi();
    /// assert!(inchi.contains("/z"));
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    pub fn add_polymer_unit(&mut self, unit: crate::polymer::PolymerUnit) {
        self.polymer_units.push(unit);
    }

    /// The number of atoms.
    ///
    /// ```
    /// # use inchi::{Molecule, Atom};
    /// let mut mol = Molecule::new();
    /// mol.add_atom(Atom::new("C"));
    /// assert_eq!(mol.atom_count(), 1);
    /// ```
    #[must_use]
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// The number of bonds.
    ///
    /// ```
    /// # use inchi::{Molecule, Atom, BondOrder};
    /// let mut mol = Molecule::new();
    /// let a = mol.add_atom(Atom::new("C"));
    /// let b = mol.add_atom(Atom::new("C"));
    /// mol.add_bond(a, b, BondOrder::Single)?;
    /// assert_eq!(mol.bond_count(), 1);
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    #[must_use]
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// The number of 0D stereo descriptors.
    ///
    /// ```
    /// # use inchi::Molecule;
    /// assert_eq!(Molecule::new().stereo_count(), 0);
    /// ```
    #[must_use]
    pub fn stereo_count(&self) -> usize {
        self.stereo.len()
    }

    /// Generates the InChI for this molecule using the given [`Options`].
    ///
    /// ```
    /// use inchi::{Molecule, Atom};
    /// // A lone oxygen atom becomes water once implicit H are added.
    /// let mut mol = Molecule::new();
    /// mol.add_atom(Atom::new("O"));
    /// assert_eq!(mol.to_inchi(())?.inchi(), "InChI=1S/H2O/h1H2");
    /// # Ok::<(), inchi::InchiError>(())
    /// ```
    pub fn to_inchi(&self, options: impl Into<Options>) -> Result<InchiOutput> {
        let options = options.into();
        let mut atoms = self.build_atoms()?;
        let mut stereo = self.build_stereo()?;
        let opts = crate::raw::to_cstring(&options.to_arg_string())?;

        let num_atoms = i16::try_from(atoms.len()).map_err(|_| InchiError::InvalidStructure {
            reason: format!("too many atoms ({})", atoms.len()),
        })?;
        let num_stereo = i16::try_from(stereo.len()).map_err(|_| InchiError::InvalidStructure {
            reason: format!("too many stereo elements ({})", stereo.len()),
        })?;

        // SAFETY: `input` borrows the `atoms`/`stereo`/`opts` allocations, all
        // of which outlive the `GetINCHI` call below. `GetINCHI` does not take
        // ownership of the input (the caller owns it), and we serialize via the
        // global lock. The output is owned by an `OutputGuard` that frees it.
        let mut input: inchi_sys::inchi_Input = unsafe { std::mem::zeroed() };
        input.atom = atoms.as_mut_ptr();
        input.num_atoms = num_atoms;
        input.stereo0D = if stereo.is_empty() {
            std::ptr::null_mut()
        } else {
            stereo.as_mut_ptr()
        };
        input.num_stereo0D = num_stereo;
        input.szOptions = opts.as_ptr() as *mut std::os::raw::c_char;

        if self.polymer_units.is_empty() {
            let _guard = crate::raw::lock();
            let mut out = crate::raw::OutputGuard::new();
            let rc = unsafe { inchi_sys::GetINCHI(&mut input, out.as_mut_ptr()) };
            // Keep the input-backing allocations alive until after the FFI call.
            drop(atoms);
            drop(stereo);
            drop(opts);
            return crate::build_output(rc, &out);
        }

        // Polymer input requires the extended `GetINCHIEx` entry point. Build
        // the polymer block, keeping every backing allocation alive across the
        // call. `inchi_InputEx` shares its first fields with `inchi_Input`.
        let mut backing = PolymerBacking::build(&self.polymer_units, self.atoms.len())?;
        let mut input_ex: inchi_sys::inchi_InputEx = unsafe { std::mem::zeroed() };
        input_ex.atom = atoms.as_mut_ptr();
        input_ex.num_atoms = num_atoms;
        input_ex.stereo0D = if stereo.is_empty() {
            std::ptr::null_mut()
        } else {
            stereo.as_mut_ptr()
        };
        input_ex.num_stereo0D = num_stereo;
        input_ex.szOptions = opts.as_ptr() as *mut std::os::raw::c_char;
        input_ex.polymer = backing.as_mut_ptr();

        let _guard = crate::raw::lock();
        let mut out = crate::raw::OutputGuard::new();
        let rc = unsafe { inchi_sys::GetINCHIEx(&mut input_ex, out.as_mut_ptr()) };
        // Keep all input-backing allocations alive until after the FFI call.
        drop(atoms);
        drop(stereo);
        drop(opts);
        drop(backing);
        crate::build_output(rc, &out)
    }

    fn build_atoms(&self) -> Result<Vec<inchi_sys::inchi_Atom>> {
        if self.atoms.is_empty() {
            return Err(InchiError::InvalidStructure {
                reason: "molecule has no atoms".to_string(),
            });
        }
        if self.atoms.len() > MAX_ATOMS {
            return Err(InchiError::InvalidStructure {
                reason: format!("too many atoms ({} > {MAX_ATOMS})", self.atoms.len()),
            });
        }

        let mut raw: Vec<inchi_sys::inchi_Atom> = Vec::with_capacity(self.atoms.len());
        for atom in &self.atoms {
            let mut a: inchi_sys::inchi_Atom = unsafe { std::mem::zeroed() };
            a.x = atom.x;
            a.y = atom.y;
            a.z = atom.z;
            write_elname(&mut a.elname, &atom.element)?;
            a.charge = atom.charge;
            a.radical = atom.radical.code();
            if let Some(mass) = atom.isotope {
                a.isotopic_mass =
                    i16::try_from(mass).map_err(|_| InchiError::InvalidStructure {
                        reason: format!("isotopic mass {mass} out of range"),
                    })?;
            }
            a.num_iso_H = match atom.implicit_h {
                ImplicitH::Auto => [-1, 0, 0, 0],
                ImplicitH::Exactly(n) => [
                    i8::try_from(n).map_err(|_| InchiError::InvalidStructure {
                        reason: format!("implicit H count {n} out of range"),
                    })?,
                    0,
                    0,
                    0,
                ],
            };
            raw.push(a);
        }

        // Populate the symmetric adjacency lists from the bond list. Each bond
        // is recorded in both endpoints, as a molfile-derived input would be.
        for &(a, b, order) in &self.bonds {
            push_neighbor(&mut raw, a, b, order)?;
            push_neighbor(&mut raw, b, a, order)?;
        }

        Ok(raw)
    }

    fn build_stereo(&self) -> Result<Vec<inchi_sys::inchi_Stereo0D>> {
        let n = self.atoms.len();
        let check = |idx: usize| -> Result<i16> {
            if idx >= n {
                return Err(InchiError::InvalidStructure {
                    reason: format!("stereo references nonexistent atom {idx} (have {n})"),
                });
            }
            i16::try_from(idx).map_err(|_| InchiError::InvalidStructure {
                reason: format!("atom index {idx} out of range"),
            })
        };

        let mut raw = Vec::with_capacity(self.stereo.len());
        for stereo in &self.stereo {
            let mut s: inchi_sys::inchi_Stereo0D = unsafe { std::mem::zeroed() };
            match *stereo {
                Stereo::Tetrahedral {
                    center,
                    neighbors,
                    parity,
                } => {
                    s.central_atom = check(center)?;
                    s.neighbor = [
                        check(neighbors[0])?,
                        check(neighbors[1])?,
                        check(neighbors[2])?,
                        check(neighbors[3])?,
                    ];
                    s.type_ = inchi_sys::INCHI_StereoType_Tetrahedral as i8;
                    s.parity = parity.code();
                }
                Stereo::DoubleBond { ends, parity } => {
                    s.central_atom = inchi_sys::NO_ATOM as i16;
                    s.neighbor = [
                        check(ends[0])?,
                        check(ends[1])?,
                        check(ends[2])?,
                        check(ends[3])?,
                    ];
                    s.type_ = inchi_sys::INCHI_StereoType_DoubleBond as i8;
                    s.parity = parity.code();
                }
                Stereo::Allene {
                    center,
                    ends,
                    parity,
                } => {
                    s.central_atom = check(center)?;
                    s.neighbor = [
                        check(ends[0])?,
                        check(ends[1])?,
                        check(ends[2])?,
                        check(ends[3])?,
                    ];
                    s.type_ = inchi_sys::INCHI_StereoType_Allene as i8;
                    s.parity = parity.code();
                }
            }
            raw.push(s);
        }
        Ok(raw)
    }
}

/// Owns every heap allocation behind an [`inchi_sys::inchi_Input_Polymer`] so
/// the C side sees stable pointers for the whole `GetINCHIEx` call.
///
/// Field order matters only for clarity; all pointers are taken after the
/// backing vectors are fully populated, and the inner buffers stay put for the
/// lifetime of the value.
struct PolymerBacking {
    // Per-unit 1-based atom lists (SAL) and crossing-bond lists (SBL).
    alists: Vec<Vec<std::os::raw::c_int>>,
    blists: Vec<Vec<std::os::raw::c_int>>,
    // The unit structs and the array of pointers to them.
    units: Vec<inchi_sys::inchi_Input_PolymerUnit>,
    unit_ptrs: Vec<*mut inchi_sys::inchi_Input_PolymerUnit>,
    polymer: inchi_sys::inchi_Input_Polymer,
}

impl PolymerBacking {
    fn build(units_in: &[crate::polymer::PolymerUnit], num_atoms: usize) -> Result<Box<Self>> {
        let one_based = |idx: usize| -> Result<std::os::raw::c_int> {
            if idx >= num_atoms {
                return Err(InchiError::InvalidStructure {
                    reason: format!(
                        "polymer unit references nonexistent atom {idx} (have {num_atoms})"
                    ),
                });
            }
            i32::try_from(idx + 1).map_err(|_| InchiError::InvalidStructure {
                reason: format!("atom index {idx} out of range"),
            })
        };

        let mut alists = Vec::with_capacity(units_in.len());
        let mut blists = Vec::with_capacity(units_in.len());
        for unit in units_in {
            let mut alist = Vec::with_capacity(unit.atoms.len());
            for &a in &unit.atoms {
                alist.push(one_based(a)?);
            }
            let mut blist = Vec::with_capacity(unit.crossing_bonds.len() * 2);
            for &[a, b] in &unit.crossing_bonds {
                blist.push(one_based(a)?);
                blist.push(one_based(b)?);
            }
            alists.push(alist);
            blists.push(blist);
        }

        // Allocate boxed so the struct's address (and thus every interior
        // pointer the C side stores) is stable even if the caller moves us.
        let mut me = Box::new(PolymerBacking {
            alists,
            blists,
            units: Vec::with_capacity(units_in.len()),
            unit_ptrs: Vec::with_capacity(units_in.len()),
            polymer: unsafe { std::mem::zeroed() },
        });

        // Build the unit structs into a local Vec, taking stable buffer
        // pointers from the already-populated `alists`/`blists` (which live in
        // `me`). Using a local Vec sidesteps simultaneous borrows of `me`.
        let mut built = Vec::with_capacity(units_in.len());
        let lists = me.alists.iter().zip(me.blists.iter());
        for (unit, (alist, blist)) in units_in.iter().zip(lists) {
            let mut raw: inchi_sys::inchi_Input_PolymerUnit = unsafe { std::mem::zeroed() };
            raw.id = unit.id;
            raw.label = unit.label;
            raw.type_ = unit.kind.code();
            raw.subtype = unit.subtype.code();
            raw.conn = unit.connection.code();
            raw.na = i32::try_from(unit.atoms.len()).unwrap_or(0);
            raw.nb = i32::try_from(unit.crossing_bonds.len()).unwrap_or(0);
            write_subscript(&mut raw.smt, &unit.subscript);
            // SAFETY: the inner buffers live as long as `me`; the C side never
            // mutates them, so casting the const buffer pointer to `*mut` is sound.
            raw.alist = if alist.is_empty() {
                std::ptr::null_mut()
            } else {
                alist.as_ptr() as *mut std::os::raw::c_int
            };
            raw.blist = if blist.is_empty() {
                std::ptr::null_mut()
            } else {
                blist.as_ptr() as *mut std::os::raw::c_int
            };
            built.push(raw);
        }
        me.units = built;

        // Record the array of pointers into the now-stable `units`.
        let mut ptrs = Vec::with_capacity(me.units.len());
        for u in me.units.iter_mut() {
            ptrs.push(u as *mut inchi_sys::inchi_Input_PolymerUnit);
        }
        me.unit_ptrs = ptrs;
        me.polymer.n = i32::try_from(me.unit_ptrs.len()).unwrap_or(0);
        me.polymer.units = me.unit_ptrs.as_mut_ptr();
        Ok(me)
    }

    fn as_mut_ptr(&mut self) -> *mut inchi_sys::inchi_Input_Polymer {
        &mut self.polymer
    }
}

/// Writes a polymer Sgroup subscript into the fixed-size `smt` field (80 bytes,
/// NUL-terminated), truncating if necessary.
fn write_subscript(dst: &mut [std::os::raw::c_char; 80], subscript: &str) {
    let max = dst.len().saturating_sub(1);
    for (slot, &b) in dst.iter_mut().zip(subscript.as_bytes().iter().take(max)) {
        *slot = b as std::os::raw::c_char;
    }
}

/// Records `to` as a neighbor of `from` in the raw atom adjacency list.
fn push_neighbor(
    atoms: &mut [inchi_sys::inchi_Atom],
    from: usize,
    to: usize,
    order: BondOrder,
) -> Result<()> {
    let to_idx = i16::try_from(to).map_err(|_| InchiError::InvalidStructure {
        reason: format!("atom index {to} out of range"),
    })?;
    let atom = atoms
        .get_mut(from)
        .ok_or_else(|| InchiError::InvalidStructure {
            reason: format!("bond references nonexistent atom {from}"),
        })?;
    let slot = atom.num_bonds as usize;
    if slot >= MAX_BONDS_PER_ATOM {
        return Err(InchiError::InvalidStructure {
            reason: format!("atom {from} exceeds the maximum of {MAX_BONDS_PER_ATOM} bonds"),
        });
    }
    if let (Some(nbr), Some(bt)) = (atom.neighbor.get_mut(slot), atom.bond_type.get_mut(slot)) {
        *nbr = to_idx;
        *bt = order.code();
        atom.num_bonds += 1;
        Ok(())
    } else {
        Err(InchiError::InvalidStructure {
            reason: format!("atom {from} bond slot {slot} out of range"),
        })
    }
}

/// Writes an element symbol into the fixed-size `elname` field, validating it.
fn write_elname(dst: &mut [std::os::raw::c_char; ELNAME_CAP], symbol: &str) -> Result<()> {
    let bytes = symbol.as_bytes();
    if bytes.is_empty() {
        return Err(InchiError::InvalidStructure {
            reason: "empty element symbol".to_string(),
        });
    }
    if !symbol.is_ascii() {
        return Err(InchiError::InvalidStructure {
            reason: format!("element symbol {symbol:?} is not ASCII"),
        });
    }
    if bytes.len() >= ELNAME_CAP {
        return Err(InchiError::InvalidStructure {
            reason: format!(
                "element symbol {symbol:?} is too long (max {} chars)",
                ELNAME_CAP - 1
            ),
        });
    }
    for (slot, &b) in dst.iter_mut().zip(bytes) {
        *slot = b as std::os::raw::c_char;
    }
    Ok(())
}
