//! User guide — concepts and worked examples.
//!
//! The item-level documentation is the API *reference*; this module is the
//! narrative *guide*. Every example here is compiled and run as a doctest, so
//! it cannot drift from the code.
//!
//! - [Standard vs. non-standard InChI](#standard-vs-non-standard-inchi)
//! - [InChIKeys](#inchikeys)
//! - [Stereochemistry (0D parities)](#stereochemistry-0d-parities)
//! - [Polymers](#polymers)
//! - [Thread safety and reuse](#thread-safety-and-reuse)
//!
//! # Standard vs. non-standard InChI
//!
//! A *standard* InChI (prefix `InChI=1S/…`) is the canonical, interoperable
//! identifier produced with default normalization. Any option that changes that
//! normalization — a fixed-hydrogen layer ([`Options::fixed_h`]), metal
//! reconnection ([`Options::reconnect_metals`]), non-default stereo
//! ([`StereoMode::Relative`] / [`StereoMode::Racemic`] / [`StereoMode::UseChiralFlag`]),
//! the tautomer extensions, or polymers — yields a *non-standard* InChI
//! (prefix `InChI=1/…`). The two are **not** interchangeable and hash to
//! different InChIKeys.
//!
//! [`Options::is_standard`] tells you ahead of time which you will get, and
//! [`check_inchi`](crate::check_inchi) classifies an existing string:
//!
//! ```
//! use inchi::{check_inchi, inchi_to_inchi, InchiValidity, Options};
//!
//! let standard = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3";
//! assert!(Options::new().is_standard());
//! assert_eq!(check_inchi(standard, false)?, InchiValidity::Standard);
//!
//! // Re-derive it with a fixed-H layer: now it is non-standard.
//! let opts = Options::new().fixed_h(true);
//! assert!(!opts.is_standard());
//! let non_standard = inchi_to_inchi(standard, opts)?;
//! assert!(non_standard.inchi().starts_with("InChI=1/"));
//! assert_eq!(check_inchi(non_standard.inchi(), false)?, InchiValidity::NonStandard);
//! # Ok::<(), inchi::InchiError>(())
//! ```
//!
//! Note that [`check_inchi`](crate::check_inchi)'s *strict* mode is a different,
//! narrower check (a forced `FixedH` round-trip) — see its own documentation;
//! for "is this a valid InChI?" use the lenient mode shown above.
//!
//! # InChIKeys
//!
//! An InChIKey is a fixed-length (27-character) hashed form of an InChI, suited
//! to database keys and web search. The hash is one-way: you cannot recover the
//! InChI from it. [`inchikey`](crate::inchikey) computes it; the optional 256-bit
//! hash extensions (rarely needed) come from
//! [`inchikey_with_hashes`](crate::inchikey_with_hashes):
//!
//! ```
//! use inchi::{check_inchikey, inchikey, inchikey_with_hashes, InchiKeyValidity};
//!
//! let key = inchikey("InChI=1S/CH4/h1H4")?;
//! assert_eq!(key, "VNWKTOKETHGBQD-UHFFFAOYSA-N");
//! assert_eq!(check_inchikey(&key)?, InchiKeyValidity::Standard);
//!
//! let extended = inchikey_with_hashes("InChI=1S/CH4/h1H4", true, false)?;
//! assert_eq!(extended.key, key);
//! assert!(extended.extra1.is_some());
//! # Ok::<(), inchi::InchiError>(())
//! ```
//!
//! # Stereochemistry (0D parities)
//!
//! When a structure has no usable coordinates, configuration is conveyed with
//! *0D parities*: each [`Stereo`](crate::Stereo) element lists the neighbor
//! atoms in a fixed order, and the [`Parity`](crate::Parity) is the handedness
//! of that ordering. The orderings (taken verbatim from the InChI specification)
//! are:
//!
//! - [`Stereo::Tetrahedral`](crate::Stereo::Tetrahedral) — `neighbors = [W, X, Y, Z]`
//!   around `center`; parity is even (`+`) if `X, Y, Z` appear clockwise when
//!   viewed from `W`. For a center with only three explicit neighbors plus an
//!   implicit hydrogen, repeat the `center` index in the slot where the
//!   hydrogen would go.
//! - [`Stereo::DoubleBond`](crate::Stereo::DoubleBond) — `ends = [X, A, B, Y]`
//!   for `X–A=B–Y`.
//! - [`Stereo::Allene`](crate::Stereo::Allene) — `ends = [X, A, B, Y]` with the
//!   central cumulene atom in `center`.
//!
//! Getting the order wrong flips the parity, so always validate against a known
//! identifier. Building bromochlorofluoromethane (`CHFClBr`, a textbook single
//! stereocenter) with the implicit H represented by the repeated center index:
//!
//! ```
//! use inchi::{Atom, BondOrder, Molecule, Parity, Stereo};
//!
//! let mut mol = Molecule::new();
//! let c = mol.add_atom(Atom::new("C"));
//! let f = mol.add_atom(Atom::new("F"));
//! let cl = mol.add_atom(Atom::new("Cl"));
//! let br = mol.add_atom(Atom::new("Br"));
//! mol.add_bond(c, f, BondOrder::Single)?;
//! mol.add_bond(c, cl, BondOrder::Single)?;
//! mol.add_bond(c, br, BondOrder::Single)?;
//! // The 4th neighbor is the center itself, standing in for the implicit H.
//! mol.add_stereo(Stereo::Tetrahedral { center: c, neighbors: [f, cl, br, c], parity: Parity::Odd });
//!
//! let inchi = mol.to_inchi(())?.into_inchi();
//! assert_eq!(inchi, "InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1"); // note the /t/m/s stereo layers
//! # Ok::<(), inchi::InchiError>(())
//! ```
//!
//! The reverse direction recovers the same elements; the neighbor indices refer
//! to the InChI's canonical atom numbering, not your input order:
//!
//! ```
//! use inchi::{struct_from_inchi, Stereo};
//!
//! let s = struct_from_inchi("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1")?; // L-alanine
//! assert!(s.stereo().iter().any(|e| matches!(e, Stereo::Tetrahedral { .. })));
//! # Ok::<(), inchi::InchiError>(())
//! ```
//!
//! # Polymers
//!
//! Polymers are an experimental, non-standard extension: a structural repeating
//! unit (SRU) is described by the atoms it spans and the bonds that cross its
//! brackets, and the open ends are capped with `"Zz"` star atoms. Enable it with
//! [`Options::polymers`] and describe the unit with [`PolymerUnit`](crate::PolymerUnit).
//! The result carries a `B` ("beta") version flag and a `/z` layer.
//!
//! Building a polyethylene repeat unit `*–CH2–CH2–*` programmatically:
//!
//! ```
//! use inchi::{Atom, BondOrder, Molecule, Options, Polymers, PolymerUnit};
//!
//! let mut mol = Molecule::new();
//! let s1 = mol.add_atom(Atom::new("Zz"));
//! let c1 = mol.add_atom(Atom::new("C"));
//! let c2 = mol.add_atom(Atom::new("C"));
//! let s2 = mol.add_atom(Atom::new("Zz"));
//! mol.add_bond(s1, c1, BondOrder::Single)?;
//! mol.add_bond(c1, c2, BondOrder::Single)?;
//! mol.add_bond(c2, s2, BondOrder::Single)?;
//! // The unit spans the two carbons; its crossing bonds are the two to the stars.
//! mol.add_polymer_unit(PolymerUnit::sru([c1, c2], [[s1, c1], [c2, s2]]));
//!
//! let inchi = mol.to_inchi(Options::new().polymers(Polymers::On))?.into_inchi();
//! assert_eq!(inchi, "InChI=1B/C2H4Zz2/c3-1-2-4/h1-2H2/z101-1-2(4-2,3-1)");
//! # Ok::<(), inchi::InchiError>(())
//! ```
//!
//! [`struct_from_inchi_ex`](crate::struct_from_inchi_ex) reverses this, recovering
//! the [`PolymerUnit`](crate::PolymerUnit)s (it returns an empty unit list for
//! ordinary InChIs):
//!
//! ```
//! use inchi::struct_from_inchi_ex;
//!
//! let ext = struct_from_inchi_ex("InChI=1B/C2H4Zz2/c3-1-2-4/h1-2H2/z101-1-2(4-2,3-1)")?;
//! assert_eq!(ext.polymer_units.len(), 1);
//! assert_eq!(ext.polymer_units[0].atoms.len(), 2);
//! # Ok::<(), inchi::InchiError>(())
//! ```
//!
//! # Queryable molecules (IXA)
//!
//! [`struct_from_inchi`](crate::struct_from_inchi) hands back a structure as
//! fixed arrays. When you want a *live* object to interrogate atom-by-atom — or
//! to load a Molfile once and both inspect it and (re)generate its InChI — reach
//! for [`IxaMolecule`](crate::IxaMolecule), built on the InChI eXtensible API.
//!
//! ```
//! use inchi::{IxaMolecule, BondOrder};
//!
//! // Reconstruct ethanol from its InChI, then walk its heavy-atom skeleton.
//! let mol = IxaMolecule::from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")?;
//! assert_eq!(mol.atom_count(), 3); // C, C, O — hydrogens stay implicit
//!
//! for bond in mol.bonds() {
//!     assert_eq!(mol.bond_order(bond)?, BondOrder::Single);
//! }
//!
//! // Regenerate the identifier (and the InChIKey) straight from the object.
//! assert_eq!(mol.to_inchi(())?.inchi(), "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3");
//! assert_eq!(mol.to_inchikey(())?, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N");
//! # Ok::<(), inchi::InchiError>(())
//! ```
//!
//! Atom positions are zero for InChI-derived molecules (an InChI stores no
//! geometry); use [`IxaMolecule::from_molfile`](crate::IxaMolecule::from_molfile)
//! when you need the original coordinates. Stereo descriptors are populated when
//! reading an InChI (encoded in the `/t` layer) but not from a Molfile, where
//! stereochemistry is perceived later during generation.
//!
//! # Thread safety and reuse
//!
//! The InChI C library keeps internal `static` state and is not re-entrant, so
//! every call in this crate is serialized behind a global lock. Calls are thus
//! safe from multiple threads but never run concurrently; for throughput,
//! parallelize the *work around* InChI generation (parsing, I/O) rather than
//! expecting the calls themselves to overlap. The public types are
//! [`Send`] + [`Sync`], and an [`Options`] value is cheap to clone and reuse
//! across many calls.
//!
//! [`Options`]: crate::Options
//! [`Options::fixed_h`]: crate::Options::fixed_h
//! [`Options::reconnect_metals`]: crate::Options::reconnect_metals
//! [`Options::is_standard`]: crate::Options::is_standard
//! [`Options::polymers`]: crate::Options::polymers
//! [`StereoMode::Relative`]: crate::StereoMode::Relative
//! [`StereoMode::Racemic`]: crate::StereoMode::Racemic
//! [`StereoMode::UseChiralFlag`]: crate::StereoMode::UseChiralFlag
