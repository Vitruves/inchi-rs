//! Integration tests validating the public API against canonical
//! InChI/InChIKey reference values for real molecules.
//!
//! The Molfile fixtures in `tests/fixtures/` were generated with RDKit (with 3D
//! coordinates for the stereochemical cases). Every expected InChI/InChIKey here
//! was cross-checked against RDKit's bundled IUPAC InChI library and, for the
//! well-known compounds, against the public PubChem record.

use inchi::{
    check_inchi, check_inchikey, from_molfile, inchi_to_inchi, inchikey, inchikey_with_hashes,
    struct_from_aux_info, struct_from_inchi, struct_from_inchi_ex, Atom, BondOrder,
    InchiKeyValidity, InchiValidity, Molecule, Options, PolymerConnection, PolymerUnit,
    PolymerUnitKind, Polymers,
};

// --- Simple molecules (small, no stereo) ------------------------------------

const METHANE: &str = include_str!("fixtures/methane.mol");
const WATER: &str = include_str!("fixtures/water.mol");
const ETHANOL: &str = include_str!("fixtures/ethanol.mol");
const SODIUM_ACETATE: &str = include_str!("fixtures/sodium_acetate.mol");

// --- Complex, real-world drug-like molecules --------------------------------

const BENZENE: &str = include_str!("fixtures/benzene.mol");
const CAFFEINE: &str = include_str!("fixtures/caffeine.mol");
const ASPIRIN: &str = include_str!("fixtures/aspirin.mol");
const IBUPROFEN: &str = include_str!("fixtures/ibuprofen.mol");

// --- Stereochemically rich molecules ----------------------------------------

const L_ALANINE: &str = include_str!("fixtures/l_alanine.mol");
const GLUCOSE: &str = include_str!("fixtures/glucose.mol");
const L_TRYPTOPHAN: &str = include_str!("fixtures/l_tryptophan.mol");

/// `(molfile, expected InChI, expected InChIKey)` for every fixture.
const CASES: &[(&str, &str, &str)] = &[
    (METHANE, "InChI=1S/CH4/h1H4", "VNWKTOKETHGBQD-UHFFFAOYSA-N"),
    (WATER, "InChI=1S/H2O/h1H2", "XLYOFNOQVPJJNP-UHFFFAOYSA-N"),
    (
        ETHANOL,
        "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
    ),
    (
        SODIUM_ACETATE,
        "InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1",
        "VMHLLURERBWHNL-UHFFFAOYSA-M",
    ),
    (
        BENZENE,
        "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",
        "UHOVQNZJYSORNB-UHFFFAOYSA-N",
    ),
    (
        CAFFEINE,
        "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
        "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
    ),
    (
        ASPIRIN,
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
        "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    ),
    (
        IBUPROFEN,
        "InChI=1S/C13H18O2/c1-9(2)8-11-4-6-12(7-5-11)10(3)13(14)15/h4-7,9-10H,8H2,1-3H3,(H,14,15)",
        "HEFNNWSXXWATRW-UHFFFAOYSA-N",
    ),
    (
        L_ALANINE,
        "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1",
        "QNAYBMKLOCPYGJ-REOHCLBHSA-N",
    ),
    (
        GLUCOSE,
        "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6+/m1/s1",
        "WQZGKKKJIJFFOK-DVKNGEFBSA-N",
    ),
    (
        L_TRYPTOPHAN,
        "InChI=1S/C11H12N2O2/c12-9(11(14)15)5-7-6-13-10-4-2-1-3-8(7)10/h1-4,6,9,13H,5,12H2,(H,14,15)/t9-/m0/s1",
        "QIVBCDIJIAJPQS-VIFPVBQESA-N",
    ),
];

fn inchi_of(molfile: &str) -> String {
    from_molfile(molfile, Options::new())
        .unwrap_or_else(|e| panic!("generation failed: {e}"))
        .into_inchi()
}

/// Every fixture must reproduce its reference InChI exactly, and the InChIKey
/// derived from that InChI must match the reference key.
#[test]
fn forward_generation_matches_references() {
    for &(molfile, expected_inchi, expected_key) in CASES {
        let inchi = inchi_of(molfile);
        assert_eq!(inchi, expected_inchi, "InChI mismatch");
        let key = inchikey(&inchi).expect("key");
        assert_eq!(key, expected_key, "InChIKey mismatch for {expected_inchi}");
    }
}

/// The InChIKey computed from each reference InChI is stable and round-trips
/// through the key validator as a valid *standard* key.
#[test]
fn every_reference_key_is_valid() {
    for &(_, inchi, expected_key) in CASES {
        let key = inchikey(inchi).expect("key");
        assert_eq!(key, expected_key);
        assert_eq!(
            check_inchikey(&key).expect("check"),
            InchiKeyValidity::Standard
        );
    }
}

/// The programmatic builder must agree with the molfile path for benzene
/// (aromatic ring) and glycine (a zwitterion-free amino acid), down to the key.
#[test]
fn molecule_builder_matches_known_identifiers() {
    // Benzene: a six-membered carbon ring with alternating bond orders.
    let mut benzene = Molecule::new();
    let c: Vec<usize> = (0..6).map(|_| benzene.add_atom(Atom::new("C"))).collect();
    for i in 0..6 {
        let order = if i % 2 == 0 {
            BondOrder::Double
        } else {
            BondOrder::Single
        };
        benzene
            .add_bond(c[i], c[(i + 1) % 6], order)
            .expect("bond");
    }
    let inchi = benzene.to_inchi(Options::new()).expect("gen").into_inchi();
    assert_eq!(inchi, "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H");
    assert_eq!(inchikey(&inchi).expect("key"), "UHOVQNZJYSORNB-UHFFFAOYSA-N");

    // Glycine: H2N-CH2-C(=O)-OH.
    let mut gly = Molecule::new();
    let n = gly.add_atom(Atom::new("N"));
    let ca = gly.add_atom(Atom::new("C"));
    let cc = gly.add_atom(Atom::new("C"));
    let od = gly.add_atom(Atom::new("O"));
    let oh = gly.add_atom(Atom::new("O"));
    gly.add_bond(n, ca, BondOrder::Single).expect("bond");
    gly.add_bond(ca, cc, BondOrder::Single).expect("bond");
    gly.add_bond(cc, od, BondOrder::Double).expect("bond");
    gly.add_bond(cc, oh, BondOrder::Single).expect("bond");
    let inchi = gly.to_inchi(Options::new()).expect("gen").into_inchi();
    assert_eq!(inchi, "InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)");
    assert_eq!(inchikey(&inchi).expect("key"), "DHMQDGOQFOQNFH-UHFFFAOYSA-N");
}

/// `struct_from_inchi` must recover the right heavy-atom and bond counts for a
/// fused-ring drug molecule (caffeine: bicyclic purine, ten heavy atoms).
#[test]
fn parse_structure_of_caffeine() {
    let s = struct_from_inchi(
        "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
    )
    .expect("parse");
    // 8 carbons + 4 nitrogens + 2 oxygens = 14 heavy atoms.
    assert_eq!(s.atoms().len(), 14);
    let carbons = s.atoms().iter().filter(|a| a.element == "C").count();
    let nitrogens = s.atoms().iter().filter(|a| a.element == "N").count();
    let oxygens = s.atoms().iter().filter(|a| a.element == "O").count();
    assert_eq!((carbons, nitrogens, oxygens), (8, 4, 2));
    // The bicyclic purine skeleton has more bonds than heavy atoms minus one.
    assert!(s.bonds().len() >= s.atoms().len());
}

/// Stereocenters survive the parse: L-tryptophan has exactly one tetrahedral
/// center, and L-alanine likewise.
#[test]
fn parse_recovers_stereo() {
    let trp = struct_from_inchi(
        "InChI=1S/C11H12N2O2/c12-9(11(14)15)5-7-6-13-10-4-2-1-3-8(7)10/h1-4,6,9,13H,5,12H2,(H,14,15)/t9-/m0/s1",
    )
    .expect("parse");
    assert_eq!(trp.stereo().len(), 1);

    let ala =
        struct_from_inchi("InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1")
            .expect("parse");
    assert_eq!(ala.stereo().len(), 1);

    // Glucose carries five stereocenters in the pyranose ring.
    let glc = struct_from_inchi(
        "InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6+/m1/s1",
    )
    .expect("parse");
    assert_eq!(glc.stereo().len(), 5);
}

/// AuxInfo round-trip: parsing the auxiliary block of a freshly generated InChI
/// recovers the full hydrogen-complete atom set (aspirin, C9H8O4 = 21 atoms).
#[test]
fn aux_info_round_trip() {
    let out = from_molfile(ASPIRIN, Options::new()).expect("gen");
    assert!(out.aux_info().starts_with("AuxInfo="));

    // The fixture carries explicit hydrogens, so the AuxInfo records all 21
    // atoms (9 C + 8 H + 4 O); both hydrogen modes therefore recover them.
    let s = struct_from_aux_info(out.aux_info(), true).expect("parse");
    assert_eq!(s.atoms().len(), 21);
    assert_eq!(s.atoms().iter().filter(|a| a.element == "C").count(), 9);
    assert_eq!(s.atoms().iter().filter(|a| a.element == "O").count(), 4);
    assert_eq!(s.atoms().iter().filter(|a| a.element == "H").count(), 8);

    let no_add = struct_from_aux_info(out.aux_info(), false).expect("parse");
    assert_eq!(no_add.atoms().len(), 21);
}

/// `aux_info()` is present by default and suppressible via options.
#[test]
fn aux_info_present_and_suppressible() {
    let with = from_molfile(ETHANOL, Options::new()).expect("gen");
    assert!(with.aux_info().starts_with("AuxInfo="));

    let without = from_molfile(ETHANOL, Options::new().aux_info(false)).expect("gen");
    assert!(without.aux_info().is_empty());
}

/// `inchi_to_inchi` reproduces a standard InChI unchanged, and converts a
/// standard InChI into a non-standard FixedH form.
#[test]
fn inchi_to_inchi_normalizes_and_converts() {
    // Identity: standard in, standard out.
    for &(_, inchi, _) in CASES {
        let out = inchi_to_inchi(inchi, Options::new()).expect("convert");
        assert_eq!(out.inchi(), inchi);
    }

    // Standard -> non-standard via FixedH. Re-deriving from the *standard*
    // InChI cannot recover the mobile-H detail needed for an `/f` layer, so the
    // result is a non-standard InChI with the same layers under the `1/` prefix.
    let std = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)";
    let fixed = inchi_to_inchi(std, Options::new().fixed_h(true)).expect("convert");
    assert_eq!(
        fixed.inchi(),
        "InChI=1/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    );
    // The result is a valid *non-standard* InChI.
    assert_eq!(
        check_inchi(fixed.inchi(), false).expect("check"),
        InchiValidity::NonStandard
    );
}

/// `check_inchi` classifies standard, non-standard, and malformed inputs in
/// lenient mode, and faithfully reproduces the reference library's strict-mode
/// semantics (a round-trip through forced `FixedH` that only non-standard
/// InChIs can satisfy).
#[test]
fn check_inchi_classifies_validity() {
    // Lenient mode: every reference InChI is a valid standard InChI.
    for &(_, inchi, _) in CASES {
        assert_eq!(
            check_inchi(inchi, false).expect("check"),
            InchiValidity::Standard
        );
    }

    // Strict mode forces a FixedH round-trip whose output is non-standard, so a
    // standard InChI can never match itself: this is the documented reference
    // behavior, reproduced faithfully here.
    for &(_, inchi, _) in CASES {
        assert_eq!(
            check_inchi(inchi, true).expect("strict check"),
            InchiValidity::FailedRoundtrip
        );
    }

    // A non-standard InChI that re-derives to itself *is* accepted under strict
    // mode (the case strict checking is actually for).
    let nonstd = "InChI=1/CH4/h1H4";
    assert_eq!(
        check_inchi(nonstd, false).expect("check"),
        InchiValidity::NonStandard
    );
    assert_eq!(
        check_inchi(nonstd, true).expect("strict check"),
        InchiValidity::NonStandard
    );

    // Missing/garbled prefix and malformed layout are rejected.
    assert_eq!(
        check_inchi("not an inchi", false).expect("check"),
        InchiValidity::InvalidPrefix
    );
    assert!(!check_inchi("InChI=1S/garbage layout", false)
        .expect("check")
        .is_valid());
}

/// `check_inchikey` accepts genuine keys and rejects malformed ones.
#[test]
fn check_inchikey_classifies_validity() {
    assert_eq!(
        check_inchikey("RYYVLZVUVIJVGH-UHFFFAOYSA-N").expect("check"),
        InchiKeyValidity::Standard
    );
    // Standard-formula key with a deliberately wrong check character.
    assert_eq!(
        check_inchikey("TOOSHORT").expect("check"),
        InchiKeyValidity::InvalidLength
    );
    assert!(!check_inchikey("RYYVLZVUVIJVGH/UHFFFAOYSA/N")
        .expect("check")
        .is_valid());
}

/// `inchikey_with_hashes` returns the same key plus the requested 256-bit
/// extension blocks, and omits the blocks that were not requested.
#[test]
fn inchikey_extra_hashes() {
    let inchi = "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3";

    let none = inchikey_with_hashes(inchi, false, false).expect("key");
    assert_eq!(none.key, "RYYVLZVUVIJVGH-UHFFFAOYSA-N");
    assert!(none.extra1.is_none() && none.extra2.is_none());

    let both = inchikey_with_hashes(inchi, true, true).expect("key");
    assert_eq!(both.key, "RYYVLZVUVIJVGH-UHFFFAOYSA-N");
    let e1 = both.extra1.expect("extra1");
    let e2 = both.extra2.expect("extra2");
    assert!(!e1.is_empty() && e1.len() <= 64 && e1.chars().all(|c| c.is_ascii_hexdigit()));
    assert!(!e2.is_empty() && e2.len() <= 64 && e2.chars().all(|c| c.is_ascii_hexdigit()));
    assert_ne!(e1, e2);

    // Only one block requested.
    let only1 = inchikey_with_hashes(inchi, true, false).expect("key");
    assert!(only1.extra1.is_some() && only1.extra2.is_none());
}

const POLYMER_SRU: &str = include_str!("fixtures/polymer_sru.mol");

/// The polyethylene SRU reference InChI, shared across the polymer tests. It is
/// a beta-flagged, non-standard InChI carrying a `/z` polymer layer.
const POLYETHYLENE_INCHI: &str = "InChI=1B/C2H4Zz2/c3-1-2-4/h1-2H2/z101-1-2(4-2,3-1)";

/// Polymer generation from a Molfile carrying a CTFile SRU Sgroup yields the
/// expected beta InChI, and only when polymer processing is enabled.
#[test]
fn polymer_generation_from_molfile() {
    let out = from_molfile(POLYMER_SRU, Options::new().polymers(Polymers::On)).expect("gen");
    assert_eq!(out.inchi(), POLYETHYLENE_INCHI);
    assert!(out.inchi().contains("/z"), "expected a polymer layer");

    // Without the polymer flag the star atoms are rejected.
    assert!(from_molfile(POLYMER_SRU, Options::new()).is_err());
}

/// The programmatic builder reproduces the Molfile polymer InChI exactly via the
/// extended `GetINCHIEx` entry point.
#[test]
fn polymer_generation_from_builder() {
    let mut mol = Molecule::new();
    let s1 = mol.add_atom(Atom::new("Zz"));
    let c1 = mol.add_atom(Atom::new("C"));
    let c2 = mol.add_atom(Atom::new("C"));
    let s2 = mol.add_atom(Atom::new("Zz"));
    mol.add_bond(s1, c1, BondOrder::Single).expect("bond");
    mol.add_bond(c1, c2, BondOrder::Single).expect("bond");
    mol.add_bond(c2, s2, BondOrder::Single).expect("bond");
    mol.add_polymer_unit(PolymerUnit::sru([c1, c2], [[s1, c1], [c2, s2]]));

    let inchi = mol
        .to_inchi(Options::new().polymers(Polymers::On))
        .expect("gen")
        .into_inchi();
    assert_eq!(inchi, POLYETHYLENE_INCHI);
}

/// `struct_from_inchi_ex` recovers the polymer repeating unit from a polymer
/// InChI, and reports none for an ordinary InChI.
#[test]
fn polymer_round_trip_parse() {
    let ext = struct_from_inchi_ex(POLYETHYLENE_INCHI).expect("parse");
    assert_eq!(ext.structure.atoms().len(), 4); // 2 carbons + 2 star atoms
    assert_eq!(ext.polymer_units.len(), 1);

    let unit = &ext.polymer_units[0];
    assert_eq!(unit.kind, PolymerUnitKind::StructureBasedRepeat);
    assert_eq!(unit.connection, PolymerConnection::HeadToTail);
    assert_eq!(unit.atoms.len(), 2); // the two backbone carbons
    assert_eq!(unit.crossing_bonds.len(), 2); // the two bonds to the star atoms

    // An ordinary InChI parses with no polymer data.
    let plain = struct_from_inchi_ex("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3").expect("parse");
    assert!(plain.polymer_units.is_empty());
    assert_eq!(plain.structure, struct_from_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3").unwrap());
}

/// Errors are surfaced, not panicked, for malformed input across entry points.
#[test]
fn invalid_input_is_error() {
    assert!(struct_from_inchi("not an inchi").is_err());
    assert!(inchikey("not an inchi").is_err());
    assert!(inchi_to_inchi("not an inchi", Options::new()).is_err());
    assert!(struct_from_aux_info("not aux info", true).is_err());
    // An empty molecule cannot be converted.
    assert!(Molecule::new().to_inchi(Options::new()).is_err());
}
