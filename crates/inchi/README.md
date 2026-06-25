# inchi

[![crates.io](https://img.shields.io/crates/v/inchi.svg)](https://crates.io/crates/inchi) [![docs.rs](https://docs.rs/inchi/badge.svg)](https://docs.rs/inchi)

Safe, idiomatic Rust bindings to the IUPAC [InChI](https://www.inchi-trust.org/) reference library, for generating and parsing **InChI** and **InChIKey** chemical identifiers. The C library is vendored and statically linked via [`inchi-sys`](https://crates.io/crates/inchi-sys) — no system install, no network at build time. No `unsafe` is exposed; native allocations are freed automatically.

## Install

```toml
[dependencies]
inchi = "0.1"
```

## Example

From a Molfile/SDF record (pass `()` for default options):

```rust
use inchi::{from_molfile, inchikey};

let out = from_molfile(molfile, ())?;
let inchi = out.inchi();        // "InChI=1S/H2O/h1H2"
let key = inchikey(inchi)?;     // "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
```

Or build a structure programmatically:

```rust
use inchi::{Molecule, Atom, BondOrder};

let mut mol = Molecule::new();
let c = mol.add_atom(Atom::new("C"));
let o = mol.add_atom(Atom::new("O"));
mol.add_bond(c, o, BondOrder::Double)?;
let inchi = mol.to_inchi(())?.into_inchi();
```

Customize generation by passing [`Options`] instead of `()`:

```rust
use inchi::{from_molfile, Options, StereoMode};

let out = from_molfile(molfile, Options::new().fixed_h(true).stereo(StereoMode::Relative))?;
```

## API

| Direction | Entry points |
| --------- | ------------ |
| Structure → InChI | [`from_molfile`], [`Molecule::to_inchi`] |
| InChI → structure | [`struct_from_inchi`], [`struct_from_inchi_ex`] (polymers), [`struct_from_aux_info`] |
| InChI → InChIKey | [`inchikey`], [`inchikey_with_hashes`] |
| Convert / validate | [`inchi_to_inchi`], [`check_inchi`], [`check_inchikey`] |

Generation is configured with [`Options`] (stereo, fixed-H, metals, tautomers, timeout, polymers, …). Polymers (CTFile SRU Sgroups) are supported via [`Options::polymers`] and [`Molecule::add_polymer_unit`].

Full reference, semantics, and runnable examples: **<https://docs.rs/inchi>**.

## Features

| Feature | Default | Description |
| ------- | ------- | ----------- |
| `key` | yes | InChIKey computation (`GetINCHIKeyFromINCHI`). |
| `regenerate-bindings` | no | Regenerate FFI bindings with `bindgen` (needs `libclang`). |

## Notes

- The InChI C library keeps `static` state and is not re-entrant; all calls are serialized behind a global lock. The public types are `Send + Sync`.
- MSRV: Rust 1.77. Raising it is a semver-breaking change.

## License

MIT. Bundles the IUPAC InChI software (MIT); see [`inchi-sys`](https://crates.io/crates/inchi-sys).

[`from_molfile`]: https://docs.rs/inchi/latest/inchi/fn.from_molfile.html
[`Molecule::to_inchi`]: https://docs.rs/inchi/latest/inchi/struct.Molecule.html#method.to_inchi
[`Molecule::add_polymer_unit`]: https://docs.rs/inchi/latest/inchi/struct.Molecule.html#method.add_polymer_unit
[`struct_from_inchi`]: https://docs.rs/inchi/latest/inchi/fn.struct_from_inchi.html
[`struct_from_inchi_ex`]: https://docs.rs/inchi/latest/inchi/fn.struct_from_inchi_ex.html
[`struct_from_aux_info`]: https://docs.rs/inchi/latest/inchi/fn.struct_from_aux_info.html
[`inchikey`]: https://docs.rs/inchi/latest/inchi/fn.inchikey.html
[`inchikey_with_hashes`]: https://docs.rs/inchi/latest/inchi/fn.inchikey_with_hashes.html
[`inchi_to_inchi`]: https://docs.rs/inchi/latest/inchi/fn.inchi_to_inchi.html
[`check_inchi`]: https://docs.rs/inchi/latest/inchi/fn.check_inchi.html
[`check_inchikey`]: https://docs.rs/inchi/latest/inchi/fn.check_inchikey.html
[`Options`]: https://docs.rs/inchi/latest/inchi/struct.Options.html
[`Options::polymers`]: https://docs.rs/inchi/latest/inchi/struct.Options.html#method.polymers
