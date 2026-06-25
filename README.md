# inchi-rs

[![CI](https://github.com/Vitruves/inchi-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/Vitruves/inchi-rs/actions/workflows/ci.yml) [![crates.io](https://img.shields.io/crates/v/inchi.svg)](https://crates.io/crates/inchi) [![docs.rs](https://docs.rs/inchi/badge.svg)](https://docs.rs/inchi) [![license](https://img.shields.io/badge/license-MIT-blue.svg)](#license)

Rust bindings to the IUPAC [InChI](https://www.inchi-trust.org/) reference library, for generating and parsing **InChI** and **InChIKey** chemical identifiers. The InChI 1.07.3 C library is vendored and statically linked — no system install, no network at build time.

## Crates

| Crate | Description | Use it when |
| ----- | ----------- | ----------- |
| [**`inchi`**](crates/inchi) [![docs.rs](https://docs.rs/inchi/badge.svg)](https://docs.rs/inchi) | Safe, idiomatic API — `Result`-based, panic-free, no `unsafe`, automatic cleanup. | Almost always. |
| [`inchi-sys`](crates/inchi-sys) [![docs.rs](https://docs.rs/inchi-sys/badge.svg)](https://docs.rs/inchi-sys) | Raw `-sys` FFI bindings + the vendored C source. | Only for an entry point the safe crate does not wrap. |

## Quick start

```toml
[dependencies]
inchi = "0.1"
```

```rust
use inchi::{from_molfile, inchikey};

let out = from_molfile(molfile, ())?;     // () = default options
let inchi = out.inchi();                  // "InChI=1S/H2O/h1H2"
let key = inchikey(inchi)?;               // "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
```

Build a structure programmatically instead:

```rust
use inchi::{Molecule, Atom, BondOrder};

let mut mol = Molecule::new();
let c = mol.add_atom(Atom::new("C"));
let o = mol.add_atom(Atom::new("O"));
mol.add_bond(c, o, BondOrder::Double)?;
let inchi = mol.to_inchi(())?.into_inchi();
```

## What's covered

The full InChI generation, parsing, and validation API:

| Direction | Functions |
| --------- | --------- |
| Structure → InChI | `from_molfile`, `Molecule::to_inchi` |
| InChI → structure | `struct_from_inchi`, `struct_from_inchi_ex` (polymers), `struct_from_aux_info` |
| InChI → InChIKey | `inchikey`, `inchikey_with_hashes` |
| Convert / validate | `inchi_to_inchi`, `check_inchi`, `check_inchikey` |

Generation is configured with `Options` (stereo, fixed-H, metals, tautomers, timeout, polymers). Polymers (CTFile SRU Sgroups) are supported in both directions via `Options::polymers` / `Molecule::add_polymer_unit`.

Full reference and runnable examples: **<https://docs.rs/inchi>**.

## Notes

- The InChI C library keeps `static` state and is not re-entrant; all calls are serialized behind a global lock. Public types are `Send + Sync`.
- MSRV: Rust 1.77.
- Vendored InChI version, included files, and the single local memory-safety patch are documented in [`crates/inchi-sys/vendor/README.md`](crates/inchi-sys/vendor/README.md).

## Development

```sh
cargo test                       # safe crate + sys crate + doctests + reference vectors
cargo clippy --all-targets       # lints (the safe crate denies unwrap/panic/indexing)
```

## License

MIT. Bundles the IUPAC InChI software, also MIT licensed; see [`crates/inchi-sys/vendor/inchi/LICENSE`](crates/inchi-sys/vendor/inchi/LICENSE).
