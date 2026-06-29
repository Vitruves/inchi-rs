# inchi-sys

[![crates.io](https://img.shields.io/crates/v/inchi-sys.svg)](https://crates.io/crates/inchi-sys) [![docs.rs](https://docs.rs/inchi-sys/badge.svg)](https://docs.rs/inchi-sys)

Low-level `-sys` FFI bindings to the vendored IUPAC [InChI](https://www.inchi-trust.org/) 1.07.3 reference C library. The C source is bundled and statically linked (`links = "inchi"`) via a `cc` build script — no system install, no network at build time.

## Use the `inchi` crate instead

This crate exposes only **raw, `unsafe` C bindings** (`GetINCHI`, `GetStructFromINCHI`, `GetINCHIKeyFromINCHI`, `MakeINCHIFromMolfileText`, the `*Ex` polymer entry points, the `IXA_*` extensible-API functions, …) with no safety, RAII, or ergonomics.

For almost all uses, depend on the safe, idiomatic wrapper:

```toml
[dependencies]
inchi = "0.1.4"
```

→ **<https://crates.io/crates/inchi>** · **<https://docs.rs/inchi>**

Reach for `inchi-sys` directly only when you need an entry point the safe crate does not yet wrap.

## Bindings

Pre-generated bindings (`src/bindings.rs`) are committed, so the default build needs no `libclang`. Regenerate them from the vendored headers with:

```sh
cargo build -p inchi-sys --features regenerate-bindings
```

## Vendored source & patches

The InChI C sources live under [`vendor/`](vendor/README.md), which documents the exact upstream version, what is included, how to update, and the single local memory-safety patch applied to `GetStructFromINCHIEx`.

## License

MIT. Bundles the IUPAC InChI software (MIT); see [`vendor/inchi/LICENSE`](vendor/inchi/LICENSE).
