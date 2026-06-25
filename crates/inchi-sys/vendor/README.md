# Vendored IUPAC InChI source

This directory contains a verbatim subset of the official IUPAC International Chemical Identifier (InChI) reference implementation, vendored so that `inchi-sys` builds with **no network access** and **no system InChI install**.

| Field    | Value                                      |
| -------- | ------------------------------------------ |
| Project  | InChI                                      |
| Upstream | <https://github.com/IUPAC-InChI/InChI>     |
| Version  | 1.07.3                                     |
| Git tag  | `v1.07.3`                                  |
| Commit   | `0b3e941d29289f3e5024c9ecfd45186285319420` |
| License  | MIT (see [`inchi/LICENSE`](inchi/LICENSE)) |

## What is included

Only the files required to build the InChI shared library (the "InChI API") are vendored, mirroring the upstream `INCHI_API/libinchi/gcc/makefile` build:

```
inchi/INCHI_BASE/src/             core algorithm + classic API header (inchi_api.h)
inchi/INCHI_API/libinchi/src/     DLL/API entry points (GetINCHI, FreeINCHI, ...)
inchi/INCHI_API/libinchi/src/ixa/ InChI eXtensible API (IXA)
```

The command-line front-end (`INCHI_EXE`), the GUI helper (`inchi_gui.c`), prebuilt binaries (`INCHI-1-BIN`), tests, and documentation from the upstream repository are intentionally **not** vendored.

## Updating

1. Download the desired release tarball from the upstream releases page.
2. Replace the three `src` trees above with the new versions.
3. Update the version/commit metadata in this file.
4. Run `cargo test -p inchi-sys` and `cargo test -p inchi` to confirm the build still links and the reference InChI/InChIKey test vectors still pass.
5. If the public C API changed, regenerate the bindings: `cargo build -p inchi-sys --features regenerate-bindings` and copy `target/.../out/bindings.rs` over `crates/inchi-sys/src/bindings.rs`.
6. Re-apply the patch listed below if the affected file changed upstream.

## Local patches

The vendored C sources are otherwise verbatim, with a single memory-safety fix. Each patched site is marked with an `inchi-rs patch:` comment so it can be located and re-applied after an upstream update.

| File | Function | Fix |
| ---- | -------- | --- |
| `inchi/INCHI_API/libinchi/src/inchi_dll.c` | `SetInChIExtInputByExtOrigAtData` | Removes a use-after-free: the upstream oss-fuzz hardening (issues #67695, #66748) freed the polymer data it had just assigned to the caller-owned `*iip`, so `GetStructFromINCHIEx` returned a dangling `polymer` pointer. The data is owned by the caller and released later by `FreeStructFromINCHIEx`/`FreeInChIExtInput`. The fix does **not** change any computed InChI string. |

This patch was reported upstream; remove it once a fixed release is vendored.
