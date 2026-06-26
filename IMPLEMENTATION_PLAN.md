# Implementation Plan — Complete InChI C API Coverage

This document maps the full IUPAC InChI 1.07.3 C library surface onto a safe, idiomatic Rust API. It explains which parts of the C API are already wrapped, which are worth adding, which should be skipped and why, and gives the exact Rust design for each addition.

---

## 1. C API inventory

The InChI C library exposes three distinct API families.

### 1.1 Classic API (~35 functions in `inchi_api.h`)

| Group | C functions | Status |
|---|---|---|
| Structure → InChI | `GetINCHI`, `GetStdINCHI`, `GetINCHIEx` | Wrapped (`Molecule::to_inchi`, `from_molfile`) |
| InChI → Structure | `GetStructFromINCHI`, `GetStructFromStdINCHI`, `GetStructFromINCHIEx` | Wrapped (`struct_from_inchi`, `struct_from_inchi_ex`) |
| InChI ↔ InChI | `GetINCHIfromINCHI` | Wrapped (`inchi_to_inchi`) |
| AuxInfo → Structure | `Get_inchi_Input_FromAuxInfo`, `Get_std_inchi_Input_FromAuxInfo` | Wrapped (`struct_from_aux_info`) |
| InChIKey | `GetINCHIKeyFromINCHI`, `GetStdINCHIKeyFromStdINCHI`, `CheckINCHIKey` | Wrapped (`inchikey`, `inchikey_with_hashes`, `check_inchikey`) |
| Validation | `CheckINCHI` | Wrapped (`check_inchi`) |
| Molfile shortcut | `MakeINCHIFromMolfileText` | Wrapped (`from_molfile`, `from_sdf`) |
| Deallocation | `FreeINCHI`, `FreeStdINCHI`, `FreeStructFromINCHI`, `FreeStructFromStdINCHI`, `FreeStructFromINCHIEx`, `Free_inchi_Input`, `Free_std_inchi_Input` | Internal only (RAII via `OutputGuard`) |
| Utility | `GetStringLength` | Internal only |
| Std* duplicates | `GetStdINCHI`, `GetStdINCHIKeyFromStdINCHI`, `GetStructFromStdINCHI`, `Get_std_inchi_Input_FromAuxInfo`, `Free_std_inchi_Input`, `FreeStdINCHI`, `FreeStructFromStdINCHI` | **Skip** — see §3.1 |

### 1.2 Modularized generator (~14 functions in `inchi_api.h`)

`INCHIGEN_Create/Setup/DoNormalization/DoCanonicalization/DoSerialization/Reset/Destroy` and their `STDINCHIGEN_*` equivalents. Not currently bound.

### 1.3 IXA extensible API (~90 functions in `ixa.h`)

Object-oriented, handle-based. Not currently bound. Handles: `IXA_STATUS_HANDLE`, `IXA_MOL_HANDLE`, `IXA_INCHIBUILDER_HANDLE`, `IXA_INCHIKEYBUILDER_HANDLE`. Atom/bond/stereo IDs: `IXA_ATOMID`, `IXA_BONDID`, `IXA_STEREOID`, `IXA_POLYMERUNITID`.

---

## 2. Design constraints (apply to all phases)

**Thread safety.** The InChI C library is not re-entrant; all calls must be serialized through the existing `static INCHI_LOCK: Mutex<()>` in `raw.rs`. Every new FFI call acquires this lock before crossing the C boundary and holds it for the full duration of the call. No exceptions.

**Memory safety.** No `unsafe` at the public boundary. All allocations owned by C must be released via their corresponding `Free*` or `Destroy` functions, wrapped in Rust RAII types whose `Drop` calls the appropriate cleanup. Raw pointers are never stored in public types.

**Error model.** Rust functions return `Result<T, InchiError>`. C integer return codes are mapped via `Status::from_code`; IXA status handles are checked immediately after each call and converted to `Err(InchiError::Failed { ... })` if `IXA_STATUS_HasError` is true.

**Panic freedom.** The existing `deny!(unwrap_used, expect_used, panic, indexing_slicing)` lints remain active for all new code.

**No redundancy.** New types do not duplicate `Molecule`'s builder API. New types fill genuine gaps.

---

## 3. Non-goals

### 3.1 `Std*` variants — skip

`GetStdINCHI`, `GetStdINCHIKeyFromStdINCHI`, `GetStructFromStdINCHI` etc. are convenience wrappers in the C library that call the same code paths as their non-`Std` counterparts with the options locked to standard mode. In Rust, `Options::is_standard()` already covers this constraint. Exposing these would double the public API surface with zero new capability.

### 3.2 IXA molecule *building* — skip

`IXA_MOL_CreateAtom`, `IXA_MOL_SetAtomElement`, `IXA_MOL_CreateBond`, `IXA_MOL_SetBondType`, etc. fully duplicate `Molecule::add_atom` / `Molecule::add_bond`. The existing `Molecule` builder is already idiomatic and well-tested. Wrapping the IXA building API would add ~40 functions of FFI noise for zero gain.

### 3.3 INCHIGEN intermediate data — skip (for now)

`INCHIGEN_DATA` contains `NORM_ATOMS` arrays which are deeply internal C structures with no stable public layout guarantees. The normalization and canonicalization steps produce no user-observable data beyond what `GetINCHI` already returns. Binding `NORM_ATOMS` would require hundreds of lines of layout-sensitive unsafe code for a use case (step-through debugging of the InChI algorithm) that has no practical cheminformatics value. Mark as a future possibility if someone files a concrete use case.

---

## 4. Phase 1 — IXA bindings in `inchi-sys`

**Scope:** Add the IXA function declarations to `crates/inchi-sys/src/bindings.rs` and the INCHIGEN declarations.

**Approach:** Manually extend the existing `extern "C"` block in `bindings.rs` rather than regenerating with `bindgen`. The IXA handles are opaque sentinel types (a pattern of `typedef struct { int dummy; } Foo_STRUCT; typedef Foo_STRUCT *Foo;`) — they map to Rust `#[repr(C)] pub struct FooStruct { pub dummy: c_int }` with `pub type Foo = *mut FooStruct`.

**New types to add to `bindings.rs`:**

```rust
// Opaque sentinel handles — values are cast from int by the C runtime
#[repr(C)] pub struct IXA_STATUS_HANDLE_STRUCT   { pub dummy: ::std::os::raw::c_int }
#[repr(C)] pub struct IXA_MOL_HANDLE_STRUCT       { pub dummy: ::std::os::raw::c_int }
#[repr(C)] pub struct IXA_INCHIBUILDER_HANDLE_STRUCT { pub dummy: ::std::os::raw::c_int }
#[repr(C)] pub struct IXA_INCHIKEYBUILDER_HANDLE_STRUCT { pub dummy: ::std::os::raw::c_int }
#[repr(C)] pub struct IXA_ATOMID_STRUCT           { pub dummy: ::std::os::raw::c_int }
#[repr(C)] pub struct IXA_BONDID_STRUCT           { pub dummy: ::std::os::raw::c_int }
#[repr(C)] pub struct IXA_STEREOID_STRUCT         { pub dummy: ::std::os::raw::c_int }
#[repr(C)] pub struct IXA_POLYMERUNITID_STRUCT    { pub dummy: ::std::os::raw::c_int }

pub type IXA_STATUS_HANDLE          = *mut IXA_STATUS_HANDLE_STRUCT;
pub type IXA_MOL_HANDLE             = *mut IXA_MOL_HANDLE_STRUCT;
pub type IXA_INCHIBUILDER_HANDLE    = *mut IXA_INCHIBUILDER_HANDLE_STRUCT;
pub type IXA_INCHIKEYBUILDER_HANDLE = *mut IXA_INCHIKEYBUILDER_HANDLE_STRUCT;
pub type IXA_ATOMID                 = *mut IXA_ATOMID_STRUCT;
pub type IXA_BONDID                 = *mut IXA_BONDID_STRUCT;
pub type IXA_STEREOID               = *mut IXA_STEREOID_STRUCT;
pub type IXA_POLYMERUNITID          = *mut IXA_POLYMERUNITID_STRUCT;

// IXA_BOOL and IXA_STATUS enum values
pub type IXA_BOOL = ::std::os::raw::c_int;
pub const IXA_FALSE: IXA_BOOL = 0;
pub const IXA_TRUE:  IXA_BOOL = 1;

pub type IXA_STATUS = ::std::os::raw::c_int;
pub const IXA_STATUS_SUCCESS: IXA_STATUS = 0;
pub const IXA_STATUS_WARNING: IXA_STATUS = 1;
pub const IXA_STATUS_ERROR:   IXA_STATUS = 2;
```

**New extern functions to add (IXA, ~50 most-useful ones):**

Status API, Mol read/query API, InChI builder API, InChIKey builder API. The full IXA building API (CreateAtom, SetAtomElement, …) is skipped per §3.2.

**INCHIGEN types and functions:**

```rust
pub type INCHIGEN_HANDLE = *mut ::std::os::raw::c_void;

// INCHIGEN_DATA — needed as an out-param for Setup/Do* calls
// Contains normalized atom data; treat as opaque from Rust (use *mut c_void)
// Only INCHIGEN_DoSerialization produces user-visible output (inchi_Output).
```

The `NORM_ATOMS` types inside `INCHIGEN_DATA` are not bound (see §3.3). `INCHIGEN_DATA` is passed as a raw pointer to a stack-allocated zeroed blob of known size; the C library writes into it. Size must be determined from `sizeof(INCHIGEN_DATA)` at build time via `build.rs` or a fixed constant derived from the header.

---

## 5. Phase 2 — `IxaMolecule`: parsing and querying

**New file:** `crates/inchi/src/ixa.rs`

### 5.1 What this adds

The IXA API provides two capabilities the classic API lacks:

1. **Round-trip from InChI to a queryable molecule.** `struct_from_inchi` already reconstructs atoms and bonds but returns a `Structure` with fixed arrays. `IxaMolecule::from_inchi` uses `IXA_MOL_ReadInChI` to populate a live mol object whose atoms and bonds can then be read back individually via the IXA getter API — including element symbol, atomic number, mass, charge, radical, H counts, and 3D coordinates.

2. **Parsing a Molfile into a queryable molecule.** `from_molfile` today returns only the InChI. `IxaMolecule::from_molfile` uses `IXA_MOL_ReadMolfile` to populate a live mol object, enabling atom/bond inspection without regenerating the InChI.

### 5.2 ID newtypes

IXA atom/bond/stereo IDs are raw C pointers cast from runtime integers. In Rust they are wrapped in `Copy` newtypes that carry no lifetime — the caller is responsible for using them only with the `IxaMolecule` that produced them. Mixing IDs from different molecules is undefined behavior at the C level; the Rust API makes this unlikely through naming conventions but cannot enforce it statically without `PhantomData` lifetimes (which would make the iterator APIs awkward).

```rust
/// An atom identifier returned by [`IxaMolecule`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomId(inchi_sys::IXA_ATOMID);

/// A bond identifier returned by [`IxaMolecule`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BondId(inchi_sys::IXA_BONDID);

/// A stereo descriptor identifier returned by [`IxaMolecule`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StereoId(inchi_sys::IXA_STEREOID);
```

### 5.3 Internal layout

```rust
// crate-internal; never exposed publicly
struct IxaHandles {
    status: inchi_sys::IXA_STATUS_HANDLE,
    mol:    inchi_sys::IXA_MOL_HANDLE,
}

impl Drop for IxaHandles {
    fn drop(&mut self) {
        let _guard = raw::lock();
        // SAFETY: handles are valid (non-null); status was created before mol.
        unsafe {
            inchi_sys::IXA_MOL_Destroy(self.status, self.mol);
            inchi_sys::IXA_STATUS_Destroy(self.status);
        }
    }
}

/// A molecule loaded from a Molfile or InChI string, queryable atom-by-atom.
pub struct IxaMolecule(IxaHandles);
```

`IxaHandles` is not `Send` or `Sync` by default (raw pointers). We add:

```rust
// SAFETY: IxaHandles contains opaque C heap pointers that are valid for the
// lifetime of the struct and accessed only under INCHI_LOCK.
unsafe impl Send for IxaHandles {}
unsafe impl Sync for IxaHandles {}
```

### 5.4 Constructor methods

```rust
impl IxaMolecule {
    /// Parses a Molfile (V2000/V3000) text string into a queryable molecule.
    pub fn from_molfile(molfile: impl AsRef<str>) -> Result<Self>;

    /// Reconstructs a molecule from an InChI string.
    ///
    /// Coordinates are not stored in InChI; all atom positions are zero.
    pub fn from_inchi(inchi: impl AsRef<str>) -> Result<Self>;
}
```

Both constructors:
1. Acquire `INCHI_LOCK`
2. Call `IXA_STATUS_Create()` — fails only on OOM, return `Err`
3. Call `IXA_MOL_Create(status)` — same
4. Call `IXA_MOL_ReadMolfile` / `IXA_MOL_ReadInChI`
5. Call `IXA_STATUS_HasError` — if true, collect messages and return `Err(InchiError::Failed)`
6. Return `Ok(IxaMolecule(IxaHandles { status, mol }))`

### 5.5 Query methods

Each getter acquires the lock, calls the C function, clears the status, and returns the result or an error.

```rust
impl IxaMolecule {
    /// Number of heavy atoms (excludes implicit hydrogens).
    pub fn atom_count(&self) -> usize;

    /// Number of bonds.
    pub fn bond_count(&self) -> usize;

    /// Number of stereo descriptors.
    pub fn stereo_count(&self) -> usize;

    /// Returns an iterator over all atom IDs.
    pub fn atoms(&self) -> impl Iterator<Item = AtomId> + '_;

    /// Returns an iterator over all bond IDs.
    pub fn bonds(&self) -> impl Iterator<Item = BondId> + '_;

    /// Element symbol (e.g. `"C"`, `"Fe"`).
    pub fn atom_element(&self, atom: AtomId) -> Result<String>;

    /// Atomic number.
    pub fn atom_atomic_number(&self, atom: AtomId) -> Result<u32>;

    /// Mass number (0 = natural abundance).
    pub fn atom_mass(&self, atom: AtomId) -> Result<i32>;

    /// Formal charge.
    pub fn atom_charge(&self, atom: AtomId) -> Result<i32>;

    /// Radical state.
    pub fn atom_radical(&self, atom: AtomId) -> Result<Radical>;

    /// Explicit hydrogen counts at each isotopic position [H, 1H, 2H, 3H, 4H].
    pub fn atom_hydrogens(&self, atom: AtomId) -> Result<[i32; 5]>;

    /// Cartesian coordinates (all zero for InChI-derived molecules).
    pub fn atom_position(&self, atom: AtomId) -> Result<(f64, f64, f64)>;

    /// Number of bonds attached to this atom.
    pub fn atom_bond_count(&self, atom: AtomId) -> Result<usize>;

    /// The bond connecting `atom` to its `n`th neighbor.
    pub fn atom_bond(&self, atom: AtomId, n: usize) -> Result<BondId>;

    /// Atoms at both ends of a bond.
    pub fn bond_atoms(&self, bond: BondId) -> Result<(AtomId, AtomId)>;

    /// Bond order.
    pub fn bond_order(&self, bond: BondId) -> Result<BondOrder>;

    /// Whether the molecule has been flagged as chiral.
    pub fn is_chiral(&self) -> bool;
}
```

### 5.6 InChI generation from `IxaMolecule`

Uses `IXA_INCHIBUILDER_*` internally. The builder handle is created, used, and destroyed within a single `to_inchi` call (not stored), so no additional RAII type is needed.

```rust
impl IxaMolecule {
    /// Generates an InChI from the current molecule state.
    pub fn to_inchi(&self, options: impl Into<Options>) -> Result<InchiOutput>;

    /// Generates an InChIKey from the current molecule state (feature `key`).
    #[cfg(feature = "key")]
    pub fn to_inchikey(&self, options: impl Into<Options>) -> Result<String>;
}
```

Implementation of `to_inchi`:
1. Lock
2. `IXA_INCHIBUILDER_Create(self.status, self.mol)`
3. For each option: `IXA_INCHIBUILDER_SetOption` / `IXA_INCHIBUILDER_SetOption_Stereo`
4. `IXA_INCHIBUILDER_GetInChI` → read the returned `const char*` into `String`
5. `IXA_INCHIBUILDER_GetAuxInfo` → same
6. `IXA_INCHIBUILDER_Destroy`
7. Check status, return `InchiOutput`

Note: the `const char*` returned by `IXA_INCHIBUILDER_GetInChI` is owned by the builder handle and becomes invalid after `Destroy`. Rust must copy it into a `String` before step 6.

### 5.7 Public re-exports

In `lib.rs`:
```rust
pub use ixa::{AtomId, BondId, IxaMolecule, StereoId};
```

Add to API overview:
```
- [`IxaMolecule::from_molfile`] — parse a Molfile into a queryable mol object.
- [`IxaMolecule::from_inchi`] — reconstruct atom/bond structure from an InChI string.
```

---

## 6. Phase 3 — `InchiGenerator`: step-through generation

**Priority:** Low. Implement only after Phase 2. Gate behind a Cargo feature flag `generator` (disabled by default).

### 6.1 Value proposition

`InchiGenerator` exposes the four internal steps of InChI generation — setup, normalization, canonicalization, serialization — as separate Rust calls. This is useful for:
- Diagnosing which step fails for a pathological molecule
- Inspecting the error message produced at each step
- Research use cases that need the normalized structure before serialization

It does **not** expose `INCHIGEN_DATA` internals (see §3.3); only the error/warning string at each step.

### 6.2 Design

```rust
/// A step-through InChI generator.
///
/// Wraps the `INCHIGEN_*` C API. Each step returns the accumulated warning/error
/// message from that stage. Obtain the final InChI by calling [`do_serialization`].
///
/// This type is only useful when you need per-step diagnostic information.
/// For normal generation use [`Molecule::to_inchi`].
#[cfg(feature = "generator")]
pub struct InchiGenerator {
    handle: *mut ::std::os::raw::c_void, // INCHIGEN_HANDLE (void*)
}

#[cfg(feature = "generator")]
impl InchiGenerator {
    pub fn new() -> Result<Self>;

    /// Converts `molecule` into the internal C input format and registers it.
    /// Must be called before normalization.
    pub fn setup(&mut self, molecule: &Molecule, options: impl Into<Options>) -> Result<StepReport>;

    /// Runs the normalization step.
    pub fn do_normalization(&mut self) -> Result<StepReport>;

    /// Runs the canonicalization step.
    pub fn do_canonicalization(&mut self) -> Result<StepReport>;

    /// Runs the serialization step and returns the final InChI.
    pub fn do_serialization(&mut self) -> Result<InchiOutput>;

    /// Resets internal state; call before reusing the generator for another molecule.
    pub fn reset(&mut self);
}

/// Diagnostic message produced by a single InChI generation step.
#[cfg(feature = "generator")]
#[non_exhaustive]
pub struct StepReport {
    pub message: String,
    pub status: Status,
}
```

`INCHIGEN_DATA` is allocated on the Rust stack as a zero-initialized byte array of `size_of::<INCHIGEN_DATA>()` and passed by mutable raw pointer to each `Do*` call. Its internal `NORM_ATOMS` fields are never read from Rust. This avoids binding `NORM_ATOMS` while still satisfying the C API's requirement for a writable `INCHIGEN_DATA*`.

The size of `INCHIGEN_DATA` must be known at compile time. Two options:
- **Option A:** Add a `build.rs` step that uses `cc` to compile a small C probe and emit `cargo:rustc-env=INCHIGEN_DATA_SIZE=N`.
- **Option B:** Derive the size from the header manually (it contains 4 pointers + one fixed-size char array = `STR_ERR_LEN + 4 * sizeof(pointer)`). Compute the constant in `bindings.rs` and use `[u8; INCHIGEN_DATA_SIZE]` as the Rust representation.

Option B is simpler and avoids an extra build step; prefer it.

`INCHIGEN_DATA` between steps is owned by `InchiGenerator` as a `Box<[u8; INCHIGEN_DATA_SIZE]>` field. This keeps it alive and pinned for the duration of the generator's lifetime.

---

## 7. Phase 4 — IXA mutation (optional, future)

If a use case emerges that requires modifying a molecule after loading it from a Molfile (e.g., change an atom's charge, add a bond), the IXA setter API (`IXA_MOL_SetAtomCharge`, `IXA_MOL_SetAtomMass`, etc.) can be exposed as `&mut self` methods on `IxaMolecule`:

```rust
impl IxaMolecule {
    pub fn set_atom_charge(&mut self, atom: AtomId, charge: i32) -> Result<()>;
    pub fn set_atom_mass(&mut self, atom: AtomId, mass: i32) -> Result<()>;
    // ...
}
```

This phase is intentionally deferred. The existing `Molecule` builder already handles construction from scratch; mutation after parsing is a narrower need. File an issue with a concrete use case before starting this work.

---

## 8. Rollout order and semver

| Phase | Feature gate | SemVer impact |
|---|---|---|
| Phase 1 — IXA bindings in `inchi-sys` | none | patch (`inchi-sys` is a `-sys` crate; adding extern declarations is non-breaking) |
| Phase 2 — `IxaMolecule` | none (always on) | minor (new public API) |
| Phase 3 — `InchiGenerator` | `generator` (opt-in) | minor |
| Phase 4 — IXA mutation | part of `IxaMolecule` | minor |

All phases add new API; none break existing callers.

---

## 9. Testing strategy

**Phase 1:** No new tests needed; binding declarations are validated by the build.

**Phase 2:**
- `IxaMolecule::from_molfile` round-trip: load each fixture `.mol`, call `to_inchi()`, assert against the expected InChI from `vectors.rs`.
- `IxaMolecule::from_inchi`: parse each expected InChI from `vectors.rs`, verify `atom_count()` matches the formula, verify `atom_element(id)` returns expected elements.
- Atom query: load `caffeine.mol` (8 heavy atoms + 10N), assert `atom_count() == 14`, spot-check element symbols and charges.
- Error propagation: pass a garbage string to `from_inchi`, assert `Err`.

**Phase 3:**
- Step-through on each fixture molecule: assert each step returns `Status::Success` and `do_serialization()` matches `Molecule::to_inchi()`.
- Partial step: call `setup` then `do_normalization` only, assert no panic or UB.

---

## 10. Files to create/modify

| File | Change |
|---|---|
| `crates/inchi-sys/src/bindings.rs` | Add IXA type defs + extern declarations; add INCHIGEN types + extern declarations |
| `crates/inchi/src/ixa.rs` | New file — `IxaMolecule`, `AtomId`, `BondId`, `StereoId` |
| `crates/inchi/src/lib.rs` | `mod ixa; pub use ixa::{...};` + API overview update |
| `crates/inchi/Cargo.toml` | Add `generator` feature flag (Phase 3) |
| `crates/inchi/src/generator.rs` | New file — `InchiGenerator`, `StepReport` (Phase 3, feature-gated) |
| `crates/inchi/tests/vectors.rs` | Add `IxaMolecule` round-trip tests (Phase 2) |
| `CLAUDE.md` | Update module list (add `ixa.rs`, `generator.rs`) |
