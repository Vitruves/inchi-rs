//! Internal helpers bridging the safe API and the `inchi-sys` FFI surface.
//!
//! Everything here is crate-private. This is the *only* module that touches
//! raw pointers or `unsafe`; the rest of the crate stays safe.

use crate::error::InchiError;
use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::sync::{Mutex, MutexGuard};

/// Global lock serializing all calls into the InChI C library.
///
/// The reference library keeps `static` internal state (e.g. a memoized
/// hydrogen element index) and is not guaranteed to be re-entrant, so the safe
/// API funnels every C call through this mutex.
static INCHI_LOCK: Mutex<()> = Mutex::new(());

/// Acquires the global InChI lock, recovering from poisoning.
///
/// Poisoning can only occur if a thread panicked while holding the guard. The
/// guarded data is `()`, so there is no invariant to protect; we simply reuse
/// the lock rather than propagate the poison (which would surface as a panic).
pub(crate) fn lock() -> MutexGuard<'static, ()> {
    match INCHI_LOCK.lock() {
        Ok(guard) => guard,
        Err(poisoned) => poisoned.into_inner(),
    }
}

/// Converts a Rust string into a NUL-terminated C string, mapping interior
/// NULs to a typed error instead of panicking.
pub(crate) fn to_cstring(s: &str) -> Result<CString, InchiError> {
    CString::new(s).map_err(|e| InchiError::InteriorNul {
        position: e.nul_position(),
    })
}

/// Reads a C string pointer into an owned `String`, treating null as empty and
/// replacing any non-UTF-8 bytes (InChI output is ASCII, so this is lossless in
/// practice). Safe to call with a null pointer.
///
/// # Safety
///
/// `ptr` must be either null or a valid pointer to a NUL-terminated C string
/// that remains valid for the duration of the call.
pub(crate) unsafe fn cstr_to_string(ptr: *const c_char) -> String {
    if ptr.is_null() {
        String::new()
    } else {
        CStr::from_ptr(ptr).to_string_lossy().into_owned()
    }
}

/// RAII wrapper around an [`inchi_sys::inchi_Output`] that guarantees
/// `FreeINCHI` runs exactly once, even on early return or panic.
pub(crate) struct OutputGuard {
    raw: inchi_sys::inchi_Output,
}

impl OutputGuard {
    /// Creates a zeroed output ready to be populated by a `Get*`/`Make*` call.
    pub(crate) fn new() -> Self {
        // `inchi_Output` is a plain-old-data struct of pointers/lengths; an
        // all-zero value is the documented "empty" state.
        OutputGuard {
            raw: unsafe { std::mem::zeroed() },
        }
    }

    /// Returns a mutable pointer to the underlying C struct for FFI calls.
    pub(crate) fn as_mut_ptr(&mut self) -> *mut inchi_sys::inchi_Output {
        &mut self.raw
    }

    /// The generated InChI string (empty if none was produced).
    pub(crate) fn inchi(&self) -> String {
        unsafe { cstr_to_string(self.raw.szInChI) }
    }

    /// The auxiliary-information string (empty if absent).
    pub(crate) fn aux_info(&self) -> String {
        unsafe { cstr_to_string(self.raw.szAuxInfo) }
    }

    /// The error/warning message (empty if none).
    pub(crate) fn message(&self) -> String {
        unsafe { cstr_to_string(self.raw.szMessage) }
    }

    /// The human-readable log string (empty if none).
    pub(crate) fn log(&self) -> String {
        unsafe { cstr_to_string(self.raw.szLog) }
    }
}

impl Drop for OutputGuard {
    fn drop(&mut self) {
        // `FreeINCHI` tolerates already-null fields and frees only what the C
        // side allocated; safe to call unconditionally.
        unsafe { inchi_sys::FreeINCHI(&mut self.raw) }
    }
}
