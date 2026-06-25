//! The [`Options`] builder for controlling InChI generation.

use std::time::Duration;

/// How stereochemistry is interpreted during InChI generation.
///
/// Any value other than [`StereoMode::Standard`] (except [`StereoMode::Ignore`],
/// which is standard-compatible) causes a **non-standard** InChI to be
/// produced.
///
/// ```
/// use inchi::{Options, StereoMode};
///
/// let opts = Options::new().stereo(StereoMode::Ignore);
/// assert!(opts.is_standard()); // /SNon keeps the InChI standard
///
/// let rel = Options::new().stereo(StereoMode::Relative);
/// assert!(!rel.is_standard());
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum StereoMode {
    /// Default stereo perception, producing a standard InChI.
    #[default]
    Standard,
    /// Ignore all stereochemistry (`/SNon`). Standard-compatible.
    Ignore,
    /// Relative stereo (`/SRel`). Produces a non-standard InChI.
    Relative,
    /// Racemic stereo (`/SRac`). Produces a non-standard InChI.
    Racemic,
    /// Use the chiral flag from the input (`/SUCF`). Produces a non-standard InChI.
    UseChiralFlag,
}

/// A builder for InChI generation options.
///
/// `Options` produces the space-delimited option string the InChI library
/// expects, with the correct per-platform flag prefix (`-` on Unix, `/` on
/// Windows), so you never construct raw flag strings by hand.
///
/// The default ([`Options::new`]) requests a **standard** InChI with auxiliary
/// information included.
///
/// ```
/// use inchi::Options;
///
/// let opts = Options::new()
///     .fixed_h(true)
///     .reconnect_metals(true)
///     .aux_info(false);
/// assert!(!opts.is_standard()); // FixedH / RecMet are non-standard
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Options {
    add_hydrogens: bool,
    aux_info: bool,
    stereo: StereoMode,
    fixed_h: bool,
    reconnect_metals: bool,
    keto_enol: bool,
    tautomer_15: bool,
    chiral_flag: Option<bool>,
    timeout: Option<Duration>,
    save_opt: bool,
    polymers: Polymers,
    no_frame_shift: bool,
    fold_sru: bool,
    no_edits: bool,
    allow_pseudo_atoms: bool,
    extra: Vec<String>,
}

/// How polymer structural-repeating-unit (SRU) data is processed.
///
/// Polymer support is an experimental, **non-standard** extension (the
/// resulting InChI carries a `B` "beta" version flag and a `/z` layer). It is
/// only meaningful for Molfiles that carry CTFile Sgroup (`STY SRU`) polymer
/// data, parsed via [`from_molfile`](crate::from_molfile).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum Polymers {
    /// Do not process polymer data (the default).
    #[default]
    Off,
    /// Process polymers with the current algorithm (`/Polymers`).
    On,
    /// Process polymers with the legacy v1.05 algorithm (`/Polymers105`).
    Legacy,
}

impl Default for Options {
    fn default() -> Self {
        Self::new()
    }
}

/// `()` converts to the default options, so the common "no special options"
/// case can be written as `from_molfile(mol, ())` instead of
/// `from_molfile(mol, &Options::new())`.
///
/// ```
/// use inchi::{from_molfile, Options};
/// let m = "\n  ex\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n\
///     \x20   0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
/// assert_eq!(from_molfile(m, ())?.inchi(), from_molfile(m, Options::new())?.inchi());
/// # Ok::<(), inchi::InchiError>(())
/// ```
impl From<()> for Options {
    fn from((): ()) -> Self {
        Options::new()
    }
}

/// Borrowed options convert by cloning, so a `&Options` still works wherever an
/// owned `impl Into<Options>` is expected (e.g. when reusing one option set).
impl From<&Options> for Options {
    fn from(options: &Options) -> Self {
        options.clone()
    }
}

impl Options {
    /// Creates a new option set requesting a standard InChI with auxiliary
    /// information.
    ///
    /// ```
    /// use inchi::Options;
    /// assert!(Options::new().is_standard());
    /// ```
    #[must_use]
    pub fn new() -> Self {
        Options {
            add_hydrogens: true,
            aux_info: true,
            stereo: StereoMode::Standard,
            fixed_h: false,
            reconnect_metals: false,
            keto_enol: false,
            tautomer_15: false,
            chiral_flag: None,
            timeout: None,
            save_opt: false,
            polymers: Polymers::Off,
            no_frame_shift: false,
            fold_sru: false,
            no_edits: false,
            allow_pseudo_atoms: false,
            extra: Vec::new(),
        }
    }

    /// Controls whether the library adds implicit hydrogens (default `true`).
    ///
    /// Setting this to `false` emits `/DoNotAddH`.
    ///
    /// ```
    /// use inchi::Options;
    /// let opts = Options::new().add_hydrogens(false);
    /// assert!(opts.is_standard());
    /// ```
    #[must_use]
    pub fn add_hydrogens(mut self, yes: bool) -> Self {
        self.add_hydrogens = yes;
        self
    }

    /// Controls whether `AuxInfo` is produced (default `true`).
    ///
    /// Setting this to `false` emits `/AuxNone` and leaves
    /// [`InchiOutput::aux_info`](crate::InchiOutput::aux_info) empty.
    ///
    /// ```
    /// use inchi::Options;
    /// let _ = Options::new().aux_info(false);
    /// ```
    #[must_use]
    pub fn aux_info(mut self, yes: bool) -> Self {
        self.aux_info = yes;
        self
    }

    /// Sets how stereochemistry is interpreted.
    ///
    /// ```
    /// use inchi::{Options, StereoMode};
    /// let _ = Options::new().stereo(StereoMode::Racemic);
    /// ```
    #[must_use]
    pub fn stereo(mut self, mode: StereoMode) -> Self {
        self.stereo = mode;
        self
    }

    /// Includes the fixed-hydrogen layer (`/FixedH`). Non-standard.
    ///
    /// ```
    /// use inchi::Options;
    /// assert!(!Options::new().fixed_h(true).is_standard());
    /// ```
    #[must_use]
    pub fn fixed_h(mut self, yes: bool) -> Self {
        self.fixed_h = yes;
        self
    }

    /// Reconnects bonds to metal atoms (`/RecMet`). Non-standard.
    ///
    /// ```
    /// use inchi::Options;
    /// assert!(!Options::new().reconnect_metals(true).is_standard());
    /// ```
    #[must_use]
    pub fn reconnect_metals(mut self, yes: bool) -> Self {
        self.reconnect_metals = yes;
        self
    }

    /// Accounts for keto-enol tautomerism (`/KET`). Non-standard.
    ///
    /// ```
    /// use inchi::Options;
    /// let _ = Options::new().keto_enol_tautomerism(true);
    /// ```
    #[must_use]
    pub fn keto_enol_tautomerism(mut self, yes: bool) -> Self {
        self.keto_enol = yes;
        self
    }

    /// Accounts for 1,5-tautomerism (`/15T`). Non-standard.
    ///
    /// ```
    /// use inchi::Options;
    /// let _ = Options::new().tautomerism_15(true);
    /// ```
    #[must_use]
    pub fn tautomerism_15(mut self, yes: bool) -> Self {
        self.tautomer_15 = yes;
        self
    }

    /// Forces the chiral flag on (`/ChiralFlagON`) or off (`/ChiralFlagOFF`),
    /// or leaves it unset (`None`). Non-standard when set.
    ///
    /// ```
    /// use inchi::Options;
    /// let _ = Options::new().chiral_flag(Some(true));
    /// ```
    #[must_use]
    pub fn chiral_flag(mut self, flag: Option<bool>) -> Self {
        self.chiral_flag = flag;
        self
    }

    /// Sets a per-structure timeout (`/WMnumber`, milliseconds).
    ///
    /// A zero duration means "unlimited", matching the library default.
    ///
    /// ```
    /// use inchi::Options;
    /// use std::time::Duration;
    /// let _ = Options::new().timeout(Duration::from_secs(30));
    /// ```
    #[must_use]
    pub fn timeout(mut self, dur: Duration) -> Self {
        self.timeout = Some(dur);
        self
    }

    /// Saves custom InChI creation options into the output (`/SaveOpt`).
    /// Non-standard.
    ///
    /// ```
    /// use inchi::Options;
    /// let _ = Options::new().save_opt(true);
    /// ```
    #[must_use]
    pub fn save_opt(mut self, yes: bool) -> Self {
        self.save_opt = yes;
        self
    }

    /// Enables polymer processing (`/Polymers` or `/Polymers105`). Non-standard.
    ///
    /// Use this with [`from_molfile`](crate::from_molfile) on a Molfile carrying
    /// CTFile polymer Sgroups; the result is a beta-flagged, non-standard InChI
    /// with a `/z` polymer layer.
    ///
    /// ```
    /// use inchi::{Options, Polymers};
    /// let opts = Options::new().polymers(Polymers::On);
    /// assert!(!opts.is_standard());
    /// ```
    #[must_use]
    pub fn polymers(mut self, mode: Polymers) -> Self {
        self.polymers = mode;
        self
    }

    /// Disables the polymer CRU frame shift (`/NoFrameShift`).
    ///
    /// Only meaningful together with [`Options::polymers`].
    ///
    /// ```
    /// use inchi::{Options, Polymers};
    /// let _ = Options::new().polymers(Polymers::On).no_frame_shift(true);
    /// ```
    #[must_use]
    pub fn no_frame_shift(mut self, yes: bool) -> Self {
        self.no_frame_shift = yes;
        self
    }

    /// Folds the polymer SRU where possible (`/FoldSRU`).
    ///
    /// Only meaningful together with [`Options::polymers`].
    ///
    /// ```
    /// use inchi::{Options, Polymers};
    /// let _ = Options::new().polymers(Polymers::On).fold_sru(true);
    /// ```
    #[must_use]
    pub fn fold_sru(mut self, yes: bool) -> Self {
        self.fold_sru = yes;
        self
    }

    /// Disables polymer CRU frame shift *and* folding (`/NoEdits`).
    ///
    /// Only meaningful together with [`Options::polymers`].
    ///
    /// ```
    /// use inchi::{Options, Polymers};
    /// let _ = Options::new().polymers(Polymers::On).no_edits(true);
    /// ```
    #[must_use]
    pub fn no_edits(mut self, yes: bool) -> Self {
        self.no_edits = yes;
        self
    }

    /// Allows non-polymer `Zz` pseudo-element placeholder atoms (`/NPZz`).
    ///
    /// ```
    /// use inchi::Options;
    /// let _ = Options::new().allow_pseudo_atoms(true);
    /// ```
    #[must_use]
    pub fn allow_pseudo_atoms(mut self, yes: bool) -> Self {
        self.allow_pseudo_atoms = yes;
        self
    }

    /// Appends a raw option token (without its `/` or `-` prefix) as an escape
    /// hatch for options not modeled by a dedicated method.
    ///
    /// The token is passed through verbatim, so consult the upstream InChI
    /// technical manual for valid values.
    ///
    /// ```
    /// use inchi::Options;
    /// let _ = Options::new().raw("NEWPSOFF");
    /// ```
    #[must_use]
    pub fn raw(mut self, token: impl Into<String>) -> Self {
        self.extra.push(token.into());
        self
    }

    /// Returns `true` if these options produce a *standard* InChI.
    ///
    /// ```
    /// use inchi::{Options, StereoMode};
    /// assert!(Options::new().is_standard());
    /// assert!(!Options::new().fixed_h(true).is_standard());
    /// ```
    #[must_use]
    pub fn is_standard(&self) -> bool {
        !self.fixed_h
            && !self.reconnect_metals
            && !self.keto_enol
            && !self.tautomer_15
            && self.chiral_flag.is_none()
            && matches!(self.stereo, StereoMode::Standard | StereoMode::Ignore)
            && matches!(self.polymers, Polymers::Off)
            && self.extra.is_empty()
    }

    /// The platform-appropriate option prefix: `/` on Windows, `-` elsewhere.
    const PREFIX: &'static str = if cfg!(windows) { "/" } else { "-" };

    /// Renders the option string passed to the InChI C library.
    ///
    /// ```
    /// use inchi::Options;
    /// // Default options request only auxiliary info, which is the library default,
    /// // so the rendered string is empty.
    /// assert_eq!(Options::new().to_arg_string(), "");
    /// ```
    #[must_use]
    pub fn to_arg_string(&self) -> String {
        let mut tokens: Vec<String> = Vec::new();
        let mut push = |t: &str| tokens.push(format!("{}{t}", Self::PREFIX));

        if !self.add_hydrogens {
            push("DoNotAddH");
        }
        if !self.aux_info {
            push("AuxNone");
        }
        match self.stereo {
            StereoMode::Standard => {}
            StereoMode::Ignore => push("SNon"),
            StereoMode::Relative => push("SRel"),
            StereoMode::Racemic => push("SRac"),
            StereoMode::UseChiralFlag => push("SUCF"),
        }
        if self.fixed_h {
            push("FixedH");
        }
        if self.reconnect_metals {
            push("RecMet");
        }
        if self.keto_enol {
            push("KET");
        }
        if self.tautomer_15 {
            push("15T");
        }
        match self.chiral_flag {
            Some(true) => push("ChiralFlagON"),
            Some(false) => push("ChiralFlagOFF"),
            None => {}
        }
        if let Some(dur) = self.timeout {
            let ms = dur.as_millis();
            push(&format!("WM{ms}"));
        }
        if self.save_opt {
            push("SaveOpt");
        }
        match self.polymers {
            Polymers::Off => {}
            Polymers::On => push("Polymers"),
            Polymers::Legacy => push("Polymers105"),
        }
        if self.no_frame_shift {
            push("NoFrameShift");
        }
        if self.fold_sru {
            push("FoldSRU");
        }
        if self.no_edits {
            push("NoEdits");
        }
        if self.allow_pseudo_atoms {
            push("NPZz");
        }
        for token in &self.extra {
            push(token);
        }

        tokens.join(" ")
    }
}
