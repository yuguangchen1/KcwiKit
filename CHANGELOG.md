# Changelog

All notable changes to KCWIKit since the 0.2 release (September 2024). Entries
name the commit on the `dev` branch and the author in brackets: YC = Yuguang
Chen, YD = Yuanze Ding, ZZ = Zhuyun Zhuang.

## Unreleased (dev, November 2026) — review fixes

Fixes from a code review of everything added since October 2024
(`0e01346`, YD). Behaviour changes are marked **(behaviour)**.

### Fixed
- **`kcwi_stack` regression (behaviour).** Commit `95de13b` (Sep 2025) had
  replaced `kcwi_stack` with a pre-merge copy and silently dropped Yuguang's
  Dec 2024–Jan 2025 fixes. Restored: per-frame absolute-astrometry offsets
  from a 3-column `.astrom.list` (mosaics; a 2-column file still applies one
  offset to all frames, and a frame without an entry is now an error instead of
  a silent mis-match); the corrected `BUNIT` strings
  (`1e-16 erg/s/cm2/Angstrom/arcsec2` and `1e-32 erg2/s2/cm4/Angstrom2/arcsec4`
  instead of the dimensionally wrong `cm3` / `cm6` forms); the Montage
  issue #80 work-around (`ISS6FIX` integer scaling), now applied **only on the
  MontagePy path** since the command-line `mProjectCube` is not affected and
  the int32 quantisation costs precision on faint pixels; and removal of the
  `independent_celestial_slices` keyword, which current `reproject` (≥ 0.10)
  no longer accepts.
- **`pyregion` → `regions` migration completed.** `setup.py` had dropped
  `pyregion` in Dec 2024 while `kcwi.py`, `kcwi_tools.py`, `kcwi_mask_cube.py`
  and `kcwi_xalign.py` still imported it, so a fresh install could not
  `import kcwikit.kcwi`. All region handling now goes through the new
  `kcwikit/kcwi/kcwi_regions.py` (`regions_to_mask`), which restores the three
  behaviours lost in the first migration: region files in sky coordinates
  (fk5/icrs) are converted with the image WCS, `-` exclusion shapes are
  subtracted, and a shape entirely outside the image contributes nothing
  instead of raising. It uses `regions` and falls back to `pyregion` if only
  that is installed. `kcwi_pyDRP.py` / `kcwi_pypeit.py` (legacy) import
  `pyregion` optionally.
- `kcwi_masksky_ds9_zap`: an empty region file masked the whole image
  (`allmask[None] = 2`); the other five `masksky` scripts crashed. All six now
  write an empty mask with a warning.
- `kcwi_stack`: `AIRMASS` in the stacked header was `0.0` when trimmed cubes
  from a previous run were reused; it is now read from each trimmed cube
  after the exposure-time correction.
- `kcwi_vachelio` **(behaviour)**: the `UNCERT` extension was resampled with
  nearest-neighbour like a mask while the science cube is cubic; it is now
  linear in variance (σ², then √), so σ and flux stay aligned. MASK/FLAGS
  unchanged (nearest-up), science unchanged (cubic).
- `empty_hdul_like`: handles table and header-only HDUs instead of setting
  `header = None`.
- `kcwi_astrometry`: the reference image is read from extension `nhdu` when it
  holds a 2-D image, otherwise from the first extension that does (MEF and
  compressed files), restoring Yuguang's `427e63c` fallback.
- `crr_grow_rel` default is `0.1` everywhere (`kcwi_stack`,
  `_crr_flag_outliers`, `kcwi_combine_stacks`); the `kcwi_stack` docstring now
  says `crr_split_mode='auto'` is the default.
- `kcrm_create_crmsk`: `--sigma` was accepted but ignored (both clips used 3);
  it is now honoured. New `--ratio` (relative-excess threshold, default 0.6,
  was hard-coded) and `--proc` (processing table, default `kcwir.proc`); the
  table is read once instead of once per frame; division by zero in gaps is
  guarded; the debugging `print` is replaced by a per-frame summary.
- `kcwi_gen_skyfile`: rewritten; no longer writes `tmp.pickle` into the working
  directory on every run, `-l` is required, the unused `-f` is gone, the table
  is written directly in the whitespace format `MakeMasterSky` reads, and the
  docstring/`-h` text (copied from `kcwi_gen_log`) is corrected.
- `kcwi_mask_cube`: rewritten; `pyregion` dropped, `-r`/`-b` required (was a
  `NameError` without them), `FLAGS` bit 32 is OR-ed in (was `+= 32`, which
  corrupted the bitmask on a second run), no `KeyError` on cubes without
  `CNAME3`, and the grating-specific white-light windows are gone (only the
  spatial WCS is needed).
- `kcwi_shiftcheck`: imports `kcwi_xalign` from the package instead of a
  `sys.path` hack (works when installed, not only from the source tree); the
  blue-vs-red summary reads the stack pixel scale from the `.par` file instead
  of assuming `0.146"`.
- `kcwi_skylevel`: the process-wide `warnings.simplefilter('ignore')` now
  silences only numpy nan-reduction messages, so astropy WCS/FITS warnings
  show again; empty wavelength windows exit with a clear message or skip the
  anchor; the dichroic step checks that the two stacks share shape **and**
  WCS; image and cube-window anchors use their own channel's WCS.
- `kcwi_xalign`: an `_rm.reg` that exists but cannot be applied is now an
  error rather than a printed warning that silently changed the alignment.
- `kcwi_automask`: `--grow`, `--maxfrac`, `--shape`, `--regfrom` help text
  matches the code; region files are read through the shared helper (sky
  coordinates, exclusions) with the built-in parser as fallback.
- Documentation: `Running_red_DRP.md` (said blue / `-b kb*`; step numbering;
  `kcwib.proc` vs `kcwir.proc`; duplicate sentence; recovery-region
  instructions replaced by the `--bcrmsk` / `-s` workflow),
  `Running_blue_DRP.md` (where `bcrmsk` goes in `kcwi.cfg`, `kcwi_automask
  --bcrmsk` workflow), `kcwi_masksky_ds9_bcrr` docstring.

### Added
- `kcwikit/kcwi/kcwi_regions.py` — shared DS9-region → mask helper (see above).
- `setup.py` entry points for `kcwi_shiftcheck` and `kcwi_mask_cube`.

### Removed
- Committed backup copies `kcwi_automask.py.bak_wlpoly`,
  `kcwi_xalign.py.bak_astrom`.

### Known remaining items
- `kcwi.py.bak_crmask`, `kcwi_deepcr_crmsk.py` (never functional) and
  `scripts/tmp_telluric.ipynb` are still tracked and slated for removal.
- The `pyDRP/` local copy of DRP primitives has diverged from the
  `KCWI_DRP` fork (`MakeCube` warp order, `FitCenter` keyword,
  `RemoveCosmicRays` file names).

## 2026-11-05 — sky-level matching across the dichroic (`0a2302c`, YD)

### Added
- `kcwi_skylevel` (`measure` / `apply` / `batch`): measures and corrects an
  additive residual sky level per channel, modelled as `β·T(λ)` with `T` the
  continuum of the sky that was subtracted from each frame
  (`SKYMODEL_ZAP`, `SKY`, or `NOSKYSUB − data`). Constraints: blue–red
  continuity across the dichroic plus optional absolute anchors from HST,
  Legacy Survey (nanomaggies), Pan-STARRS or another IFU cube; weighted
  least squares; `--groups` for one β per pointing; writes `<frame>_skl_icubes.fits`
  and new `.list` files for restacking.
- `kcwi_dichroic_step`: model-free check of the additive step between the
  blue and red stacks (`F_blue = k F_red + c`).
- `kcwi_cube_windows`: collapses a large cube (e.g. MUSE) into windowed median
  images for `kcwi_skylevel --cubewin`.
- `docs/scripts_instruction.md` entries for the three tools.

## 2026-11-02 — frame alignment without a point source (`7d23dd9`, YD)

### Added
- `kcwi_xalign` (`shifts` / `wcs` / `batch`): registers every blue and red
  frame on a common 0.1" grid from white-light images next to the dichroic,
  using leave-one-out templates with the (N−1)/N bias correction, quadratic
  peak refinement with Hessian errors, an overlap-fraction guard and
  coarse-to-fine search; writes both `.shift.list` files (with timestamped
  backups), a JSON record and a diagnostic PNG; `wcs` checks the two stacks
  against each other and writes `<stack>_icubes_wcs.fits`.
- `kcwi_shiftcheck`: judges one or two `.shift.list` sets by placing frames
  exactly as `kcwi_stack` would and measuring per-frame residuals.
- `kcwi_align`: several `align_box` lines in the `.par` file are
  cross-correlated together (`box_combine` = standardised / median / sum);
  single-box files behave as before. `kcwi_astrometry` gained `nhdu`.

## 2026-10-07 … 10-22 — automatic masks and split-frame cosmic-ray rejection (`f9646da`, `f466a57`, `b9124e3`, YD)

### Added
- `kcwi_automask`: automatic sky masks (`_intf` → `.reg` → `_smsk.fits`),
  median-filter masks (`thum`), ZAP source masks for KSkyWizard (`--wl`) and
  astroscrappy protection masks for the DRP (`--bcrmsk` →
  `<base>_intf_bcrmsk.fits`); polygon or ellipse regions, `--convert` to rebuild
  a mask from a hand-edited `.reg`, `--no-clobber`, `--json` grouping, `--plot`.
- `kcwi_stack(crr_split=True)`: cosmic-ray rejection on the common stacking
  grid by comparing groups of frames (`crr_split_mode` = `'auto'`, `'frame'`,
  `'half'`, `'alternate'`, N or custom labels), with per-frame
  `.crmsk.fits` / `.premsk.fits` masks validated by a grid signature so a
  stale mask is never applied, tile-compressed mask output, memory-aware
  wavelength blocking and a `.crrsplit.fits` diagnostic.
- `kcwikit/kcwi/kcwi_combine.py` — `kcwi_combine_stacks`: aligns and co-adds
  stacked cubes (NCC alignment through WCS, union or reference grid,
  exposure or inverse-variance weighting, optional flux scaling, the same
  CRR as above).

## 2025-09-19 (`95de13b`, YD)

- Removed the duplicated `kcwi_stack_yd` / `kcwi_check_samewave` /
  `kcwi_resample_wave` copies introduced in `7875693`. (This commit also
  reverted Yuguang's `kcwi_stack` fixes by accident; restored in November
  2026, see above.)
- `kcwi_stack`: mask threshold `hdu_m > 0.5` against interpolation artefacts;
  `plot=` keyword on `kcwi_resample_wave`.
- `kcwi_medfilter`: spline fit needs > 20 good points (was 4), linear instead
  of nearest interpolation elsewhere, prints which masks were found.

## 2025-06-05 … 06-06 — documentation (YD)

- `Running_blue_DRP.md`: step 3 on CR false positives and the `bcrmsk`
  protection mask. `Running_red_DRP.md` updates. Merge of `master` into `dev`.

## 2025-05-12 — red-channel masking tools (`7875693`, YD)

### Added
- `kcwi_masksky_ds9_bcrr`: DS9 regions → astroscrappy protection mask
  (`_intf_bcrmsk.fits`).
- `kcwi_mask_cube`: DS9 regions (`<frame>_rm.reg`) → `FLAGS` bit 32 in a copy
  of the cube, for removing pixels from the final stack.
- `kcwi_gen_skyfile`: builds the DRP `kcwi.sky` table from the frame-grouping
  JSON.
- `kcrm_group_frames --blue`; `kcwi_makemask_medfilter --red -f <json>` reuses
  one region file per pointing group; `kcwi_combinestd --notelluric`.
- `empty_hdul_like` (deep copy of an HDUList's structure).

### Changed
- `kcwi_vachelio`: default `method='barycentric'`, skips frames that already
  carry `VCORR`, no extrapolation of science spectra beyond the grid,
  `nearest-up` resampling for MASK/FLAGS.
- `kcrm_create_crmsk`: a cosmic-ray pixel must also exceed 60 % relative
  excess over the median frame.
- `kcwi_tools.cart2cyli`: masked voxels are NaN instead of 0.

### Fixed
- `kcwi_combinestd`: `np.zerors` typo.

## 2024-12-09 … 2025-04-29 — Python 3.12 / astropy 6 compatibility and stacking fixes (YC)

### Changed
- `pyregion` → astropy `regions` in `kcwi_masksky_ds9`, `_2d`, `_thum`,
  `_zap` and `setup.py` (`3d3f2c2`, `00bdeb5`, `6e4d6bc`); installation
  instructions for `python=~3.12`; the `astropy<6.0` pin removed.
- Tracked `__pycache__` files removed and ignored (`490cd6b`, `fabb9e5`,
  `052a930`).
- `pyDRP/`: local copy of the DRP modifications (`kcwi_file_primitives.py`,
  `kcwi_pipeline2.py` re-enabling calibration `process_*` entries,
  `MakeCube` warping MASK/FLAGS with `order=0`, `MakeInvsens` prompts at
  `plot_level >= 2`, `reduce_kcwi -st cr`) (`07eb186`, `2f1504e`).
- `pyDRP/primitives/FitCenter.py`: `.copy()` on the reference-spectrum slices
  before the in-place Tukey taper (fixes corruption of the shared atlas
  spectrum), `TUKEYALPHA` parameter (`0505ccd`).
- `pyDRP/primitives/RemoveCosmicRays.py`: `<frame>_crmsk_recover.reg` regions
  un-mask pixels of a `_crmsk.fits` (`ebfdaf1`, `e2c799d`). (Superseded in the
  `KCWI_DRP` fork by the protection mask.)

### Added
- `kcwi_stack(use_astrom=True)` reads a per-frame (3-column) `.astrom.list`
  for mosaics (`df92496`, `7d43e0d`).
- `kcwi_astrometry`: `background_subtraction`, `background_kcwi`,
  `background_ref` options (`485adca`); more robust reference-image reading
  (`427e63c`, `0505ccd`).
- `kcrm_create_crmsk --sigma` (`23e24d8`; only effective since Nov 2026).
- `kcwi_makemask` wrapper (globs `*.reg`).

### Fixed
- Montage issue #80 work-around in `kcwi_stack` (`ISS6FIX`, `37b6dd9`).
- `BUNIT` strings unfolded to `erg/s/cm2/Angstrom/arcsec2` (`4171316`).
- `(crls != 0)` guard so a zero plateau is not picked as the correlation peak
  in `kcwi_align` / `kcwi_astrometry` (`8e3408b`).
- Deprecated `independent_celestial_slices` removed from `reproject_interp`
  (`1a65f06`).
- `kcwi_check_flux` and `kcwi_astrometry` diagnostic figures are shown
  (`024e019`, `572d15b`).

## 2024-10-21 (`2ac9e17`, YC)

- `kcwi_gen_log`: the blue grating column held a literal list instead of
  `BGRATNAM`.

## 0.2 — 2024-09-11 (YC, ZZ)

- First packaged release: `setup.py` with console entry points, directory
  reorganisation (`kcwi/`, `scripts/`, `pyDRP/`, `arch/`), documentation set
  (installation, blue/red DRP instructions, post-DRP stacking, command-line
  tools, parameter-file reference), examples refreshed.
- GUI split off into the separate KSkyWizard tool; `kcwi_combinestd` supports
  KSkyWizard `invsens` files; `kcwi_makemask_medfilter` consolidated for
  `icube` / `icubed`; streamlined KCRM cosmic-ray rejection
  (`kcrm_group_frames`, `kcrm_create_crmsk`).
