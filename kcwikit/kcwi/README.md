# kcwikit.kcwi

Post-DRP stacking and combination of KCWI data cubes.

- `kcwi.py` - `kcwi_align`, `kcwi_stack` (incl. `crr_split` cosmic-ray rejection), `kcwi_astrometry`, `kcwi_vachelio`, `kcwi_stack_readpar` and friends.
- `kcwi_combine.py` - `kcwi_combine_stacks`: align and co-add several stacked cubes.
- `kcwi_regions.py` - `regions_to_mask`: DS9 region files (image, physical or sky coordinates, with `-` exclusions) to boolean masks; used by every script that reads a `.reg` file.
- `kcwi_tools.py` - cube manipulation helpers (collapse, cylindrical projection, ...).
- `kcwi_pyDRP.py`, `kcwi_pypeit.py` - legacy variants kept for reference; not maintained.

Requirements: Python >= 3.11 (3.12 recommended), `Montage` (command-line tools or `MontagePy`), and the packages installed by `pip install .` at the repository root.

Usage is described in [docs/KCWI_post-DRP_stacking.md](../docs/KCWI_post-DRP_stacking.md); the parameter file in [docs/parameter_file_reference.md](../docs/parameter_file_reference.md).
