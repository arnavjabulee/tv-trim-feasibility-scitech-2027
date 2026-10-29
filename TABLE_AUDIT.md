# Embedded lookup-table audit

Reviewed 29 October 2026. Scope: MATLAB and Python source in the public upload package, compared against the six inherited model-data/setup scripts in the original ICE model folder.

## Changes

- Removed `ukicematlab.m`, `enginedata.m`, `aerolookupdataall.m`, `constants.m`, `airdensity.m` and `airtemperature.m` from this package only. The original research copies are unchanged.
- Found one literal sideslip-breakpoint vector in `ice_init.m` matching the original lookup setup. Replaced it with a read of `BETA` from the separately supplied model. No aerodynamic coefficient values were embedded in that initialiser.
- Adapted both initialisers to use `load_ice_model_data.m`. It checks for all six external files before running them and reports an actionable missing-data error. Setup is described in EXTERNAL_MODEL.md.
- Added Git ignore rules for the inherited filenames and the local private_ice folder. No model data was copied into a replacement file or encoded into the loader.

## Inspection findings

The package now contains 24 MATLAB files and one Python helper. The aerodynamic routines `ice_block_aero_data_6dof.m` and `ice_block_aero_data_6dof_ctrl.m` interpolate coefficient tables from the externally populated `db.aero` structure. The engine routine `ice_block_engine.m` reads externally populated `db.eng` and `db.atm` tables. The initialisers collect variables from the external scripts rather than defining numerical coefficient tables.

The other inspected arrays define controller gains/weights and perturbations, operating-point cases, solver/search grids, state indices, coordinate transforms, or plot colours. No additional embedded aerodynamic, engine or atmosphere lookup table was identified. Standard physical constants and study settings are retained.

Checks included:

- Inspection of numeric array literals and table/initialisation code.
- Comparison of complete literal numeric arrays containing at least three numbers against the omitted source files; no remaining matches were found.
- Comparison of distinctive consecutive eight-number sequences against the omitted files, including differently formatted scalar sequences; no matches remained. The sequence check requires at least five distinct values and a non-integer value with magnitude above one to avoid common zero vectors and index lists.
- Search for literal coefficient/engine/atmosphere table assignments and copied model parameter values.
- Inspection of external-data loading and interpolation paths.

These are source-level checks, not proof of independent authorship or a legal clearance for every source line. No additional coefficient table was found by these checks. MATLAB was unavailable, so initialisation/error handling and numerical analyses have not been executed in MATLAB.

## Outputs retained

The manuscript figures, its 45 CSV inputs, current MATLAB output counterparts and selected supporting simulation datasets are retained as requested. These contain computed trim, contour and time-history results; they are not the omitted raw aerodynamic/engine coefficient database. This audit does not infer an output-redistribution restriction merely from use of an external model.

The LaTeX figures and datasets are retained; LQR weights and eigenvalues in the source were subsequently populated from the supplied final report. The six external model files are needed to regenerate MATLAB results, but not to read archived results or build the manuscript.
