# Supplying the ICE model locally

The inherited aerodynamic lookup tables, engine tables, atmosphere data and original constants/lookup setup scripts are intentionally absent from this repository. Obtain them separately from an authorised source. This repository does not supply or download them, and does not grant rights to them.

The required files are:

- `ukicematlab.m`
- `enginedata.m`
- `aerolookupdataall.m`
- `constants.m`
- `airdensity.m`
- `airtemperature.m`

Place the six compatible files together in a local folder. From the repository's `matlab/` directory, configure MATLAB with:

```matlab
setenv('ICE_MODEL_DIR', 'C:\path\to\authorised\ICE-model');
db = ice_init();
% Or initialise the aerodynamic-control branch:
db_aero = ice_init_aeroctrl();
```

Use the actual folder on your computer. A POSIX absolute path works on Linux/macOS. Alternatively, create `matlab/private_ice/` and place the files there; this is the default if `ICE_MODEL_DIR` is unset. Both the folder and the inherited filenames are ignored by Git.

The loader checks that all six files exist before executing any of them, reports missing dependencies with error identifier `ICE:MissingModelData`, and runs the separately supplied scripts in MATLAB's base workspace as the original initialisation did. It temporarily adds their folder to the MATLAB path for their internal script calls and restores the previous path afterward.

The files must define the variables expected by `ice_init.m` and `ice_init_aeroctrl.m`, including aerodynamic and effector tables, engine tables, geometry, mass/inertia, and atmosphere vectors. An arbitrary ICE model release is not guaranteed to be compatible. No replacement or fabricated coefficients are supplied.

The included manuscript figures and archived simulation results are retained. Reading those outputs or building the manuscript does not require the omitted MATLAB scripts. Regenerating the analyses requires the separately obtained model and the MATLAB toolboxes described in the README.
