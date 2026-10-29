# Paper input provenance

`latex/main.tex` directly reads **45 CSV files** and **two PDF illustrations**. MATLAB files are not executed by LaTeX: their relevance is as generators, related analyses, export helpers or dependencies. Relationships below distinguish explicit filename evidence from inferred workflow relationships.

All input CSVs have counterparts in `ICEmodel - TVC`. The MATLAB counterparts are copied to `matlab/paper_outputs/`; the exact manuscript versions are in `latex/`. `OUTPUT_COMPARISON.csv` records their comparison. `FILE_MANIFEST.csv` records the original location and hash of each copied file.

## Data-to-code relationships

| Paper input | Related saved code | Evidence and limits |
| --- | --- | --- |
| `branch_445_main.csv`, `branch_445_island.csv` | `continue_trim_branches_445_tvc.m` | Both output filenames are explicitly assigned in the script. |
| `lqr_tvc_main_envelope_vt_250_alt_6000.csv`, `lqr_tvc_island_vt_480_alt_18000.csv`, `lqr_tvc_unstable_1_vt_450_alt_16000.csv` | `run_lqr_tvc_cases.m` | Export naming pattern and case names/speeds/altitudes match. Island and unstable cases are currently commented out. The current main-envelope output differs from the paper copy. |
| `lqr_turn_vt_300_alt_0.csv` | `simulate_saved_turn_trim_lqr_6dof.m` | Explicit export pattern matches. Pass target speed 300 and altitude 0; current defaults select a different case. The saved turn-trim input is included. |
| `aero_pla_lower.csv`, `aero_pla_upper.csv` | `make_clipped.py` | Both output names and input `aero_pla_127.csv` are explicit. This is a Python processing step, included alongside MATLAB because the paper uses its named outputs. Exact historical execution is not recorded. |
| `pla_80.csv`, `pla_100.csv`, `pla_110.csv`, `pla_127.csv`, `pla_132.csv`, `pla_140.csv`, `pla_150.csv`, `pla_200.csv`, `pla_300.csv`, `pla_500.csv`, `pla_1000.csv` | `export_contours_to_csv.m`, `split_contour_levels.m`, `sweep_trim_vt_6dof_equilibrium_map.m` | Export helpers generate contour tables and `<prefix>_<level>.csv` filenames; the map script supplies the related TVC equilibrium data. Calls, selected contour levels and assembly of the full grid were not saved as a complete driver. |
| `dzneg_m1.csv`, `dzneg_m5.csv`, `dzneg_m10.csv`, `dzneg_m15.csv`, `dzneg_m20.csv`, `dzneg_m25.csv`, `dzneg_m28.csv` | Same contour helpers and TVC equilibrium map | Split helper explicitly uses `<prefix>_m<abs(level)>.csv` for negative contours. The historical `dzneg` helper invocation is not saved. |
| `aero_pla_80.csv`, `aero_pla_100.csv`, `aero_pla_110.csv`, `aero_pla_127.csv`, `aero_pla_132.csv`, `aero_pla_140.csv`, `aero_pla_150.csv`, `aero_pla_200.csv`, `aero_pla_300.csv`, `aero_pla_500.csv` | Same contour helpers and `sweep_trim_vt_6dof_aero_equilibrium_map.m` | Related aerodynamic equilibrium-map and contour workflow. Historical export calls and full-grid assembly are not saved. |
| `final_m15.csv`, `final_m20.csv`, `final_m25.csv` | `sweep_trim_vt_6dof.m` and TVC equilibrium/contour workflow | Related vector-limit sweeps and raw datasets exist. No saved MATLAB or Python source names these final output files; their final clipping/splicing/assembly is untraced. They are preserved as archived paper data. |
| `pla_upper.csv` | TVC equilibrium/contour workflow | Curve is present, but no saved source explicitly names this output. Upper-branch extraction is untraced. |
| `deltaz_alpha_altitude_thinned.csv` | TVC equilibrium-map workflow | Raw map contains related trim quantities. No saved source explicitly writes the thinned file. Thinning and selection are untraced. |
| `gamma_vt_tau_max_results_sorted.csv` | `run_gamma_vt_tau_max_plots.m` | Script explicitly writes the unsorted `gamma_vt_tau_max_results.csv`; the sorted paper counterpart exists. Sorting/filtering steps are untraced, and the current raw file is larger than the paper input. |
| `vn_curve_0m.csv`, `vn_curve_10000m.csv`, `vn_curve_20000m.csv` | `tvc_turn_envelope_6dof_tas.m` | Script produces the related raw coordinated-turn envelope. A processed combined file with `EAS_mps` and `n_max_smoothed` columns exists and is included. No saved source writes the individual `vn_curve_*` files or records the smoothing procedure. |

## Other paper content

- `ICE_silhouette.pdf` and `ICE_roll.pdf` are directly included illustrations. Both are copied with the manuscript. No MATLAB generator was found for them.
- `make_table3_eigs.m` targets the four operating points in the manuscript's eigenvalue table. The LQR weights and table values in `main.tex` were checked against the supplied final report, *IEEE_Conference_Template (4).pdf*, pages 5-6; they were not regenerated during packaging.
- `sweep_trim_vt_6dof_aero.m` supports the nominal aerodynamic-control methodology. It does not directly write the plotted `aero_pla_*` filenames.
- The class, bibliography style and bibliography are retained through the manuscript's direct and transitive dependencies.

## Dependency selection

The twelve related entry points/export/table helpers were followed through references to local MATLAB filenames and function/script names, including initialisation via `run`/`evalin`. The initial selection contained 29 `.m` files. The six inherited model-data/setup scripts were subsequently removed from the upload package, leaving 23 original MATLAB source files plus a new external-data loader (24 total). `ice_init.m` and `ice_init_aeroctrl.m` were adapted to call that loader; the former now also reads the sideslip grid from external data. See EXTERNAL_MODEL.md and TABLE_AUDIT.md. Numerical methods and settings are retained; comments and descriptive text use British spelling.

The selected scripts use MATLAB functions directly; no reference to the saved Simulink model was needed by this subset. `icemodel.mdl`, generated `slprj`/`.slxc` content, `.asv` files, emergency workspaces and the large thrust-bridge checkpoint were therefore excluded. The original and duplicate TVC model directories were not copied.

Eleven supporting datasets are retained separately: the three nominal vector-limit sweeps; full TVC and aerodynamic equilibrium maps; PLA, aerodynamic PLA and negative vector-angle contour tables; the current raw recovery-metric CSV; and raw/processed coordinated-turn envelopes. These support investigation of the processing chain; they are not additional files read directly by LaTeX.

## Verification performed

- All initial copies were hash-checked. The current manifest records current hashes, source hashes and whether each file is unchanged, modified or newly authored; excluded files are listed separately in EXCLUDED_FILES.csv.
- Every directly referenced CSV and PDF in the manuscript is present in `latex/`.
- All 45 paper CSVs were compared to their MATLAB-directory counterparts; one differs, and both versions are preserved.
- Local MATLAB source dependencies were checked. The six omitted scripts are intentional external dependencies; the loader checks for them before use. No remaining copied lookup-table arrays were identified in the source audit.
- MATLAB was not available on the execution path, so no numerical regeneration was performed. A complete end-to-end reproduction and a successful manuscript build have not been established.
