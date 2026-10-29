# Publication provenance and exclusions

Updated 29 October 2026 after removing the inherited model data from the public upload package.

## Current package

The package contains the manuscript, archived results and analysis code with an external model-data dependency. The six inherited files are absent: ukicematlab.m, enginedata.m, aerolookupdataall.m, constants.m, airdensity.m and airtemperature.m. Their names and hashes are recorded in EXCLUDED_FILES.csv without their contents. No model data is bundled in an alternative format. TABLE_AUDIT.md describes the inspection and the one embedded grid removed from ice_init.m. EXTERNAL_MODEL.md explains local setup.

The researcher states that the figures are already published through AIAA and asks to retain figures and outputs. They remain in this package. No blanket licence has been added or third-party ownership asserted.

## Model provenance

The researcher received the original model through an internal university email from research supervisors and identifies C. D. C. Jones as a previous holder. No explicit redistribution licence for the inherited data was identified.

The original ICEreport file is a PostScript report titled Description of the ICE 101-TV, attributed to C. D. C. Jones and dated 13 December 2003. The local thesis ChrisJones_thesis_2005.pdf identifies its author as Christopher D. C. Jones, MEng, with submission in September 2005 and revision in November 2005. Neither document is included in the upload package.

The thesis acknowledgements (printed page ii, PDF page 3) credit Bill Blake at the Air Force Research Laboratory, Wright-Patterson Air Force Base, with providing the data and information from which the ICE Simulink model was created. Appendix A.7 (printed page 238, PDF page 267) explains that only a selection of the aerodynamic data is included in the thesis. This is evidence of provenance, not a licence for the full model database.

[Bristol's publication record](https://research-information.bris.ac.uk/en/publications/tailored-dynamic-gain-scheduled-control/) lists CDC Jones, MH Lowenberg and TS Richardson as coauthors of Tailored Dynamic Gain-Scheduled Control. The [ICE 101-TV paper](https://www.cambridge.org/core/journals/aeronautical-journal/article/abs/dynamic-gainscheduled-control-of-the-ice-101tv/9B7CE75D81E24138DCD6A281F740F47C) documents related research. These publications do not establish a redistribution licence for the exact scripts.

Exact-filename web searches and GitHub code-search queries for ukicematlab and aerolookupdataall did not establish a licensed public copy. Search coverage is limited and absence of a result does not prove that no public release exists. The excluded files should remain external unless appropriate redistribution permission is established.

## Scope and remaining technical limits

- The remaining analysis code reads model coefficients from external data; no further embedded coefficient tables were found by the documented audit.
- The original research and LaTeX folders remain untouched by this exclusion step.
- One archived paper CSV differs from its current MATLAB counterpart; both versions remain identified in OUTPUT_COMPARISON.csv.
- Current script settings and missing historical processing steps prevent a claim of exact automated reproduction. MATLAB analyses have not been rerun.
- The LQR weights and eigenvalue table have been populated from the supplied final report. The revised equation and table compile successfully in isolation. The full-document check stops at the existing newtx/amssymb conflict (`\Bbbk` already defined); bibliography-style declarations also remain inconsistent.
- The bibliography style retains its LaTeX Project Public License notice; the class credits Overleaf. Existing notices are preserved. No claim of ownership over those components is made.
- A targeted scan found no common credential patterns or absolute user-home paths; this was not a comprehensive security audit.

[GitHub licensing guidance](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository) and [UK IPO guidance](https://www.gov.uk/using-somebody-elses-intellectual-property/copyright) explain general permission principles. This document records provenance and the chosen exclusions; it does not determine legal ownership of every component.
