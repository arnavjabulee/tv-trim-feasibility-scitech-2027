# Trim Feasibility Analysis of Thrust Vectoring and Aerodynamic Control in a Tailless Delta Wing Aircraft

LaTeX manuscript, MATLAB analysis and results for *Trim Feasibility Analysis of Thrust Vectoring and Aerodynamic Control in a Tailless Delta Wing Aircraft*.

- `latex/`: manuscript (`main.tex`), references and figure data.
- `matlab/`: trim searches, control analysis and supporting outputs.

Build the manuscript from `latex/` with LaTeX and TikZ shell escape enabled.

MATLAB analyses require Optimization Toolbox and Control System Toolbox. The original ICE lookup tables are excluded; see [model setup](EXTERNAL_MODEL.md). Run analyses from `matlab/` after configuring the model.

See [provenance](PROVENANCE.md) for the data-to-code mapping and reproduction notes.
