# Matched-flat controls for shortest-path diagnostics

This repository contains the Wolfram Language notebook and script used to generate the original matched-flat control benchmark for shortest-path anisotropy diagnostics in geometric graphs from black-hole embeddings.

The main diagnostic is `CLog`, a graph-shell statistic based on the dispersion of logarithmic shortest-path multiplicities. The code compares black-hole embedding graphs with matched-flat controls and exports the tables used in the manuscript.

## Repository structure

```text
.
├── notebooks/
│   └── original_matched_flat_control_step4_clean.nb
├── scripts/
│   └── original_matched_flat_control_step4_clean.wl
├── results/
│   ├── figures/
│   └── tables/
├── README.md
└── .gitignore
```

## Requirements

- Wolfram Mathematica / Wolfram Language 14.3 or newer recommended.
- No external Wolfram packages are required for the core run.

## How to run

### Option 1: Notebook

Open

```text
notebooks/original_matched_flat_control_step4_clean.nb
```

and evaluate the notebook from top to bottom.

The notebook has been cleaned for GitHub: generated output cells, print cells, graphics boxes, cell labels, expression UUIDs, and notebook cache metadata have been removed.

### Option 2: Wolfram script

From the repository root, run:

```bash
wolframscript -file scripts/original_matched_flat_control_step4_clean.wl
```

The script exports generated outputs into directories named like

```text
step4_outputs_originalMatchedFlatControl_randomN1000_k16_epsFactor1p15_rg3/
step4_checkpoints_originalMatchedFlatControl_randomN1000_k16_epsFactor1p15_rg3/
```

The checkpoint directory stores `.mx` files and is ignored by Git.

## Main settings

The default run uses:

```wolfram
M = 1.0;
rMin = 2.2;
rMax = 12.0;
nSamples = 1000;
randomSeed = 1234;
kNN = 16;
epsilonFactor = 1.15;
fixedGraphShellRadius = 3;
nRadialBins = 12;
maxRefsPerBin = All;
```

For a quick smoke test, set:

```wolfram
maxRefsPerBin = 5;
rgSweepMaxRefsPerBin = 20;
sensitivityMaxRefsPerBin = 20;
```

For the final manuscript values, use `All`.

## Generated manuscript tables

The script exports the main correlation and sensitivity tables, including:

```text
graph_summary_original_control.csv
all_geometry_correlation_gap_table_original_control.csv
all_geometry_amplitude_shape_gap_table_original_control.csv
extra_sensitivity_checks/rg_sweep_clog_only_schwarzschild_original_control.csv
extra_sensitivity_checks/control_protocol_sensitivity_CLog_only.csv
extra_sensitivity_checks/control_protocol_sensitivity_table_schwarzschild.csv
```

The final full-reference control-protocol sensitivity result for `CLog` is:

```text
Original matched-flat KNN:             |Δρr| ≈ 0.781
Strict radius-threshold matched flat:  |Δρr| ≈ 0.001
```

This is the key methodological sensitivity result: the apparent black-hole/control separation is strongly dependent on the graph-construction and matched-control convention.

## Notes for contributors

- Do not commit generated checkpoint files (`*.mx`).
- Do not commit large generated output folders unless they are intentionally archived as a release artifact.
- If you regenerate figures or tables for a paper submission, place final curated versions under `results/figures/` and `results/tables/`.

## Citation

A manuscript based on this code is in preparation. If you use the code before publication, please cite the repository and the associated arXiv/preprint when available.
