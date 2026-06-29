# journal_benchmark_version
# Outputs

This folder stores generated diagnostic tables.

Large generated files are not tracked by Git.

Expected Step 4 outputs:

- `step4_source_vertex_diagnostics.csv`
- `step4_radial_bin_diagnostics.csv`
- `step4_source_vertex_diagnostics.mx`
- `step4_radial_bin_diagnostics.mx`

The source-vertex table contains one row per reference vertex.
The radial-bin table aggregates diagnostics by geometry and radial bin.

# Figures

This folder stores generated figures for the benchmark extension.

Planned Step 4 figures:

1. `comparison_schwarzschild_flamm_baselines.pdf`
2. `comparison_matched_flat_controls.pdf`
3. `diagnostic_correlation_summary.pdf`

These figures compare:

- the main shortest-path anisotropy statistic \(C_{\log}\)
- shell-size diagnostics
- mean log shortest-path multiplicity
- local degree
- local clustering coefficient

# Benchmark Extension Plan

The current extension reframes the shortest-path anisotropy statistic as part of a benchmark protocol for graph-based curvature diagnostics.

## Goal

Compare the main statistic \(C_{\log}\) against simple graph baselines to test whether the observed radial signal is explained by:

1. shell size,
2. average shortest-path multiplicity,
3. local graph density,
4. local clustering.

## Graphs

The benchmark includes:

- Schwarzschild/Flamm
- Reissner–Nordström
- Bardeen
- Hayward
- matched-flat controls for each geometry

## Diagnostics

For each graph and each radial bin, compute:

- \(C_{\log}\)
- mean shell size
- standard deviation of shell size
- mean log shell size
- mean log shortest-path multiplicity
- standard deviation of log shortest-path multiplicity
- degree
- local clustering coefficient

## Planned next step

Compute radial correlations between each diagnostic and the logarithmic Kretschmann profile, then compare curved geometries against their matched-flat controls.


## Current extension: benchmark diagnostics

This repository is being extended to benchmark the shortest-path anisotropy statistic against simpler graph diagnostics. The goal is to test whether the observed radial signal in black-hole graph geometries is specific to the proposed statistic or can be explained by simpler quantities such as shell size, shortest-path multiplicity, degree, or clustering.

The current benchmark includes:

- Schwarzschild/Flamm graph geometry
- Reissner–Nordström graph geometry
- Bardeen graph geometry
- Hayward graph geometry
- matched-flat controls

The main notebook for this extension is:

```text
notebooks/step4_baseline_diagnostics.nb
