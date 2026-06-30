# Manifest

## Cleaned files

- `notebooks/original_matched_flat_control_step4_clean.nb`  
  Clean Mathematica notebook with outputs, print cells, graphics boxes, labels, UUIDs, and cache removed.

- `scripts/original_matched_flat_control_step4_clean.wl`  
  Reproducible Wolfram Language script version. Includes the main original matched-flat benchmark and the final sensitivity checks used in the manuscript.

## Generated output policy

Generated output directories are ignored by default:

- `step4_outputs_*`
- `step4_checkpoints_*`
- `extra_sensitivity_checks/`

Curated final outputs for the paper can be copied manually into `results/tables/` and `results/figures/` if you want them tracked in Git.
