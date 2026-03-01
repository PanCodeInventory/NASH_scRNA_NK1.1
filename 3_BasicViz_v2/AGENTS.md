# MODULE KNOWLEDGE BASE

## OVERVIEW
R-based visualization stage that consumes edited `.h5ad` data and exports publication-style PDFs.
Entry point is `main.R`; it dispatches plot scripts using a Snakemake-like runtime object.

## INPUT DATA PROFILE
Expected default input: `../2_ClusterEditor/results/renamed/merged.processed.pruned.renamed.h5ad`.

| Slot | Required keys for stage 3 |
|------|----------------------------|
| `obs` | `cluster_annot`, `time`, `sample`, `group` |
| `obsm` | `X_umap` |
| `X`/`layers` | normalized matrix in `X`; raw counts in `layers['counts']` if needed |
| `var` | gene index compatible with marker lookup from `config.yaml` |

Observed handoff profile (current dataset): shape `17550 x 21500`, `cluster_annot` has 5 classes, `X_umap` is `17550 x 2`.

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| Pipeline entry | `main.R` | loads `config.yaml`, builds `snakemake` object, runs scripts |
| Plot config | `config.yaml` | plot list, marker sets, palettes, theme defaults |
| Plot implementations | `scripts/plot_umap.R`, `plot_proportions.R`, `plot_dotplot.R`, `plot_feature.R`, `plot_river.R` | one script per plot type |
| Shared utilities | `scripts/utils.R` | plotting helpers and data loading helpers |
| Environment setup | `setup_env.sh` | creates `basicviz` conda env and export instructions |

## CONVENTIONS
- Run with working directory `3_BasicViz_v2/` so `config.yaml` and `scripts/` resolve.
- Keep plot behavior config-driven; add/modify plot variants in `config.yaml` before editing script logic.
- Preserve the `snakemake` object contract (`input`, `output`, `params`, `config`, `log`) used by sourced scripts.
- Maintain reticulate environment handoff from `setup_env.sh` (`RETICULATE_PYTHON`, optional `LD_LIBRARY_PATH`).

## ANTI-PATTERNS
- Do not edit generated PDFs or `results/colors.yaml` as source-of-truth.
- Do not call scripts directly without `main.R` context unless you recreate required `snakemake` slots.
- Do not assume system Python; pin to the configured conda python for reticulate compatibility.
- Do not move script paths without updating `main.R` dispatcher and `config.yaml` references.

## COMMANDS
```bash
# setup compatible runtime once
cd "3_BasicViz_v2" && bash setup_env.sh

# run full visualization pipeline
cd "3_BasicViz_v2" && Rscript main.R

# inspect output artifacts
cd "3_BasicViz_v2" && ls results
```

## NOTES
- Stage-3 input is normally the renamed output from `2_ClusterEditor/results/renamed/`.
- This module has no standalone test suite; validation is successful `Rscript main.R` completion plus expected PDFs.
