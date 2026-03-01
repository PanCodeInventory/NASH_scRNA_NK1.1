# PROJECT KNOWLEDGE BASE

**Generated:** 2026-02-26 15:27:14 CST  
**Branch:** `main`  
**Commit:** `d163621`

## OVERVIEW
Single-cell RNA-seq analysis pipeline for NK1.1 with three ordered modules: matrix-to-markers, cluster editing, and visualization.
Repository is artifact-heavy; active code lives mainly in stage/workflow script directories.

## STRUCTURE
```text
NK1.1/
├── 1_Matrix2Markers/      # Stage 1: 10x ingest, QC, merge, marker discovery
├── 2_ClusterEditor/       # Stage 2: prune/rename/reprocess clusters
├── 3_BasicViz_v2/         # Stage 3: R plotting and figure export
└── Files/                 # data artifacts (gitignored)
```

## DATA PROFILE
Primary handoff dataset: `Files/merged.processed.pruned.renamed.h5ad`.

| Item | Value |
|------|-------|
| Shape (`n_obs x n_vars`) | `17550 x 21500` |
| `obs` key columns | `sample`, `group`, `time`, `leiden`, `cluster_annot`, QC fields |
| `var` columns | `highly_variable`, `highly_variable_rank`, `means`, `variances`, `variances_norm` |
| `obsm` keys | `X_pca` (`17550 x 50`), `X_pca_harmony` (`17550 x 50`), `X_umap` (`17550 x 2`) |
| `layers` keys | `counts` (`int64`, `17550 x 21500`) |
| `uns` keys | `neighbors`, `leiden`, `hvg`, `log1p`, color maps |
| `raw` present | `true` |

Cluster annotation distribution in `cluster_annot`: `Cytotoxic NK` 6103, `ILC1` 4741, `Circulating NK` 3745, `High-Ribosomal NK` 1666, `Homeostatic NK` 1295.
Sample distribution in `sample`: `NCD_0W` 6525, `MCD_6W` 3935, `MCD_1W` 3907, `MCD_2W` 3183.

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| Stage 1 workflow logic | `1_Matrix2Markers/scripts/Snakefile` | orchestrates stages `01`-`05` |
| Stage 1 parameters | `1_Matrix2Markers/scripts/config.yaml` | samples, QC thresholds, processing params |
| Stage 1 utilities | `1_Matrix2Markers/scripts/utils/io_and_qc_utils.py` | central helper module |
| Stage 2 workflow logic | `2_ClusterEditor/workflow/Snakefile` | includes env check + 3 editing steps |
| Stage 2 parameters | `2_ClusterEditor/config/workflow_config.yaml` | prune/rename/reprocess config |
| Stage 2 utilities | `2_ClusterEditor/workflow/scripts/utils.py` | config + AnnData manipulation |
| Stage 3 entrypoint | `3_BasicViz_v2/main.R` | constructs runtime object and dispatches plots |
| Stage 3 plotting config | `3_BasicViz_v2/config.yaml` | plot list, markers, colors |
| Stage 3 setup | `3_BasicViz_v2/setup_env.sh` | conda env + reticulate export guidance |

## CODE MAP
| Symbol/Rule | Type | Location | Role |
|-------------|------|----------|------|
| `read_10x`/`qc_filter`/`merge_and_process` | Snakemake rules | `1_Matrix2Markers/scripts/Snakefile` | core processing DAG |
| `prune_clusters`/`rename_clusters`/`reprocess_clusters` | Snakemake rules | `2_ClusterEditor/workflow/Snakefile` | cluster post-processing DAG |
| `run_script()` | R function | `3_BasicViz_v2/main.R` | orchestrates per-plot script execution |
| `load_10x_mouse()` | Python function | `1_Matrix2Markers/scripts/utils/io_and_qc_utils.py` | robust 10x loading |
| `apply_cluster_mapping()` | Python function | `2_ClusterEditor/workflow/scripts/utils.py` | cluster label remap |

## CONVENTIONS
- Execute each module from its own working directory (`scripts/`, `workflow/`, or `3_BasicViz_v2/`).
- Keep workflow parameters in YAML; scripts consume config keys rather than hardcoded constants.
- Stage scripts use numbered filenames (`01_...`, `02_...`) to mirror pipeline order.
- Outputs are expected under `results/` and logs under `logs/`; these are generated artifacts.
- 2_ClusterEditor keeps utility logic centralized in `workflow/scripts/utils.py` and thin stage scripts.

## ANTI-PATTERNS (THIS PROJECT)
- Do not edit generated artifacts in `results/`, `logs/`, `.snakemake/`, or `Files/`.
- Do not run workflow commands from repo root; relative config paths assume module-local working dirs.
- Do not break YAML key structure in `workflow_config.yaml` and `config.yaml` files.
- Do not treat duplicated nested path `1_Matrix2Markers/scripts/scripts/` as canonical unless explicitly needed.
- Do not assume CI or unit tests exist; validation is workflow dry-run + stage execution.

## UNIQUE STYLES
- Mixed Python + R pipeline with Snakemake-driven preprocessing and R-based figure generation.
- Artifact-first bioinformatics layout: small code footprint, very large `.h5ad` outputs.
- Cluster editing explicitly separates prune, rename, and reprocess as independent stages.

## COMMANDS
```bash
# Stage 1
cd "1_Matrix2Markers/scripts" && snakemake --cores 4

# Stage 2
cd "2_ClusterEditor/workflow" && snakemake --cores 4

# Stage 3
cd "3_BasicViz_v2" && Rscript main.R

# Safe validation pass
cd "1_Matrix2Markers/scripts" && snakemake --dryrun
cd "2_ClusterEditor/workflow" && snakemake --dryrun
```

## NOTES
- No repository CI workflow files were found.
- No dedicated unit/integration test suite was found.
- `.gitignore` intentionally excludes primary data/output paths; rebuild artifacts from workflows.
