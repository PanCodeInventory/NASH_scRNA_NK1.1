# MODULE KNOWLEDGE BASE

## OVERVIEW
Snakemake workflow for post-clustering edits: prune clusters, rename labels, and optional reprocessing.
Canonical execution directory is `2_ClusterEditor/workflow/`.

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| Workflow orchestration | `Snakefile` | rules: `check_env`, `prune_clusters`, `rename_clusters`, `reprocess_clusters` |
| Script entrypoints | `scripts/00_check_env.py`, `01_prune_clusters.py`, `02_rename_clusters.py`, `03_reprocess_clusters.py` | staged execution flow |
| Shared logic | `scripts/utils.py` | config load, AnnData I/O, mapping/reprocess helpers |
| Runtime config | `../config/workflow_config.yaml` | prune/rename/reprocess sections and plot options |
| Dependencies | `../requirements.txt` | scanpy/anndata/pandas/matplotlib/pyyaml |

## CONVENTIONS
- Run from `workflow/`; Snakefile assumes config path `config/workflow_config.yaml` (or fallback via `../`).
- Keep section names stable: `prune`, `rename`, `reprocess`.
- Respect mapping precedence: `merge_rules` overrides direct `cluster_mapping` when both apply.
- Keep env check as first gate; downstream rules depend on `logs/env_checked.flag`.

## ANTI-PATTERNS
- Do not remove config keys unless intentionally refactoring parser logic.
- Do not hardcode output paths in scripts; rely on config-driven paths under `../results`.
- Do not edit generated outputs (`../results`, `../logs`, `../.snakemake`) directly.
- Do not run workflow from repo root; relative config/input paths will break.

## COMMANDS
```bash
# full cluster-edit workflow
cd "2_ClusterEditor/workflow" && snakemake --cores 4

# dry-run validation
cd "2_ClusterEditor/workflow" && snakemake --dryrun

# run single rule for iterative edits
cd "2_ClusterEditor/workflow" && snakemake rename_clusters --cores 2
```

## NOTES
- Inputs typically originate from stage 1 merged output; keep handoff paths aligned across modules.
- `scripts/utils.py` contains compatibility cleanup before AnnData save; keep that path intact for cross-version safety.
