# MODULE KNOWLEDGE BASE

## OVERVIEW
Snakemake-driven stage-1 pipeline: 10x ingest, single-sample QC, merge/process, marker discovery, report generation.
Canonical working directory is `1_Matrix2Markers/scripts/`.

## WHERE TO LOOK
| Task | Location | Notes |
|------|----------|-------|
| Workflow DAG | `Snakefile` | rules `read_10x`, `qc_filter`, `merge_and_process`, `find_markers`, `generate_report` |
| Runtime parameters | `config.yaml` | sample paths, filtering, processing, marker params |
| Stage scripts | `stages/01_*.py` to `stages/05_*.py` | one script per pipeline stage |
| Shared utilities | `utils/io_and_qc_utils.py` | loading, QC helpers, plotting save helpers |
| Environment lock | `environment.yaml` and `requirements.txt` | dependencies for snakemake + scanpy stack |

## CONVENTIONS
- Treat `scripts/` as canonical; duplicated `scripts/scripts/` exists but should not be primary edit target.
- Keep stage ordering aligned with filename prefixes `01`-`05`.
- Preserve YAML key schema in `config.yaml`; stage scripts expect exact key paths.
- Let Snakefile route CLI args into stages; avoid hardcoding thresholds in Python scripts.
- Write generated outputs only to configured `results/`/`reports/`/`logs/` paths.

## ANTI-PATTERNS
- Do not edit `.h5ad`, `.png`, `.csv`, `.tsv`, or generated HTML under `../results` and `../reports`.
- Do not run `snakemake` from repo root; run from this directory so relative paths resolve.
- Do not switch canonical edits to `scripts/scripts/` unless explicitly performing duplication cleanup.
- Do not bypass config by embedding sample-specific absolute paths into stage scripts.

## COMMANDS
```bash
# full stage-1 run
cd "1_Matrix2Markers/scripts" && snakemake --cores 4

# inspect planned jobs only
cd "1_Matrix2Markers/scripts" && snakemake --dryrun

# run a single stage rule
cd "1_Matrix2Markers/scripts" && snakemake find_markers --cores 2
```

## NOTES
- Config comments are bilingual; keep keys and value structure stable even when comment language changes.
- `io_and_qc_utils.py` is high-centrality; changes ripple across multiple stage scripts.
