# 4_GenesetScore

This module is reorganized into two independent tracks:

- `geneset_scoring/`: score cells/clusters using multi-gene signatures.
- `single_gene_scoring/`: visualize or score individual genes independently.

## Structure

- `4_GenesetScore/geneset_scoring/config/genesets.yaml`
- `4_GenesetScore/geneset_scoring/scripts/run_geneset_score.py`
- `4_GenesetScore/geneset_scoring/docs/`
- `4_GenesetScore/geneset_scoring/results/`
- `4_GenesetScore/single_gene_scoring/xcl1_violin_cluster_annot/`

## Quick Run

```bash
# Gene-set scoring
python "4_GenesetScore/geneset_scoring/scripts/run_geneset_score.py"

# Single-gene plotting
python "4_GenesetScore/single_gene_scoring/xcl1_violin_cluster_annot/plot_xcl1_violin.py"
```
