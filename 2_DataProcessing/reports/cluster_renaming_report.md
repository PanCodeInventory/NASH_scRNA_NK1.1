# Cluster Renaming Report

## Renaming Strategy
- **Cluster 5 → Cluster 6**: B cell contamination cluster (Iglc3, Cd79a, Iglc2, Ebf1, Ms4a1, etc.)
- **Cluster 6 → Cluster 5**: Proliferation cluster (H2afx, Cdca3, Hist1h3c, Mki67, etc.)

## Cluster Summary

| Original Cluster | Renamed Cluster | Cell Count |
|------------------|-----------------|------------|
| 0 | 0 | 4704 |
| 1 | 1 | 0 |
| 2 | 2 | 0 |
| 3 | 3 | 0 |
| 4 | 4 | 0 |
| 5 | 5 | 0 |
| 6 | 0 | 0 |
| 0 | 1 | 0 |
| 1 | 2 | 4631 |
| 2 | 3 | 0 |
| 3 | 4 | 0 |
| 4 | 5 | 0 |
| 5 | 0 | 0 |
| 6 | 1 | 0 |
| 0 | 2 | 0 |
| 1 | 3 | 0 |
| 2 | 4 | 3918 |
| 3 | 5 | 0 |
| 4 | 0 | 0 |
| 5 | 1 | 0 |
| 6 | 2 | 0 |
| 0 | 3 | 0 |
| 1 | 4 | 0 |
| 2 | 5 | 0 |
| 3 | 0 | 3885 |
| 4 | 1 | 0 |
| 5 | 2 | 0 |
| 6 | 3 | 0 |
| 0 | 4 | 0 |
| 1 | 5 | 0 |
| 2 | 0 | 0 |
| 3 | 1 | 0 |
| 4 | 2 | 1094 |
| 5 | 3 | 0 |
| 6 | 4 | 0 |
| 0 | 5 | 0 |
| 1 | 0 | 0 |
| 2 | 1 | 0 |
| 3 | 2 | 0 |
| 4 | 3 | 0 |
| 5 | 4 | 180 |
| 6 | 5 | 118 |

## Files Generated
- **Renamed RDS**: `2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.renamed.rds`
- **Original RDS**: `2_DataProcessing/RDS/nk.integrated.singleR_annotated.noCluster6.tuned.rds`
- **Cluster Summary**: `cluster_renaming_summary.csv`

## Next Steps
1. Update analysis result files (markers, proportions)
2. Generate visualizations with hidden cluster 6
3. Validate downstream analysis compatibility

---
*Generated: 2025-10-23 18:13:58*
