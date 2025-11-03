#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(Matrix)
  library(DESeq2)
  library(dplyr)
  library(data.table)
})

# ------------------------------------------------------------
# 中文内联注释版（在原脚本基础上补充注释并支持“整体/簇分析”）：
# 目的：按时间点进行样本比对，既可做整体分析（不分簇与细胞类型），也可做簇分析（按 Seurat 的 seurat_clusters 分层）。
#
# 配置约定（scripts/config.yaml）：
# - columns.sample_id：建议设为 "timepoint"（无 orig.ident 时，以时间点作为样本ID用于伪样本聚合）
# - columns.cluster  ：建议新增并设为 "seurat_clusters"；本脚本的“簇分析”以此列做分层，不使用 singleR.label
# - columns.timepoint："timepoint"
# - columns.condition：如无条件列，可留空/缺失；本脚本将仅用 timepoint 构造分组标签
# - group_label_format：若无条件列，应设为 "{timepoint}"（例如 0W、1W、2W、6W）
# - deg.scope：分析范围，取值 "overall"（仅整体）、"by_cluster"（仅簇分析）、"both"（两者都做）。
# - comparisons：control/case 应与 group_label_format 生成的标签一致（如 "0W" vs "1W"）。
# 结果输出结构同原脚本，新增整体分析的 all_cells.* 文件。
# ------------------------------------------------------------

# ---------------------------
# CLI args
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)
cfg_path <- NULL
if (length(args) >= 2) {
  for (i in seq_along(args)) {
    if (args[i] == "--config" && i + 1 <= length(args)) {
      cfg_path <- args[i + 1]
      break
    }
  }
}
if (is.null(cfg_path)) {
  cfg_path <- "3_Analysis/2_DifferetialAnalysis/scripts/config.yaml"
  message(sprintf("No --config provided. Falling back to default: %s", cfg_path))
}

# ---------------------------
# Helpers
# ---------------------------
ensure_dir <- function(p) {
  # 若目录不存在则递归创建
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

now_str <- function() format(Sys.time(), "%Y%m%d-%H%M%S")

write_tsv <- function(df, path) {
  # 统一以制表符分隔输出，便于查看/复用
  fwrite(df, file = path, sep = "\t", quote = FALSE, na = "NA")
}

# define `%||%` early (used in config parsing)
`%||%` <- function(a, b) if (!is.null(a)) a else b

aggregate_counts_by <- function(count_mat, groups) {
  # 伪样本聚合：将同一样本（例如同一时间点）的细胞计数按列求和
  # count_mat: dgCMatrix genes x cells
  # groups: character vector for cells (length = ncol(count_mat))
  stopifnot(ncol(count_mat) == length(groups))
  samples <- unique(groups)
  # Use rowSums for each group (sparse-friendly)
  res_list <- lapply(samples, function(s) {
    idx <- which(groups == s)
    if (length(idx) == 1L) {
      # rowSums 需要矩阵；当仅一列时特殊处理
      as.numeric(count_mat[, idx, drop = FALSE])
    } else {
      Matrix::rowSums(count_mat[, idx, drop = FALSE])
    }
  })
  res <- do.call(cbind, res_list)
  rownames(res) <- rownames(count_mat)
  colnames(res) <- samples
  methods::as(res, "dgCMatrix")
}

# Assign pseudo-replicates within each sample_id group (e.g., timepoint)
assign_pseudoreps <- function(meta_df, sample_col, k = 3, seed = 42) {
  stopifnot(sample_col %in% colnames(meta_df))
  set.seed(seed)
  lab <- as.character(meta_df[[sample_col]])
  reps <- integer(length(lab))
  u <- unique(lab)
  for (s in u) {
    idx <- which(lab == s)
    k_use <- min(k, length(idx))
    ord <- sample(idx, length(idx))
    # round-robin assignment
    rep_ids <- rep(seq_len(k_use), length.out = length(ord))
    reps[ord] <- rep_ids
  }
  paste0(lab, "_rep", reps)
}

# ---------------------------
# Load config
# ---------------------------
cfg <- yaml::read_yaml(cfg_path)

# paths and columns
out_root <- cfg$paths$out_root
seurat_rds <- cfg$paths$seurat_rds

cols <- cfg$columns
sample_col <- cols$sample_id            # 建议使用 timepoint 作为样本ID（没有 orig.ident 时）
ctype_col  <- cols$cell_type            # 旧逻辑中的“细胞类型”列（本任务不使用该列进行分层）
cluster_col <- cols$cluster %||% "seurat_clusters"  # 新增：用于簇分析的列，默认 seurat_clusters

time_col   <- cols$timepoint            # e.g., timepoint
cond_col   <- cols$condition            # e.g., condition（可能缺失）
group_col  <- cols$group %||% NULL      # optional，若缺失则由模板生成

# 若无条件列，建议 group_label_format 仅使用 {timepoint}
group_label_format <- cfg$group_label_format %||% "{timepoint}_{condition}"

# DEG config
deg_cfg <- cfg$deg
lfc_threshold <- as.numeric(deg_cfg$lfc_threshold %||% 0.25)
padj_cutoff   <- as.numeric(deg_cfg$p_adj_cutoff %||% 0.05)
# 分析范围（整体/簇/两者）
scope <- deg_cfg$scope %||% "by_cluster"

# Pseudobulk config
pb_cfg <- deg_cfg$pseudobulk
agg_method <- pb_cfg$aggregation %||% "sum"
min_cells_per_cluster <- as.integer(pb_cfg$min_cells_per_cluster %||% 50)

# 低表达过滤阈值
filt_cfg <- deg_cfg$filter_expr
min_count   <- as.integer(filt_cfg$min_count %||% 10)
min_samples <- as.integer(filt_cfg$min_samples %||% 2)

# 收缩设置
shrink_cfg <- deg_cfg$shrinkage
do_shrink   <- isTRUE(shrink_cfg$enabled)
shrink_type <- shrink_cfg$type %||% "apeglm"

comparisons <- cfg$comparisons

# ---------------------------
# Prepare dirs and logging
# ---------------------------
ensure_dir(file.path(out_root, "data"))
ensure_dir(file.path(out_root, "results", "tables"))
ensure_dir(file.path(out_root, "results", "data"))
ensure_dir(file.path(out_root, "logs", "deg"))
ensure_dir(file.path(out_root, "reports"))

log_path <- file.path(out_root, "logs", "deg", "1_run_deg_analysis.log")
zz <- file(log_path, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")

cat(sprintf("[INFO] %s | Start DEG analysis\n", now_str()))
cat(sprintf("[INFO] Config: %s\n", cfg_path))
cat(sprintf("[INFO] Seurat RDS: %s\n", seurat_rds))
cat(sprintf("[INFO] Analysis scope: %s\n", scope))

# ---------------------------
# Load Seurat object
# ---------------------------
obj <- readRDS(seurat_rds)
assay_names <- tryCatch(names(obj@assays), error = function(e) character(0))
prefer_assays <- unique(c("RNA", "SCT", DefaultAssay(obj)))
assay_use <- NULL
layer_use <- NULL
for (a in prefer_assays) {
  if (!a %in% assay_names) next
  lyr <- tryCatch(Layers(obj[[a]]), error = function(e) character(0))
  if ("counts" %in% lyr) { assay_use <- a; layer_use <- "counts"; break }
}
if (is.null(assay_use)) {
  # fallback: 扫描所有 assay，选择首个包含 counts 层的
  for (a in assay_names) {
    lyr <- tryCatch(Layers(obj[[a]]), error = function(e) character(0))
    if ("counts" %in% lyr) { assay_use <- a; layer_use <- "counts"; break }
  }
}
if (is.null(assay_use)) {
  stop(sprintf("No 'counts' layer found in assays: %s", paste(assay_names, collapse = ", ")))
}
cat(sprintf("[INFO] Using assay '%s' with layer '%s' for counts\n", assay_use, layer_use))

# metadata
meta <- obj@meta.data
meta$cells <- rownames(meta)

# 若配置的 sample_id 列在 meta 中不存在，而 timepoint 存在，则回退用 timepoint 作为样本ID
if (!sample_col %in% colnames(meta) && time_col %in% colnames(meta)) {
  cat(sprintf("[WARN] sample_id column '%s' not found. Falling back to timepoint.\n", sample_col))
  sample_col <- time_col
}

# Build group column if missing
# 说明：当没有条件列时，仅使用时间点生成分组标签；若模板包含{condition}但cond_col不存在，将改为仅用{timepoint}
if (is.null(group_col) || !(group_col %in% colnames(meta))) {
  if (!time_col %in% colnames(meta)) {
    stop("Cannot construct group: timepoint not found in meta.data")
  }
  fmt <- group_label_format
  if (!is.null(cond_col) && cond_col %in% colnames(meta) && grepl("\\{condition\\}", fmt)) {
    fmt <- gsub("\\{timepoint\\}", meta[[time_col]], fmt)
    fmt <- gsub("\\{condition\\}", meta[[cond_col]], fmt)
  } else {
    # 无条件列或模板未包含{condition}：仅用 timepoint 原值作为分组标签
    fmt <- as.character(meta[[time_col]])
  }
  meta[[".__group"]] <- fmt
  group_col_use <- ".__group"
  cat("[INFO] Constructed group column using timepoint (condition missing or ignored)\n")
} else {
  group_col_use <- group_col
  cat(sprintf("[INFO] Using existing group column: %s\n", group_col_use))
}

# Snapshot metadata and cell counts
snapshot_path <- file.path(out_root, "data", "metadata_snapshot.tsv")
write_tsv(as.data.frame(lapply(meta, function(x) if (is.factor(x)) as.character(x) else x)), snapshot_path)
cat(sprintf("[INFO] Wrote metadata snapshot -> %s\n", snapshot_path))

# Count cells per sample(=timepoint) x cluster（按 seurat_clusters 统计，用于簇分析的参考）
if (!all(c(sample_col, cluster_col) %in% colnames(meta))) {
  stop("Missing sample_id (timepoint) or cluster (seurat_clusters) columns in meta.data")
}
cc <- meta %>%
  count(!!as.name(sample_col), !!as.name(cluster_col), name = "n_cells") %>%
  arrange(!!as.name(sample_col), !!as.name(cluster_col))
counts_by_sample_cluster_path <- file.path(out_root, "data", "counts_by_sample_cluster.tsv")
write_tsv(cc, counts_by_sample_cluster_path)
cat(sprintf("[INFO] Wrote counts_by_sample_cluster -> %s\n", counts_by_sample_cluster_path))

# Prepare count matrix (genes x cells, sparse)
counts_mat <- GetAssayData(obj, layer = layer_use, assay = assay_use)
# 对齐计数矩阵与元数据的细胞顺序/集合
cells_counts <- colnames(counts_mat)
cells_meta <- meta$cells
if (!setequal(cells_counts, cells_meta)) {
  cat(sprintf("[WARN] Cells mismatch between counts (%d) and meta (%d). Aligning intersection.\n",
              length(cells_counts), length(cells_meta)))
  common_cells <- intersect(cells_counts, cells_meta)
  counts_mat <- counts_mat[, common_cells, drop = FALSE]
  meta <- meta[match(common_cells, meta$cells), , drop = FALSE]
}
stopifnot(ncol(counts_mat) == nrow(meta))

# ---------------------------
# Iterate comparisons (整体 + 簇分析)
# ---------------------------
skip_records <- list()

for (cmp in comparisons) {
  cmp_id <- cmp$id
  control_lab <- cmp$control
  case_lab    <- cmp$case

  cat(sprintf("\n[INFO] Comparison: %s | control=%s | case=%s\n", cmp_id, control_lab, case_lab))

  # Ensure output dirs
  out_dir_tables <- file.path(out_root, "results", "tables", cmp_id)
  out_dir_data   <- file.path(out_root, "results", "data", cmp_id)
  ensure_dir(out_dir_tables)
  ensure_dir(out_dir_data)

  # Filter meta rows for this comparison（仅保留当前比较的两组细胞）
  sel_meta <- meta %>% filter(.data[[group_col_use]] %in% c(control_lab, case_lab))
  if (nrow(sel_meta) == 0) {
    cat(sprintf("[WARN] No cells found for comparison %s. Skipping.\n", cmp_id))
    next
  }

  # -----------------------
  # 整体分析（不分簇/细胞类型）：将当前比较下的所有细胞按 timepoint 聚合
  # -----------------------
  if (scope %in% c("overall", "both")) {
    cat("[INFO]   Overall analysis (all cells)\n")
    sel_cells_overall <- sel_meta$cells
    mat_overall <- counts_mat[, sel_cells_overall, drop = FALSE]
    meta_overall <- sel_meta[match(colnames(mat_overall), sel_meta$cells), , drop = FALSE]

    # 依据 sample_id (=timepoint) 进行伪样本聚合；若配置了伪重复则在时间点内拆分为K个子样本
    k_reps <- as.integer(pb_cfg$replicates_per_group %||% 0L)
    cell_samples <- if (isTRUE(k_reps > 1)) {
      assign_pseudoreps(meta_overall, sample_col, k = k_reps)
    } else {
      as.character(meta_overall[[sample_col]])
    }
    if (any(is.na(cell_samples))) {
      msg <- sprintf("Overall | %s: NA in sample_id (%s). Skip.", cmp_id, sample_col)
      cat("[WARN]  ", msg, "\n", sep = "")
      skip_records[[length(skip_records) + 1]] <- list(comparison = cmp_id, cluster = "all_cells", reason = msg)
    } else {
      pb_counts <- aggregate_counts_by(mat_overall, groups = cell_samples)
      # 确定每个伪样本的分组标签（control/case）
      sample_groups <- tapply(meta_overall[[group_col_use]], cell_samples, function(x) unique(as.character(x))[1])
      sample_groups <- unname(sample_groups[colnames(pb_counts)])

      # 仅保留属于当前比较的两类样本列
      keep_cols <- which(sample_groups %in% c(control_lab, case_lab))
      pb_counts <- pb_counts[, keep_cols, drop = FALSE]
      sample_groups <- sample_groups[keep_cols]

      # 检查对照与实验是否各至少有 required_per_group 个伪样本列
      required_per_group <- if (isTRUE(k_reps > 1)) 2L else 1L
      n_ctrl <- sum(sample_groups == control_lab)
      n_case <- sum(sample_groups == case_lab)
      if (n_ctrl < required_per_group || n_case < required_per_group) {
        msg <- sprintf("Overall | %s: insufficient sample columns/replicates (need >=%d per group; ctrl=%d, case=%d). Skip.", cmp_id, required_per_group, n_ctrl, n_case)
        cat("[WARN]  ", msg, "\n", sep = "")
        skip_records[[length(skip_records) + 1]] <- list(comparison = cmp_id, cluster = "all_cells", reason = msg)
      } else {
        # 低表达过滤：至少 min_samples 个伪样本中 count>=min_count
        keep_genes <- Matrix::rowSums(pb_counts >= min_count) >= min_samples
        if (sum(keep_genes) < 5) {
          msg <- sprintf("Overall | %s: too few genes passed filter. Skip.", cmp_id)
          cat("[WARN]  ", msg, "\n", sep = "")
          skip_records[[length(skip_records) + 1]] <- list(comparison = cmp_id, cluster = "all_cells", reason = msg)
        } else {
          pb_counts <- pb_counts[keep_genes, , drop = FALSE]
          # 构建 DESeq2 对象并运行
          coldata <- data.frame(
            row.names = colnames(pb_counts),
            group = factor(ifelse(sample_groups == case_lab, "case", "control"), levels = c("control", "case"))
          )
          dds <- DESeqDataSetFromMatrix(countData = as.matrix(pb_counts), colData = coldata, design = ~ group)
          dds <- DESeq(dds, quiet = TRUE)

          # 结果与收缩
          res <- results(dds, contrast = c("group", "case", "control"))
          res <- lfcShrink(dds, coef = "group_case_vs_control", type = shrink_type, res = res, quiet = TRUE)

          res_df <- as.data.frame(res)
          res_df$gene <- rownames(res_df)
          res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]
          res_df <- res_df %>%
            mutate(sig = ifelse(!is.na(padj) & abs(log2FoldChange) >= lfc_threshold & padj <= padj_cutoff, TRUE, FALSE))

          out_res <- file.path(out_dir_tables, "all_cells.DESeq2.results.csv")
          fwrite(res_df, out_res)

          # 另导出未收缩版本
          res_raw <- results(DESeq(dds, quiet = TRUE), contrast = c("group", "case", "control"))
          res_raw_df <- as.data.frame(res_raw)
          res_raw_df$gene <- rownames(res_raw_df)
          res_raw_df <- res_raw_df[, c("gene", setdiff(colnames(res_raw_df), "gene"))]
          out_res_raw <- file.path(out_dir_tables, "all_cells.DESeq2.results.raw.csv")
          fwrite(res_raw_df, out_res_raw)

          # 标准化计数
          ncounts <- counts(dds, normalized = TRUE)
          ncounts_df <- as.data.frame(ncounts)
          ncounts_df$gene <- rownames(ncounts_df)
          ncounts_df <- ncounts_df[, c("gene", setdiff(colnames(ncounts_df), "gene"))]
          out_nc <- file.path(out_dir_data, "all_cells.normalized_counts.csv")
          fwrite(ncounts_df, out_nc)

          # 排名分数
          if (!is.null(res_df$stat)) {
            ranked <- res_df[, c("gene", "stat")]
          } else if (!is.null(res_raw_df$stat)) {
            ranked <- res_raw_df[, c("gene", "stat")]
          } else {
            sc <- with(res_df, sign(log2FoldChange) * -log10(pvalue))
            ranked <- data.frame(gene = res_df$gene, stat = sc)
          }
          ranked <- ranked[!is.na(ranked$stat), ]
          out_ranked <- file.path(out_dir_data, "all_cells.ranked_genes.csv")
          fwrite(ranked, out_ranked)

          cat(sprintf("[INFO]   Saved overall: %s | %s | %s\n", out_res, out_nc, out_ranked))
        }
      }
    }
  }

  # -----------------------
  # 簇分析（按 seurat_clusters 分层）：不使用 SingleR 的细胞类型
  # -----------------------
  if (scope %in% c("by_cluster", "both")) {
    all_clusters <- unique(meta[[cluster_col]])
    for (cl in all_clusters) {
      cat(sprintf("[INFO]   Cluster: %s\n", cl))
      sel_cells <- sel_meta$cells[sel_meta[[cluster_col]] == cl]
      if (length(sel_cells) < min_cells_per_cluster) {
        msg <- sprintf("Cluster %s | %s: insufficient cells (%d < %d). Skip.",
                       cl, cmp_id, length(sel_cells), min_cells_per_cluster)
        cat("[WARN]  ", msg, "\n", sep = "")
        skip_records[[length(skip_records) + 1]] <- list(comparison = cmp_id, cluster = cl, reason = msg)
        next
      }

      # 子集计数与元数据
      mat_sub <- counts_mat[, sel_cells, drop = FALSE]
      meta_sub <- sel_meta[match(colnames(mat_sub), sel_meta$cells), , drop = FALSE]

      # 伪样本聚合（按 sample_id/timepoint）；若配置了伪重复则拆分为K个子样本
      k_reps <- as.integer(pb_cfg$replicates_per_group %||% 0L)
      cell_samples <- if (isTRUE(k_reps > 1)) {
        assign_pseudoreps(meta_sub, sample_col, k = k_reps)
      } else {
        as.character(meta_sub[[sample_col]])
      }
      if (any(is.na(cell_samples))) {
        msg <- sprintf("Cluster %s | %s: NA in sample_id (%s). Skip.",
                       cl, cmp_id, sample_col)
        cat("[WARN]  ", msg, "\n", sep = "")
        skip_records[[length(skip_records) + 1]] <- list(comparison = cmp_id, cluster = cl, reason = msg)
        next
      }
      pb_counts <- aggregate_counts_by(mat_sub, groups = cell_samples)
      # 每个伪样本的分组标签
      sample_groups <- tapply(meta_sub[[group_col_use]], cell_samples, function(x) unique(as.character(x))[1])
      sample_groups <- unname(sample_groups[colnames(pb_counts)])

      # 仅保留当前比较的两类样本列
      keep_cols <- which(sample_groups %in% c(control_lab, case_lab))
      pb_counts <- pb_counts[, keep_cols, drop = FALSE]
      sample_groups <- sample_groups[keep_cols]

      # 检查对照与实验是否各至少有1个伪样本列
      required_per_group <- if (isTRUE(k_reps > 1)) 2L else 1L
      n_ctrl <- sum(sample_groups == control_lab)
      n_case <- sum(sample_groups == case_lab)
      if (n_ctrl < required_per_group || n_case < required_per_group) {
        msg <- sprintf("Cluster %s | %s: insufficient sample columns/replicates (need >=%d per group; ctrl=%d, case=%d). Skip.",
                       cl, cmp_id, required_per_group, n_ctrl, n_case)
        cat("[WARN]  ", msg, "\n", sep = "")
        skip_records[[length(skip_records) + 1]] <- list(comparison = cmp_id, cluster = cl, reason = msg)
        next
      }

      # 低表达过滤
      keep_genes <- Matrix::rowSums(pb_counts >= min_count) >= min_samples
      if (sum(keep_genes) < 5) {
        msg <- sprintf("Cluster %s | %s: too few genes passed filter. Skip.", cl, cmp_id)
        cat("[WARN]  ", msg, "\n", sep = "")
        skip_records[[length(skip_records) + 1]] <- list(comparison = cmp_id, cluster = cl, reason = msg)
        next
      }
      pb_counts <- pb_counts[keep_genes, , drop = FALSE]

      # 构建 DESeq2 数据集
      coldata <- data.frame(
        row.names = colnames(pb_counts),
        group = factor(ifelse(sample_groups == case_lab, "case", "control"), levels = c("control", "case"))
      )
      dds <- DESeqDataSetFromMatrix(countData = as.matrix(pb_counts), colData = coldata, design = ~ group)

      # 运行 DESeq2
      dds <- DESeq(dds, quiet = TRUE)

      # 结果与收缩
      res <- results(dds, contrast = c("group", "case", "control"))
      res <- lfcShrink(dds, coef = "group_case_vs_control", type = shrink_type, res = res, quiet = TRUE)

      # 输出表
      res_df <- as.data.frame(res)
      res_df$gene <- rownames(res_df)
      res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]
      res_df <- res_df %>%
        mutate(sig = ifelse(!is.na(padj) & abs(log2FoldChange) >= lfc_threshold & padj <= padj_cutoff, TRUE, FALSE))

      out_res <- file.path(out_dir_tables, sprintf("%s.DESeq2.results.csv", cl))
      fwrite(res_df, out_res)

      # 非收缩版本
      res_raw <- results(DESeq(dds, quiet = TRUE), contrast = c("group", "case", "control"))
      res_raw_df <- as.data.frame(res_raw)
      res_raw_df$gene <- rownames(res_raw_df)
      res_raw_df <- res_raw_df[, c("gene", setdiff(colnames(res_raw_df), "gene"))]
      out_res_raw <- file.path(out_dir_tables, sprintf("%s.DESeq2.results.raw.csv", cl))
      fwrite(res_raw_df, out_res_raw)

      # 标准化计数
      ncounts <- counts(dds, normalized = TRUE)
      ncounts_df <- as.data.frame(ncounts)
      ncounts_df$gene <- rownames(ncounts_df)
      ncounts_df <- ncounts_df[, c("gene", setdiff(colnames(ncounts_df), "gene"))]
      out_nc <- file.path(out_dir_data, sprintf("%s.normalized_counts.csv", cl))
      fwrite(ncounts_df, out_nc)

      # 排名分数
      if (!is.null(res_df$stat)) {
        ranked <- res_df[, c("gene", "stat")]
      } else if (!is.null(res_raw_df$stat)) {
        ranked <- res_raw_df[, c("gene", "stat")]
      } else {
        sc <- with(res_df, sign(log2FoldChange) * -log10(pvalue))
        ranked <- data.frame(gene = res_df$gene, stat = sc)
      }
      ranked <- ranked[!is.na(ranked$stat), ]
      out_ranked <- file.path(out_dir_data, sprintf("%s.ranked_genes.csv", cl))
      fwrite(ranked, out_ranked)

      cat(sprintf("[INFO]   Saved cluster %s: %s | %s | %s\n", cl, out_res, out_nc, out_ranked))
    }
  }
}

# Write QC summary for skipped items
if (length(skip_records) > 0) {
  skip_df <- rbindlist(lapply(skip_records, as.data.frame))
  qc_path <- file.path(out_root, "data", "cell_counts_threshold_failures.tsv")
  write_tsv(skip_df, qc_path)
  cat(sprintf("[INFO] Wrote skip/failure summary -> %s\n", qc_path))
}

# Session info
sess_path <- file.path(out_root, "logs", sprintf("sessionInfo_%s.txt", now_str()))
writeLines(c(capture.output(sessionInfo())), con = sess_path)
cat(sprintf("[INFO] Wrote sessionInfo -> %s\n", sess_path))

cat(sprintf("[INFO] %s | DEG analysis finished\n", now_str()))

# close log sinks
sink(type = "message")
sink(type = "output")
close(zz)

# Return invisibly
invisible(TRUE)
