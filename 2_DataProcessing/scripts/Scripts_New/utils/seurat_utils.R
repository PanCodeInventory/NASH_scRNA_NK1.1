#!/usr/bin/env Rscript

#' Seurat对象操作工具函数
#' 用于scRNA-seq数据处理流程中的Seurat对象操作

# 加载必要的包
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(yaml)
  library(stringr)
})

#' 安全读取10x数据
#' @param data_dir 10x数据目录路径
#' @return 表达矩阵
read_10x_data_safe <- function(data_dir) {
  features_path <- file.path(data_dir, "features.tsv.gz")
  genes_path <- file.path(data_dir, "genes.tsv.gz")
  
  # 确定使用哪个文件
  target_path <- if (file.exists(features_path)) features_path else genes_path
  if (!file.exists(target_path)) {
    stop(sprintf("目录缺少features.tsv.gz/genes.tsv.gz：%s", data_dir))
  }
  
  # 检测列数
  first_line <- tryCatch({
    con <- gzfile(target_path, "rt")
    on.exit(close(con), add = TRUE)
    readLines(con, n = 1)
  }, error = function(e) "")
  
  ncols <- if (length(first_line) == 1 && nzchar(first_line)) {
    length(strsplit(first_line, "\t", fixed = TRUE)[[1]])
  } else {
    1
  }
  
  gene_col <- if (ncols >= 2) 2 else 1
  message(sprintf("Read10X: %s using gene.column=%d (detected %d columns)", 
                basename(target_path), gene_col, ncols))
  
  return(Read10X(data.dir = data_dir, gene.column = gene_col))
}

#' 创建Seurat对象并添加元数据
#' @param counts 表达矩阵
#' @param project 项目名称
#' @param timepoint 时间点
#' @return Seurat对象
create_seurat_with_metadata <- function(counts, project, timepoint) {
  obj <- CreateSeuratObject(counts = counts, project = project, min.cells = 0, min.features = 0)
  obj$timepoint <- timepoint
  return(obj)
}

#' SCTransform归一化
#' @param obj Seurat对象
#' @param vst_flavor VST版本
#' @return 归一化后的Seurat对象
normalize_with_sctransform <- function(obj, vst_flavor = "v2") {
  message("Running SCTransform normalization...")
  obj <- SCTransform(obj, vst.flavor = vst_flavor, verbose = FALSE)
  return(obj)
}

#' 数据整合
#' @param obj_list Seurat对象列表
#' @param nfeatures 整合特征数
#' @return 整合后的Seurat对象
integrate_seurat_objects <- function(obj_list, nfeatures = 3000) {
  message("Integrating Seurat objects...")
  
  # 选择整合特征
  features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures)
  
  # 准备整合
  obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features, verbose = FALSE)
  
  # 查找锚点
  anchors <- FindIntegrationAnchors(object.list = obj_list, 
                                  normalization.method = "SCT",
                                  anchor.features = features, 
                                  verbose = FALSE)
  
  # 整合数据
  integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
  DefaultAssay(integrated) <- "integrated"
  
  return(integrated)
}

#' 运行PCA、邻居图、聚类和UMAP
#' @param obj Seurat对象
#' @param dims 主成分数量
#' @param resolution 聚类分辨率
#' @return 处理后的Seurat对象
run_dimensionality_reduction <- function(obj, dims = 30, resolution = 0.5) {
  message(sprintf("Running PCA (dims=1:%d), clustering (res=%.3f), and UMAP...", dims, resolution))
  
  # PCA
  obj <- RunPCA(obj, verbose = FALSE, npcs = dims)
  
  # 邻居图
  obj <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
  
  # 聚类
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  
  # UMAP
  obj <- RunUMAP(obj, dims = 1:dims, verbose = FALSE)
  
  return(obj)
}

#' scDblFinder双胞检测
#' @param obj Seurat对象
#' @param sample_col 样本列名
#' @param expected 期望双胞比例
#' @return 添加双胞标签的Seurat对象
detect_doublets_scDblFinder <- function(obj, sample_col = NULL, expected = 0.06) {
  message("Running scDblFinder for doublet detection...")
  
  # 转换为SingleCellExperiment对象
  sce <- as.SingleCellExperiment(obj)
  
  # 运行scDblFinder
  if (!is.null(sample_col) && sample_col %in% colnames(colData(sce))) {
    sce <- scDblFinder::scDblFinder(sce, samples = colData(sce)[[sample_col]], 
                                   expected = expected)
  } else {
    sce <- scDblFinder::scDblFinder(sce, expected = expected)
  }
  
  # 将结果写回Seurat对象
  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score
  
  return(obj)
}

#' SingleR细胞类型注释
#' @param obj Seurat对象
#' @param ref_dataset 参考数据集名称
#' @return 添加注释的Seurat对象
annotate_cells_singleR <- function(obj, ref_dataset = "ImmGenData") {
  message("Running SingleR cell type annotation...")
  
  # 获取参考数据
  ref <- get(ref_dataset)()
  
  # 提取counts数据
  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  
  # 创建SingleCellExperiment对象
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  sce <- scater::logNormCounts(sce)
  
  # 运行SingleR
  pred <- SingleR::SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.main,
    assay.type.test = "logcounts",
    assay.type.ref = "logcounts"
  )
  
  # 添加注释到Seurat对象
  obj$singleR.label <- pred$labels
  if (!is.null(pred$pruned.labels)) {
    obj$singleR.pruned <- pred$pruned.labels
  }
  
  return(obj)
}

#' 基于细胞类型过滤
#' @param obj Seurat对象
#' @param remove_types 要去除的细胞类型列表
#' @param keep_types 要保留的细胞类型列表
#' @param label_col 标签列名
#' @return 过滤后的Seurat对象
filter_cells_by_type <- function(obj, remove_types = NULL, keep_types = NULL, 
                              label_col = "singleR.pruned") {
  message("Filtering cells by cell type...")
  
  if (is.null(label_col) || !(label_col %in% colnames(obj@meta.data))) {
    stop("Label column not found in Seurat object")
  }
  
  labels <- tolower(as.character(obj@meta.data[[label_col]]))
  
  # 创建过滤向量
  keep_cells <- rep(TRUE, length(labels))
  
  # 处理去除规则
  if (!is.null(remove_types)) {
    for (type in remove_types) {
      pattern <- tolower(type)
      keep_cells <- keep_cells & !grepl(pattern, labels, ignore.case = TRUE)
    }
  }
  
  # 处理保留规则（优先级更高）
  if (!is.null(keep_types)) {
    keep_vector <- rep(FALSE, length(labels))
    for (type in keep_types) {
      pattern <- tolower(type)
      keep_vector <- keep_vector | grepl(pattern, labels, ignore.case = TRUE)
    }
    keep_cells <- keep_vector
  }
  
  # 应用过滤
  obj_filtered <- subset(obj, cells = colnames(obj)[keep_cells])
  
  message(sprintf("Filtered from %d to %d cells (%.1f%% removed)", 
                ncol(obj), ncol(obj_filtered), 
                100 * (1 - ncol(obj_filtered) / ncol(obj))))
  
  return(obj_filtered)
}

#' 参数调优网格搜索
#' @param obj Seurat对象
#' @param dims_grid 主成分候选
#' @param res_grid 分辨率候选
#' @return 包含评估指标的数据框
tune_clustering_parameters <- function(obj, dims_grid = c(10, 15, 20, 25, 30),
                                  res_grid = seq(0.2, 1.2, by = 0.1)) {
  message("Running parameter tuning grid search...")
  
  metrics_list <- list()
  has_cluster_pkg <- requireNamespace("cluster", quietly = TRUE)
  
  for (dims in dims_grid) {
    message(sprintf("[dims=%d] Running FindNeighbors/RunUMAP", dims))
    obj_temp <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
    obj_temp <- RunUMAP(obj_temp, dims = 1:dims, verbose = FALSE)
    
    for (res in res_grid) {
      message(sprintf("  └─ FindClusters(res=%.1f)", res))
      obj_temp <- FindClusters(obj_temp, resolution = res, verbose = FALSE)
      
      # 计算指标
      cluster_col <- paste0("integrated_snn_res.", as.character(res))
      if (!(cluster_col %in% colnames(obj_temp@meta.data))) {
        warning(sprintf("Cluster column %s not found", cluster_col))
        next
      }
      
      comp <- obj_temp@meta.data %>%
        dplyr::count(!!rlang::sym(cluster_col), name = "n") %>%
        dplyr::mutate(frac = n / sum(n))
      
      n_clusters <- nrow(comp)
      min_frac <- if (n_clusters > 0) min(comp$frac) else NA_real_
      
      # 轮廓系数
      median_sil <- NA_real_
      if (has_cluster_pkg && n_clusters >= 2) {
        tryCatch({
          emb <- Embeddings(obj_temp, "pca")[, 1:dims, drop = FALSE]
          cls <- as.integer(factor(obj_temp@meta.data[[cluster_col]]))
          dist_mat <- stats::dist(emb)
          sil <- cluster::silhouette(cls, dist_mat)
          median_sil <- stats::median(sil[, "sil_width"], na.rm = TRUE)
        }, error = function(e) {
          message(sprintf("Silhouette calculation failed: %s", e$message))
        })
      }
      
      metrics_list[[length(metrics_list) + 1]] <- data.frame(
        dims = dims,
        res = res,
        n_clusters = n_clusters,
        min_cluster_frac = min_frac,
        median_silhouette = median_sil,
        stringsAsFactors = FALSE
      )
    }
  }
  
  return(dplyr::bind_rows(metrics_list))
}

#' 应用簇重命名
#' @param obj Seurat对象
#' @param cluster_mapping 簇映射配置
#' @return 重命名后的Seurat对象
apply_cluster_renaming <- function(obj, cluster_mapping) {
  message("Applying cluster renaming...")
  
  original_clusters <- as.character(obj$seurat_clusters)
  renamed_clusters <- original_clusters
  
  # 应用映射
  for (mapping in cluster_mapping) {
    old_name <- as.character(mapping$old_cluster)
    new_name <- as.character(mapping$new_cluster)
    
    renamed_clusters[original_clusters == old_name] <- new_name
    message(sprintf("Cluster %s → %s", old_name, new_name))
  }
  
  # 更新对象
  obj$seurat_clusters <- factor(renamed_clusters)
  Idents(obj) <- factor(renamed_clusters)
  obj$original_cluster <- original_clusters
  
  return(obj)
}

#' 验证Seurat对象
#' @param obj Seurat对象
#' @param check_items 检查项目列表
#' @return 验证结果列表
validate_seurat_object <- function(obj, check_items = c("cells", "genes", "assays", "metadata")) {
  validation_results <- list()
  
  if ("cells" %in% check_items) {
    validation_results$cells <- list(
      count = ncol(obj),
      status = if (ncol(obj) > 0) "PASS" else "FAIL"
    )
  }
  
  if ("genes" %in% check_items) {
    validation_results$genes <- list(
      count = nrow(obj),
      status = if (nrow(obj) > 0) "PASS" else "FAIL"
    )
  }
  
  if ("assays" %in% check_items) {
    assays <- Assays(obj)
    validation_results$assays <- list(
      available = assays,
      count = length(assays),
      status = if (length(assays) > 0) "PASS" else "FAIL"
    )
  }
  
  if ("metadata" %in% check_items) {
    metadata_cols <- colnames(obj@meta.data)
    validation_results$metadata <- list(
      columns = metadata_cols,
      count = length(metadata_cols),
      status = if (length(metadata_cols) > 0) "PASS" else "FAIL"
    )
  }
  
  return(validation_results)
}

#' 加载配置文件
#' @param config_path 配置文件路径
#' @return 配置列表
load_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop(sprintf("Configuration file not found: %s", config_path))
  }
  
  config <- yaml::read_yaml(config_path)
  message(sprintf("Loaded configuration from: %s", config_path))
  return(config)
}

message("Seurat utility functions loaded successfully!")