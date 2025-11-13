#!/usr/bin/env Rscript

#' 主控制脚本
#' 
#' 功能：
#' 1. 流程协调
#' 2. 参数传递
#' 3. 进度监控
#' 4. 错误处理
#' 
#' 输入：配置文件路径、可选的阶段覆盖参数
#' 输出：完整的执行日志、各阶段的输出文件、最终状态报告、错误和警告汇总

# 加载必要的包
suppressPackageStartupMessages({
  library(yaml)
  library(optparse)
})

#' 选项配置
option_list <- list(
  make_option(c("-c", "--config"), type = "character", 
              help = "配置文件路径 [默认: config/parameters.yaml]",
              default = "config/parameters.yaml"),
  make_option(c("-p", "--paths"), type = "character", 
              help = "路径配置文件路径 [默认: config/paths.yaml]",
              default = "config/paths.yaml"),
  make_option(c("-m", "--mapping"), type = "character", 
              help = "簇映射配置文件路径 [默认: config/cluster_mapping.yaml]",
              default = "config/cluster_mapping.yaml"),
  make_option(c("-s", "--stages"), type = "character", 
              help = "要执行的阶段 (逗号分隔) [默认: all]",
              default = "all"),
  make_option(c("-r", "--resume"), type = "character", 
              help = "从指定阶段恢复执行 [默认: 从头开始]"),
  make_option(c("-v", "--verbose"), action = "store_true", 
              help = "详细输出 [默认: FALSE]",
              default = FALSE),
  make_option(c("-d", "--dry-run"), action = "store_true", 
              help = "试运行模式，只显示将要执行的操作 [默认: FALSE]",
              default = FALSE)
)

#' 解析命令行参数
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#' 阶段定义
STAGES <- list(
  data_import = list(
    script = "01_data_import.R",
    description = "数据导入与预处理",
    input = "10x数据目录",
    output = "nk.integrated.imported.rds"
  ),
  basic_qc = list(
    script = "02_basic_qc.R",
    description = "基础质量控制",
    input = "nk.integrated.imported.rds",
    output = "nk.integrated.qc_filtered.rds"
  ),
  cell_annotation = list(
    script = "03_cell_annotation.R",
    description = "细胞标注",
    input = "nk.integrated.qc_filtered.rds",
    output = "nk.integrated.annotated.rds"
  ),
  cell_filtering = list(
    script = "04_cell_filtering.R",
    description = "细胞过滤",
    input = "nk.integrated.annotated.rds",
    output = "nk.integrated.filtered.rds"
  ),
  parameter_tuning = list(
    script = "05_parameter_tuning.R",
    description = "参数优化",
    input = "nk.integrated.filtered.rds",
    output = "nk.integrated.tuned.rds"
  ),
  visualization = list(
    script = "06_visualization.R",
    description = "可视化",
    input = "nk.integrated.tuned.rds",
    output = "nk.integrated.final.rds"
  )
)

#' 日志记录函数
log_message <- function(message, level = "INFO", timestamp = TRUE) {
  if (timestamp) {
    msg <- sprintf("[%s] %s: %s", 
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, message)
  } else {
    msg <- sprintf("%s: %s", level, message)
  }
  
  cat(msg, "\n")
  flush.console()
  
  # 写入日志文件
  if (exists("log_file") && !is.null(log_file)) {
    cat(msg, "\n", file = log_file, append = TRUE)
  }
}

#' 执行阶段脚本
execute_stage <- function(stage_name, stage_config, config_path, paths_path, mapping_path) {
  log_message(sprintf("开始执行阶段: %s", stage_config$description), "INFO")
  
  script_path <- stage_config$script
  
  if (opt$dry_run) {
    log_message(sprintf("[试运行] Rscript %s %s %s %s", 
                    script_path, config_path, paths_path, mapping_path), "INFO")
    return(TRUE)
  }
  
  # 构建命令
  cmd <- sprintf("Rscript %s %s %s %s", 
                script_path, config_path, paths_path, mapping_path)
  
  log_message(sprintf("执行命令: %s", cmd), "DEBUG")
  
  # 执行脚本
  result <- system(cmd, intern = TRUE)
  
  # 检查执行结果
  exit_code <- attr(result, "status")
  if (is.null(exit_code)) exit_code <- 0
  
  if (exit_code == 0) {
    log_message(sprintf("阶段 %s 执行成功", stage_config$description), "INFO")
    return(TRUE)
  } else {
    log_message(sprintf("阶段 %s 执行失败 (退出码: %d)", 
                    stage_config$description, exit_code), "ERROR")
    if (length(result) > 0) {
      log_message("错误输出:", "ERROR")
      for (line in result) {
        log_message(sprintf("  %s", line), "ERROR")
      }
    }
    return(FALSE)
  }
}

#' 检查阶段依赖
check_dependencies <- function(stage_name, completed_stages) {
  # 定义阶段依赖关系
  dependencies <- list(
    data_import = c(),
    basic_qc = c("data_import"),
    cell_annotation = c("basic_qc"),
    cell_filtering = c("cell_annotation"),
    parameter_tuning = c("cell_filtering"),
    visualization = c("parameter_tuning")
  )
  
  stage_deps <- dependencies[[stage_name]]
  if (is.null(stage_deps)) return(TRUE)
  
  for (dep in stage_deps) {
    if (!(dep %in% completed_stages)) {
      log_message(sprintf("阶段 %s 依赖 %s，但后者未完成", stage_name, dep), "ERROR")
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#' 主函数
main <- function() {
  # 创建日志文件
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file <- sprintf("pipeline_log_%s.txt", timestamp)
  
  log_message("=== scRNA-seq 数据处理流程开始 ===", "INFO")
  log_message(sprintf("日志文件: %s", log_file), "INFO")
  log_message(sprintf("配置文件: %s", opt$config), "INFO")
  log_message(sprintf("路径配置: %s", opt$paths), "INFO")
  log_message(sprintf("簇映射配置: %s", opt$mapping), "INFO")
  
  # 确定要执行的阶段
  if (opt$stages == "all") {
    stages_to_run <- names(STAGES)
  } else {
    stages_to_run <- strsplit(opt$stages, ",")[[1]]
    stages_to_run <- trimws(stages_to_run)
    
    # 验证阶段名称
    invalid_stages <- setdiff(stages_to_run, names(STAGES))
    if (length(invalid_stages) > 0) {
      log_message(sprintf("无效的阶段名称: %s", paste(invalid_stages, collapse = ", ")), "ERROR")
      quit(status = 1)
    }
  }
  
  log_message(sprintf("将执行阶段: %s", paste(stages_to_run, collapse = ", ")), "INFO")
  
  # 确定起始阶段
  if (!is.null(opt$resume)) {
    start_index <- which(names(STAGES) == opt$resume)
    if (length(start_index) == 0) {
      log_message(sprintf("无效的恢复阶段: %s", opt$resume), "ERROR")
      quit(status = 1)
    }
    stages_to_run <- names(STAGES)[start_index:length(STAGES)]
    log_message(sprintf("从阶段 %s 恢复执行", opt$resume), "INFO")
  }
  
  # 执行阶段
  completed_stages <- c()
  failed_stages <- c()
  
  for (stage_name in stages_to_run) {
    stage_config <- STAGES[[stage_name]]
    
    # 检查依赖
    if (!check_dependencies(stage_name, completed_stages)) {
      failed_stages <- c(failed_stages, stage_name)
      next
    }
    
    # 执行阶段
    success <- execute_stage(stage_name, stage_config, opt$config, opt$paths, opt$mapping)
    
    if (success) {
      completed_stages <- c(completed_stages, stage_name)
    } else {
      failed_stages <- c(failed_stages, stage_name)
      
      # 询问是否继续
      if (!opt$dry_run) {
        log_message("阶段执行失败，是否继续执行后续阶段？(y/N)", "QUESTION")
        response <- readline(prompt = ": ")
        if (tolower(response) != "y") {
          log_message("用户选择停止执行", "INFO")
          break
        }
      }
    }
  }
  
  # 生成执行报告
  log_message("=== 执行报告 ===", "INFO")
  log_message(sprintf("成功阶段: %s", paste(completed_stages, collapse = ", ")), "INFO")
  log_message(sprintf("失败阶段: %s", paste(failed_stages, collapse = ", ")), "INFO")
  log_message(sprintf("总阶段数: %d", length(stages_to_run)), "INFO")
  log_message(sprintf("成功率: %.1f%%", 100 * length(completed_stages) / length(stages_to_run)), "INFO")
  
  # 保存执行报告
  report_path <- sprintf("pipeline_report_%s.txt", timestamp)
  report_lines <- c(
    "# scRNA-seq 数据处理流程执行报告",
    "",
    sprintf("执行时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("日志文件: %s", log_file),
    "",
    "## 配置参数",
    sprintf("- 配置文件: %s", opt$config),
    sprintf("- 路径配置: %s", opt$paths),
    sprintf("- 簇映射配置: %s", opt$mapping),
    sprintf("- 详细输出: %s", if (opt$verbose) "是" else "否"),
    sprintf("- 试运行模式: %s", if (opt$dry_run) "是" else "否"),
    "",
    "## 执行阶段",
    sprintf("- 计划执行: %s", paste(stages_to_run, collapse = ", ")),
    sprintf("- 成功完成: %s", paste(completed_stages, collapse = ", ")),
    sprintf("- 执行失败: %s", paste(failed_stages, collapse = ", ")),
    sprintf("- 成功率: %.1f%%", 100 * length(completed_stages) / length(stages_to_run)),
    "",
    "## 输出文件",
    "各阶段生成的文件请参考对应的阶段报告",
    "",
    "## 建议",
    if (length(failed_stages) > 0) {
      "检查失败的阶段日志，修复问题后重新运行"
    } else {
      "所有阶段执行成功，可以进行下游分析"
    }
  )
  
  writeLines(report_lines, report_path)
  log_message(sprintf("执行报告保存到: %s", report_path), "INFO")
  
  # 最终状态
  if (length(failed_stages) == 0) {
    log_message("=== 流程执行完成 ===", "INFO")
    log_message("所有阶段执行成功!", "INFO")
    quit(status = 0)
  } else {
    log_message("=== 流程执行中断 ===", "ERROR")
    log_message(sprintf("有 %d 个阶段执行失败", length(failed_stages)), "ERROR")
    quit(status = 1)
  }
}

# 错误处理
tryCatch({
  main()
}, error = function(e) {
  log_message(sprintf("致命错误: %s", e$message), "ERROR")
  quit(status = 1)
})