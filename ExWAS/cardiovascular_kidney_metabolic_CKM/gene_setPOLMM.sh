#!/bin/bash
# 初始化 micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate r-grab-env

Rscript - <<'EOF'
library(GRAB)

# ---- 覆写 GRAB 包命名空间里的 CCT 函数 ----
fixCCT <- function() {
  assignInNamespace("CCT", function(pvals, weights = NULL) {
    if (sum(is.na(pvals)) > 0) {
      return(NA)   # 修改：遇到 NA 返回 NA
    }

    if ((sum(pvals < 0) + sum(pvals > 1)) > 0) {
      stop("All p-values must be between 0 and 1!")
    }

    is.zero <- (sum(pvals == 0) >= 1)
    which.one <- which(pvals == 1)
    is.one <- (length(which.one) >= 1)
    if (is.zero && is.one) stop("Cannot have both 0 and 1 p-values!")
    if (is.zero) return(0)
    if (is.one) {
      warning("There are p-values that are exactly 1!")
      pvals[which.one] <- 0.999
    }

    if (is.null(weights)) {
      weights <- rep(1 / length(pvals), length(pvals))
    } else if (length(weights) != length(pvals)) {
      stop("The length of weights should be the same as that of the p-values!")
    } else if (sum(weights < 0) > 0) {
      stop("All the weights must be positive!")
    } else {
      weights <- weights / sum(weights)
    }

    is.small <- (pvals < 1e-16)
    if (sum(is.small) == 0) {
      cct.stat <- sum(weights * tan((0.5 - pvals) * pi))
    } else {
      cct.stat <- sum((weights[is.small] / pvals[is.small]) / pi)
      cct.stat <- cct.stat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
    }

    if (cct.stat > 1e+15) {
      pval <- (1 / cct.stat) / pi
    } else {
      pval <- 1 - pcauchy(cct.stat)
    }
    return(pval)
  }, ns = "GRAB")
}
fixCCT()

# ---- 主循环 ----

# 全基因组 Sparse GRM 文件
SparseGRMFile <- "/home/user2/jzhang_data/EWAS/ukb_ckm/results/SparseGRM_allchr.txt"

# 输出目录
OutputDir <- "/home/user2/jzhang_data/EWAS/ukb_ckm/results/region_assoc"
if (!dir.exists(OutputDir)) {
  dir.create(OutputDir, recursive = TRUE)
}

# 分染色体循环
for (chr in 1:22) {
  message(">>> Running POLMM-GENE set-based test on chromosome ", chr, " ...")
  
  # 载入对应的 Null Model
  NullModelFile <- paste0("/home/user2/jzhang_data/EWAS/ukb_ckm/results/objPOLMM_chr", chr, "_sparseGRM.RData")
  load(NullModelFile)   # 加载成 obj.POLMM
  
  # 分染色体基因型文件
  GenoFile <- paste0("/home/user2/jzhang_data/ukb_wes_plink/ukb_wes_eur_chr", chr, ".bed")
  
  # 分染色体注释 GroupFile
  GroupFile <- paste0("/home/user2/jzhang_data/EWAS/annotation_files/SnpEff_gene_group_chr", chr, "_nochr_nodup.txt")
  
  # 输出结果文件
  OutputFile <- file.path(OutputDir, paste0("gene_assoc_chr", chr, ".txt"))
  
  # 运行基因集检验
  GRAB.Region(
    objNull        = obj.POLMM,
    GenoFile       = GenoFile,
    GenoFileIndex  = NULL,
    OutputFile     = OutputFile,
    OutputFileIndex = NULL,
    GroupFile      = GroupFile,
    SparseGRMFile  = SparseGRMFile,
    MaxMAFVec      = "0.0001,0.001,0.01",
    annoVec        = "lof,lof:missense"
  )
}
EOF