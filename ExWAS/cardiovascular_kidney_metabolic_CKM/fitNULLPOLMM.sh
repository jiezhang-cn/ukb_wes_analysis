#!/bin/bash
# 初始化 micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate r-grab-env

Rscript - <<'EOF'
library(GRAB)
library(data.table)

# phenotype data
ukb_CKM <- fread("/home/user2/jzhang_data/EWAS/ukb_ckm/data/ukb_CKM", header = TRUE)
ukb_CKM$ckm_stage <- factor(ukb_CKM$ckm_stage, levels = c(0, 1, 2, 3, 4))

# 要标准化的协变量列表
covariates <- c("age", "bmi", paste0("pc", 1:10))

# 对这些协变量做 z-score 标准化
ukb_CKM[, (covariates) := lapply(.SD, scale), .SDcols = covariates]

# 全基因组 Sparse GRM 文件
SparseGRMFile <- "/home/user2/jzhang_data/EWAS/ukb_ckm/results/SparseGRM_allchr.txt"

# 分染色体循环
for (chr in 1:22) {
  message("Processing chromosome ", chr, " ...")
  
  # 每条染色体的 plink bed 文件
  GenoFile <- paste0("/home/user2/jzhang_data/ukb_wes_plink/ukb_wes_eur_chr", chr, ".bed")
  
  obj.POLMM <- GRAB.NullModel(
    formula = ckm_stage ~ age + sex + bmi + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
    data = ukb_CKM,
    subjData = ukb_CKM$IID,
    method = "POLMM",
    traitType = "ordinal",
    GenoFile = GenoFile,
    SparseGRMFile = SparseGRMFile,   
    control = list(showInfo = TRUE, AlleleOrder = "alt-first", LOCO = FALSE, tolTau = 0.2, tolBeta = 0.1)
  )
  
  save(obj.POLMM,
       file = paste0("/home/user2/jzhang_data/EWAS/ukb_ckm/results/objPOLMM_chr", chr, "_sparseGRM.RData"))
}
EOF