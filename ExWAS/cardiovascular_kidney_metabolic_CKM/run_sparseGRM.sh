#!/bin/bash

# 初始化 micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate r-grab-env

Rscript - <<'EOF'
library(GRAB)

# gcta64 路径
gcta64File <- "/home/user2/jzhang_data/EWAS/ukb_ckm/codes/gcta-1.94.1-linux-kernel-4-x86_64/gcta64"

# UKB 建议 nPartsGRM=250
nPartsGRM <- 250  

# 输出目录
outDir <- "/home/user2/jzhang_data/EWAS/ukb_ckm/results"

# 指定中间结果目录
tempDirPath <- "/home/user2/jzhang_data/EWAS/ukb_ckm/codes/temp"

# 如果目录不存在，先创建
if (!dir.exists(tempDirPath)) {
  dir.create(tempDirPath, recursive = TRUE)
}

# 全基因组合并后的 PLINK 文件前缀（不带 .bed/.bim/.fam 后缀）
PlinkFile <- "/home/user2/jzhang_data/ukb_wes_qc_data/relationship_ancestry/ukb_wes_chr_all_king_sample_qc_final_unrelated"

# 输出文件
SparseGRMFile <- file.path(outDir, "SparseGRM_allchr.txt")

message(">>> Step 1: 生成全基因组 GRM 分块文件")

for (part in 1:nPartsGRM) {
  message("  - 全基因组 Part ", part, " ...")
  getTempFilesFullGRM(
    PlinkFile    = PlinkFile,
    nPartsGRM    = nPartsGRM,
    partParallel = part,
    gcta64File   = gcta64File,
    threadNum    = 36,
    tempDir      = tempDirPath
  )
}

message(">>> Step 2: 合并为全基因组 Sparse GRM")

getSparseGRM(
  PlinkFile        = PlinkFile,
  nPartsGRM        = nPartsGRM,
  SparseGRMFile    = SparseGRMFile,
  tempDir          = tempDirPath,
  relatednessCutoff = 0.05,   
  minMafGRM        = 0.01,
  maxMissingGRM    = 0.1,
  rm.tempFiles     = TRUE
)

message(">>> 全基因组 Sparse GRM 已完成: ", SparseGRMFile)

EOF