#!/bin/bash
# 初始化 micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate r-grab-env

Rscript - <<'EOF'
library(GRAB)

# 输出目录
OutputDir <- "/home/user2/jzhang_data/EWAS/ukb_ckm/results/assoc"
if (!dir.exists(OutputDir)) {
  dir.create(OutputDir, recursive = TRUE)
}

# 分染色体循环
for (chr in 1:22) {
  message(">>> Running association on chromosome ", chr, " ...")
  
  # 载入对应染色体的 Null Model
  NullModelFile <- paste0("/home/user2/jzhang_data/EWAS/ukb_ckm/results/objPOLMM_chr", chr, "_sparseGRM.RData")
  load(NullModelFile)   # 会加载成对象 obj.POLMM
  
  # 对应染色体的基因型文件
  GenoFile <- paste0("/home/user2/jzhang_data/ukb_wes_plink/ukb_wes_eur_chr", chr, ".bed")
  
  # 输出结果文件
  OutputFile <- file.path(OutputDir, paste0("variant_assoc_chr", chr, ".txt"))
  
  # 运行单变量检验
  GRAB.Marker(
    obj.POLMM,
    GenoFile   = GenoFile,
    OutputFile = OutputFile
  )
}
EOF