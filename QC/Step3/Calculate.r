# 加载包
library(data.table)
library(dplyr)

# 清理环境
rm(list = ls())

# 调试信息
cat("当前工作目录:", getwd(), "\n")

# 获取文件路径
file_dir <- "/opt/notebooks/QCStep3"
file_pattern <- "\\.scount$"
cat("查找目录:", file_dir, "\n")

# 检查目录是否存在
if(!dir.exists(file_dir)) {
  stop("目录不存在:", file_dir)
}

# 获取文件列表
files <- list.files(path = file_dir, pattern = file_pattern, full.names = TRUE)
cat("找到文件数:", length(files), "\n")

if(length(files) == 0) {
  stop("未找到匹配的文件")
}

# 显示前几个文件
cat("前几个文件:", paste(head(files), collapse="\n"), "\n")

# 检查第一个文件是否可访问
if(!file.exists(files[1])) {
  stop("文件不存在:", files[1])
}

# 尝试读取第一个文件
tryCatch({
  first_data <- fread(files[1])
  cat("成功读取第一个文件，行数:", nrow(first_data), "列数:", ncol(first_data), "\n")
  cat("列名:", paste(colnames(first_data), collapse=", "), "\n")
}, error = function(e) {
  stop("读取文件失败:", e$message)
})

# 读取所有文件到列表
sc_all <- list()
for(i in 1:length(files)) {
  cat("读取文件:", i, "/", length(files), "\n")
  tryCatch({
    data <- fread(files[i])
    sc_all[[i]] <- data
  }, error = function(e) {
    cat("文件读取错误:", files[i], "-", e$message, "\n")
    # 继续处理下一个文件
  })
}

# 检查是否所有文件都读取成功
if(length(sc_all) == 0) {
  stop("没有成功读取任何文件")
}

# 获取列名
col_names <- colnames(sc_all[[1]])
cat("列名:", paste(col_names, collapse=", "), "\n")

# 验证所有文件都包含相同的样本ID
sample_count <- nrow(sc_all[[1]])
for(i in 2:length(sc_all)) {
  if(nrow(sc_all[[i]]) != sample_count) {
    warning("文件 ", i, " 的样本数 (", nrow(sc_all[[i]]), ") 与第一个文件 (", sample_count, ") 不匹配")
  }
  
  # 检查样本ID是否匹配
  if(!identical(sc_all[[1]][[1]], sc_all[[i]][[1]])) {
    warning("文件 ", i, " 的样本ID与第一个文件不完全匹配")
    
    # 计算匹配率
    match_percent <- sum(sc_all[[1]][[1]] %in% sc_all[[i]][[1]]) / length(sc_all[[1]][[1]]) * 100
    cat("样本ID匹配率: ", round(match_percent, 2), "%\n")
  }
}

# 假设每个染色体文件都包含相同的样本ID
# 使用第一个文件的样本ID和结构创建结果数据框
sc_all_df <- as.data.frame(sc_all[[1]])

# 从第2个文件开始，按行(样本ID)添加统计量
if(length(sc_all) > 1) {
  for(i in 2:length(sc_all)) {
    current_df <- as.data.frame(sc_all[[i]])
    # 确保ID列匹配
    if(!identical(sc_all_df[,1], current_df[,1])) {
      warning("文件 ", i, " 的样本ID与结果数据框不匹配，尝试重新排序")
      current_df <- current_df[match(sc_all_df[,1], current_df[,1]),]
      # 检查是否有NA值(如果有表示某些ID在当前文件中不存在)
      if(any(is.na(current_df[,1]))) {
        warning("某些样本ID在文件 ", i, " 中不存在")
      }
    }
    # 累加数值列(从第2列开始，第1列是ID)
    sc_all_df[,2:ncol(sc_all_df)] <- sc_all_df[,2:ncol(sc_all_df)] + current_df[,2:ncol(current_df)]
  }
}

# 使用汇总结果数据框
sc_all_matrix <- sc_all_df

# 计算质控指标
cat("计算质控指标...\n")
sc_all_matrix$Ti_Tv <- sc_all_matrix$TRANSITION_CT/sc_all_matrix$TRANSVERSION_CT
sc_all_matrix$Het_Hom <- sc_all_matrix$HET_CT/sc_all_matrix$HOM_CT
sc_all_matrix$singletons <- sc_all_matrix$SINGLETON_CT
sc_all_matrix$call_rate <- 1-(sc_all_matrix$MISSING_CT/25520123)

# 检查计算的指标
cat("Ti/Tv 比率统计: 均值=", mean(sc_all_matrix$Ti_Tv, na.rm=TRUE),
    "中位数=", median(sc_all_matrix$Ti_Tv, na.rm=TRUE),
    "SD=", sd(sc_all_matrix$Ti_Tv, na.rm=TRUE), "\n")
cat("Het/Hom 比率统计: 均值=", mean(sc_all_matrix$Het_Hom, na.rm=TRUE),
    "中位数=", median(sc_all_matrix$Het_Hom, na.rm=TRUE),
    "SD=", sd(sc_all_matrix$Het_Hom, na.rm=TRUE), "\n")
cat("Singletons 统计: 均值=", mean(sc_all_matrix$singletons, na.rm=TRUE),
    "中位数=", median(sc_all_matrix$singletons, na.rm=TRUE),
    "SD=", sd(sc_all_matrix$singletons, na.rm=TRUE), "\n")
cat("Call rate 统计: 均值=", mean(sc_all_matrix$call_rate, na.rm=TRUE),
    "中位数=", median(sc_all_matrix$call_rate, na.rm=TRUE),
    "SD=", sd(sc_all_matrix$call_rate, na.rm=TRUE), "\n")

# 计算统计值
Ti_Tv_mean <- mean(sc_all_matrix$Ti_Tv, na.rm=TRUE)
Ti_Tv_sd <- sd(sc_all_matrix$Ti_Tv, na.rm=TRUE)
Het_Hom_mean <- mean(sc_all_matrix$Het_Hom, na.rm=TRUE)
Het_Hom_sd <- sd(sc_all_matrix$Het_Hom, na.rm=TRUE)
singletons_mean <- mean(sc_all_matrix$singletons, na.rm=TRUE)
singletons_sd <- sd(sc_all_matrix$singletons, na.rm=TRUE)

# 过滤样本
Ti_Tv_pass <- subset(sc_all_matrix, (sc_all_matrix$Ti_Tv > (Ti_Tv_mean-8*Ti_Tv_sd)) & 
                      (sc_all_matrix$Ti_Tv < (Ti_Tv_mean+8*Ti_Tv_sd)), select = 1, drop = FALSE)
Het_Hom_pass <- subset(sc_all_matrix, (sc_all_matrix$Het_Hom > (Het_Hom_mean-8*Het_Hom_sd)) & 
                        (sc_all_matrix$Het_Hom < (Het_Hom_mean+8*Het_Hom_sd)), select = 1, drop = FALSE)
singletons_pass <- subset(sc_all_matrix, (sc_all_matrix$singletons > (singletons_mean-8*singletons_sd)) & 
                           (sc_all_matrix$singletons < (singletons_mean+8*singletons_sd)), select = 1, drop = FALSE)
call_rate_90_pass <- subset(sc_all_matrix, sc_all_matrix$call_rate >= 0.9, select = 1, drop = FALSE)

# 合并过滤结果使用基本R函数
id_ti_tv <- Ti_Tv_pass[[1]]
id_het_hom <- Het_Hom_pass[[1]]
id_singletons <- singletons_pass[[1]]
id_call_rate <- call_rate_90_pass[[1]]

# 找到所有条件都满足的ID
common_ids <- Reduce(intersect, list(id_ti_tv, id_het_hom, id_singletons, id_call_rate))
cat("通过所有过滤的样本数:", length(common_ids), "\n")

# 创建最终保留ID表
sample_qc_final_keep <- data.frame(
  FID = common_ids,
  IID = common_ids
)

# 输出统计信息
cat("\n过滤统计：\n")
cat("总样本数:", nrow(sc_all_matrix), "\n")
cat("Ti/Tv 过滤:", nrow(sc_all_matrix) - length(id_ti_tv), "样本 (", 
    round((nrow(sc_all_matrix) - length(id_ti_tv))/nrow(sc_all_matrix)*100, 2), "%)\n")
cat("Het/Hom 过滤:", nrow(sc_all_matrix) - length(id_het_hom), "样本 (",
    round((nrow(sc_all_matrix) - length(id_het_hom))/nrow(sc_all_matrix)*100, 2), "%)\n")
cat("Singletons 过滤:", nrow(sc_all_matrix) - length(id_singletons), "样本 (",
    round((nrow(sc_all_matrix) - length(id_singletons))/nrow(sc_all_matrix)*100, 2), "%)\n")
cat("Call Rate 过滤:", nrow(sc_all_matrix) - length(id_call_rate), "样本 (",
    round((nrow(sc_all_matrix) - length(id_call_rate))/nrow(sc_all_matrix)*100, 2), "%)\n")
cat("最终保留:", length(common_ids), "样本 (",
    round(length(common_ids)/nrow(sc_all_matrix)*100, 2), "%)\n")
cat("总共过滤:", nrow(sc_all_matrix) - length(common_ids), "样本 (",
    round((nrow(sc_all_matrix) - length(common_ids))/nrow(sc_all_matrix)*100, 2), "%)\n")

# 写入结果
output_file <- "/opt/notebooks/sample_qc_final_keep.txt"
write.table(sample_qc_final_keep, output_file, sep="\t", row.names=FALSE, quote=FALSE)
cat("结果已保存到:", output_file, "\n")

# 可选：输出被过滤掉的样本及原因
failed_samples <- data.frame(
  IID = sc_all_matrix[[1]],
  Failed_Ti_Tv = !sc_all_matrix[[1]] %in% id_ti_tv,
  Failed_Het_Hom = !sc_all_matrix[[1]] %in% id_het_hom,
  Failed_Singletons = !sc_all_matrix[[1]] %in% id_singletons,
  Failed_Call_Rate = !sc_all_matrix[[1]] %in% id_call_rate
)

# 添加总体失败标志
failed_samples$Failed_Any = with(failed_samples, 
                                Failed_Ti_Tv | Failed_Het_Hom | 
                                  Failed_Singletons | Failed_Call_Rate)

# 只保留失败的样本
failed_samples_only <- subset(failed_samples, Failed_Any)

# 写入失败样本信息
if(nrow(failed_samples_only) > 0) {
  output_failed_file <- "/opt/notebooks/sample_qc_failed_reasons.txt"
  write.table(failed_samples_only, output_failed_file, sep="\t", row.names=FALSE, quote=FALSE)
  cat("失败样本信息已保存到:", output_failed_file, "\n")
}