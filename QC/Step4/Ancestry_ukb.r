# 加载必要的库
library(data.table)
library(dplyr)

# 读取合并后的数据
het_data <- fread("C:/Users/张杰/Desktop/whole_exome_sequencing_analyses/codes/QC/Step4/ukb_wes_sample_qc_final_all.het_check.csv")

# 设置概率阈值
prob_threshold <- 0.95

# 1. 处理欧洲血统(EUR)样本，基于主成分分析
eur_samples <- het_data %>% 
  filter(`ancestry-prediction` == "EUR")

# 计算前四个主成分的均值和标准差
pc_stats <- data.frame(
  PC = c("PC1", "PC2", "PC3", "PC4"),
  mean = c(mean(eur_samples$PC1), mean(eur_samples$PC2), mean(eur_samples$PC3), mean(eur_samples$PC4)),
  sd = c(sd(eur_samples$PC1), sd(eur_samples$PC2), sd(eur_samples$PC3), sd(eur_samples$PC4))
)

# 根据±4标准差筛选欧洲样本
restricted_eur <- eur_samples %>%
  filter(
    PC1 >= (pc_stats$mean[1] - 4 * pc_stats$sd[1]) & 
      PC1 <= (pc_stats$mean[1] + 4 * pc_stats$sd[1]) &
      PC2 >= (pc_stats$mean[2] - 4 * pc_stats$sd[2]) & 
      PC2 <= (pc_stats$mean[2] + 4 * pc_stats$sd[2]) &
      PC3 >= (pc_stats$mean[3] - 4 * pc_stats$sd[3]) & 
      PC3 <= (pc_stats$mean[3] + 4 * pc_stats$sd[3]) &
      PC4 >= (pc_stats$mean[4] - 4 * pc_stats$sd[4]) & 
      PC4 <= (pc_stats$mean[4] + 4 * pc_stats$sd[4])
  )

# 打印EUR筛选结果
cat("EUR 原始样本数:", nrow(eur_samples), "\n")
cat("EUR ±4标准差筛选后的样本数:", nrow(restricted_eur), "\n")
cat("EUR 移除的样本数:", nrow(eur_samples) - nrow(restricted_eur), "\n\n")

# 保存EUR结果
write.csv(restricted_eur, "C:/Users/张杰/Desktop/whole_exome_sequencing_analyses/codes/QC/Step4/ukb_wes_sample_qc_final_EUR_4sd.csv", row.names = FALSE)

# 提取EUR样本ID列表 - 修改后的代码
eur_ids <- gsub("_.*$", "", restricted_eur$sample_id)
eur_sample_ids <- data.frame(FID = eur_ids, IID = eur_ids)

# 保存EUR样本ID列表
write.table(eur_sample_ids, "C:/Users/张杰/Desktop/whole_exome_sequencing_analyses/codes/QC/Step4/ukb_wes_sample_qc_final_EUR_4sd.sample_list.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)






# 2. 处理其他血统(SAS, AFR, EAS, AMR)，仅基于概率阈值
ancestries <- c("SAS", "AFR", "EAS", "AMR")
prob_threshold_alt <- 0.6

for (ancestry in ancestries) {
  # 筛选特定血统且概率>阈值的样本
  ancestry_samples <- het_data %>% 
    filter(`ancestry-prediction` == ancestry & `ancestry-prob` > prob_threshold_alt)
  
  # 打印结果
  cat(ancestry, "血统，概率>", prob_threshold_alt, "的样本数:", nrow(ancestry_samples), "\n")
  
  # 保存结果
  out_csv <- paste0("C:/Users/张杰/Desktop/whole_exome_sequencing_analyses/codes/QC/Step4/ukb_wes_sample_qc_final_", ancestry, "_prob", prob_threshold_alt, ".csv")
  write.csv(ancestry_samples, out_csv, row.names = FALSE)
  
  # 提取样本ID列表 - 修改后的代码
  ancestry_ids <- gsub("_.*$", "", ancestry_samples$sample_id)
  ancestry_sample_ids <- data.frame(FID = ancestry_ids, IID = ancestry_ids)
  
  # 保存样本ID列表
  out_txt <- paste0("C:/Users/张杰/Desktop/whole_exome_sequencing_analyses/codes/QC/Step4/ukb_wes_sample_qc_final_", ancestry, "_prob", prob_threshold_alt, ".sample_list.txt")
  write.table(ancestry_sample_ids, out_txt, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# 3. 保存所有筛选后样本的总数据集
all_filtered_samples <- rbind(
  restricted_eur,
  het_data %>% filter(`ancestry-prediction` == "SAS" & `ancestry-prob` > prob_threshold_alt),
  het_data %>% filter(`ancestry-prediction` == "AFR" & `ancestry-prob` > prob_threshold_alt),
  het_data %>% filter(`ancestry-prediction` == "EAS" & `ancestry-prob` > prob_threshold_alt),
  het_data %>% filter(`ancestry-prediction` == "AMR" & `ancestry-prob` > prob_threshold_alt)
)

# 保存所有筛选后的样本
write.csv(all_filtered_samples, "C:/Users/张杰/Desktop/whole_exome_sequencing_analyses/codes/QC/Step4/ukb_wes_sample_qc_final_all_filtered.csv", row.names = FALSE)

# 打印总体结果
cat("\n总结：\n")
cat("原始数据样本总数:", nrow(het_data), "\n")
cat("筛选后样本总数:", nrow(all_filtered_samples), "\n")
cat("总移除样本数:", nrow(het_data) - nrow(all_filtered_samples), "\n")
