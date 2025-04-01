#本地服务器运行祖先成分预测分析
conda activate UKB_R_4_2

cd /home/user2/jzhang_data/ukb_wes_qc_data/relationship_ancestry/peddy_sample_qc_final

# 批次大小设置
BATCH_SIZE=5000
TOTAL_SAMPLES=434503
BATCH_COUNT=$((TOTAL_SAMPLES / BATCH_SIZE + 1))

# 循环处理每个批次
for ((i=1; i<=$BATCH_COUNT; i++)); do
    START_LINE=$((($i-1)*$BATCH_SIZE + 1))
    
    echo "正在处理批次 $i (从行 $START_LINE 开始)"
    
    # 创建当前批次的FAM文件
    awk -v start=$START_LINE -v size=$BATCH_SIZE 'NR>=start && NR<start+size' \
        ../ukb_wes_chr_all_king_sample_qc_final_unrelated_modified.fam > peddy_batch$i.fam
    
    # 创建当前批次的样本ID列表
    awk -v start=$START_LINE -v size=$BATCH_SIZE 'NR>=start && NR<start+size {print $1, $2}' \
        ../ukb_wes_chr_all_king_sample_qc_final_unrelated.fam > batch$i.txt
    
    # 提取当前批次的样本并转换为VCF格式
    ../plink2/plink2 --bfile ../ukb_wes_chr_all_king_sample_qc_final_unrelated \
        --keep batch$i.txt --recode vcf --out batch$i
    
    # 压缩和索引VCF文件
    bgzip batch$i.vcf
    tabix -p vcf batch$i.vcf.gz
    
    # 运行peddy预测祖先成分
    python -m peddy -p 72 --prefix ukb_wes_sample_qc_final_batch$i --sites hg38 \
        ./batch$i.vcf.gz ./peddy_batch$i.fam
    
    echo "批次 $i 处理完成"
done


#!/bin/bash

# 设置工作目录
cd /home/user2/jzhang_data/ukb_wes_qc_data/relationship_ancestry/peddy_sample_qc_final

# 输出文件名
OUTPUT_FILE="ukb_wes_sample_qc_final_all.het_check.csv"

# 获取第一个文件，用于提取标题行
FIRST_FILE=$(ls ukb_wes_sample_qc_final_batch*.het_check.csv | sort -V | head -n 1)

# 提取标题行并写入输出文件
head -n 1 "$FIRST_FILE" > "$OUTPUT_FILE"

# 逐个处理每个批次文件，跳过标题行
for file in $(ls ukb_wes_sample_qc_final_batch*.het_check.csv | sort -V); do
    echo "处理文件: $file"
    # 跳过第一行（标题行），将数据追加到输出文件
    tail -n +2 "$file" >> "$OUTPUT_FILE"
done

# 计算合并后的行数
TOTAL_LINES=$(wc -l < "$OUTPUT_FILE")
DATA_LINES=$((TOTAL_LINES - 1))  # 减去标题行

echo "合并完成! 总共处理了87个批次文件"
echo "合并文件 $OUTPUT_FILE 包含 $DATA_LINES 个样本数据记录（共 $TOTAL_LINES 行，包括标题行）"