#!/bin/bash

# 设置输入目录和参考基因组
INPUT_DIR="/mnt/project/Bulk/Exome sequences/Population level exome OQFE variants, pVCF format - final release"
REF_GENOME="GRCh38_full_analysis_set_plus_decoy_hla.fa"
BCFTOOLS="./bcftools/bcftools"
PLINK="./plink/plink"
MAX_JOBS=20  # 每批并行处理20个
BATCH_SIZE=20  # 每批文件数量

# 设置输出和上传目录
VCF_OUTPUT_DIR="filtered_vcfs"
PLINK_OUTPUT_DIR="plink_files"
VCF_UPLOAD_DESTINATION="project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/pvcf_filtered"
PLINK_UPLOAD_DESTINATION="project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/plink_qc_step1"

# 创建输出目录
mkdir -p $VCF_OUTPUT_DIR
mkdir -p $PLINK_OUTPUT_DIR

# 获取所有vcf.gz文件列表
find "$INPUT_DIR" -name "*.vcf.gz" > vcf_files.txt
TOTAL_FILES=$(wc -l < vcf_files.txt)
echo "找到 $TOTAL_FILES 个文件待处理"

# 定义处理单个文件的函数
process_file() {
    local input_file=$1
    local filename=$(basename "$input_file")
    local prefix="${filename%.vcf.gz}"
    local output_file="$VCF_OUTPUT_DIR/${prefix}.filtered.vcf.gz"
    
    echo "处理文件: $filename"
    
    $BCFTOOLS norm -f $REF_GENOME --threads 4 -m -any -Oz "$input_file" | $BCFTOOLS filter --threads 4 -i  '(GT="RR" &&  FORMAT/GQ>=20 && FORMAT/DP>=10 && FORMAT/DP<=200) || (GT="RA" && (FORMAT/AD[*:0]+FORMAT/AD[*:1])/(FORMAT/DP)>=0.9 && FORMAT/AD[*:1]/(FORMAT/DP)>=0.2 && FORMAT/PL[*:0]>=20 && FORMAT/DP>=10 && FORMAT/DP<=200) || (GT="AA" && (FORMAT/AD[*:0]+FORMAT/AD[*:1])/(FORMAT/DP)>=0.9 && FORMAT/AD[*:1]/(FORMAT/DP)>=0.9 && FORMAT/PL[*:0]>=20 && FORMAT/DP>=10 && FORMAT/DP<=200)' -Oz -o "$output_file" 
    
    # 为输出文件创建索引
    $BCFTOOLS index -t "$output_file"
    
    echo "完成处理: $filename"
}

# 定义PLINK处理函数
process_plink() {
    local filename=$1
    local prefix="${filename%.filtered.vcf.gz}"
    
    echo "使用PLINK处理文件: $filename"
    
    $PLINK --vcf "$VCF_OUTPUT_DIR/$filename" \
        --keep-allele-order --vcf-idspace-to _ --double-id \
        --allow-extra-chr 0 --make-bed --vcf-half-call m \
        --out "$PLINK_OUTPUT_DIR/$prefix"
    
    echo "PLINK处理完成: $prefix"
}

# 导出函数以便在子shell中使用
export -f process_file
export -f process_plink
export BCFTOOLS
export PLINK
export REF_GENOME
export VCF_OUTPUT_DIR
export PLINK_OUTPUT_DIR

# 按批次处理文件
BATCH=1
PROCESSED=0

while [ $PROCESSED -lt $TOTAL_FILES ]; do
    CURRENT_BATCH_SIZE=$((BATCH_SIZE < (TOTAL_FILES - PROCESSED) ? BATCH_SIZE : (TOTAL_FILES - PROCESSED)))
    
    echo "开始处理第 $BATCH 批 (共 $CURRENT_BATCH_SIZE 个文件)"
    
    # 提取当前批次的文件
    sed -n "$((PROCESSED + 1)),$((PROCESSED + CURRENT_BATCH_SIZE))p" vcf_files.txt > current_batch.txt
    
    # 并行处理当前批次
    cat current_batch.txt | xargs -I{} -P $MAX_JOBS bash -c 'process_file "$@"' _ {}
    
    # 获取当前批次处理后的所有文件
    find "$VCF_OUTPUT_DIR" -name "*.filtered.vcf.gz" > processed_files.txt
    
    # 使用PLINK处理所有过滤后的VCF文件
    echo "开始PLINK处理第 $BATCH 批文件..."
    cat processed_files.txt | xargs -n 1 basename | xargs -I{} -P $MAX_JOBS bash -c 'process_plink "$@"' _ {}
    
    echo "第 $BATCH 批处理完成，上传数据..."
    
    # 上传处理完的PLINK文件
    echo "上传PLINK文件到 $PLINK_UPLOAD_DESTINATION..."
    dx upload -r $PLINK_OUTPUT_DIR --destination "$PLINK_UPLOAD_DESTINATION"
    
    # 清理输出目录
    rm -r $VCF_OUTPUT_DIR/*
    rm -r $PLINK_OUTPUT_DIR/*
    
    # 更新计数器
    PROCESSED=$((PROCESSED + CURRENT_BATCH_SIZE))
    BATCH=$((BATCH + 1))
    
    echo "已处理 $PROCESSED / $TOTAL_FILES 个文件"
    echo "----------------------------------------------------"
done

# 清理临时文件
rm vcf_files.txt current_batch.txt processed_files.txt
rmdir $VCF_OUTPUT_DIR
rmdir $PLINK_OUTPUT_DIR

echo "所有文件处理完毕!"
