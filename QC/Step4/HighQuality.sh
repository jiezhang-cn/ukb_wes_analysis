#!/bin/bash
# 设置工作目录
cd /opt/notebooks/

# 创建输出目录
mkdir -p ./QCStep4

# 下载并安装PLINK 2.0
mkdir -p ./plink2
wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip -O ./plink2/plink2.zip
cd ./plink2
unzip ./plink2.zip
chmod +x ./plink2
cd /opt/notebooks/

# 处理每条染色体并立即上传
for i in {1..22}; do
    echo "Processing chromosome ${i}..."
    
    # 运行PLINK2命令 - 添加线程参数
    ./plink2/plink2 \
    --bfile /mnt/project/Jzhang_data/exome_data/qc_step3/ukb_wes_step3_chr${i} \
    --maf 0.0001 \
    --geno 0.01 \
    --hwe 1e-6 \
    --indep-pairwise 200 100 0.1 \
    --make-bed \
    --threads 70 \
    --out ./QCStep4/HQ_chr${i}_v1
    
    # 运行第二次PLINK2命令
    ./plink2/plink2 \
    --bfile ./QCStep4/HQ_chr${i}_v1 \
    --indep-pairwise 200 100 0.05 \
    --make-bed \
    --threads 70 \
    --out ./QCStep4/HQ_chr${i}_v2
    
    # 检查PLINK2是否成功完成
    if [ $? -eq 0 ]; then
        echo "Chromosome ${i} processing completed. Uploading results..."
        
        # 上传结果
        dx upload -r ./QCStep4/HQ_chr${i}_v2.* --destination "project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/qc_step4_HQ"
        
        # 检查上传是否成功
        if [ $? -eq 0 ]; then
            echo "Upload for chromosome ${i} completed successfully."
            
            # 上传成功后删除中间文件以节省空间
            echo "Removing temporary files for chromosome ${i}..."
            rm -f ./QCStep4/HQ_chr${i}_v1.*
            
            echo "Cleanup for chromosome ${i} completed."
        else
            echo "Upload for chromosome ${i} failed."
        fi
    else
        echo "PLINK2 processing for chromosome ${i} failed."
    fi
done  # 添加了这个关键字

echo "All chromosomes processed and uploaded."