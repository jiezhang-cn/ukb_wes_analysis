#!/bin/bash
# 设置工作目录
cd /opt/notebooks/

# 创建输出目录
mkdir -p ./QCStep3

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
    
    # 运行PLINK2命令
    ./plink2/plink2 \
    --bfile /mnt/project/Jzhang_data/exome_data/qc_step2/ukb_wes_step2_chr${i} \
    --keep ./non_retracted_and_sex_qc_filtered.txt \
    --sample-counts 'cols'=hom,het,ts,tv,dipsingle,single,missing \
    --make-just-bim \
    --make-just-fam \
    --threads 30\
    --out ./QCStep3/ukb_wes_step3_outlier_chr${i}
    
    # 检查PLINK2是否成功完成
    if [ $? -eq 0 ]; then
        echo "Chromosome ${i} processing completed. Uploading results..."
        
        # 上传结果
        dx upload -r ./QCStep3/ukb_wes_step3_outlier_chr${i}.* --destination "project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/qc_step3_outlier"
        
        # 检查上传是否成功
        if [ $? -eq 0 ]; then
            echo "Upload for chromosome ${i} completed successfully."
        else
            echo "Upload for chromosome ${i} failed."
        fi
    else
        echo "PLINK2 processing for chromosome ${i} failed."
    fi
done

echo "All chromosomes processed and uploaded."