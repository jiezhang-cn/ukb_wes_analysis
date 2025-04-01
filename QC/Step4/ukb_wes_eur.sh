#!/bin/bash
# 设置工作目录
cd /opt/notebooks/

# 创建输出目录
mkdir -p ./QCStep4

# 下载并安装PLINK 2.0（如果不存在）
if [ ! -f "./plink2/plink2" ]; then
    echo "Installing PLINK 2.0..."
    mkdir -p ./plink2
    wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip -O ./plink2/plink2.zip
    cd ./plink2
    unzip ./plink2.zip
    chmod +x ./plink2
    cd /opt/notebooks/
    echo "PLINK 2.0 installed successfully."
else
    echo "PLINK 2.0 already installed."
fi

# 确保样本QC文件存在
if [ ! -f "ukb_wes_sample_qc_final_EUR_4sd.sample_list.txt" ]; then
    echo "Error: sample_qc_final_keep.txt not found!"
    exit 1
fi

# 处理每条染色体并立即上传，然后删除
for i in {1..22}; do
    echo "======================="
    echo "Processing chromosome ${i}..."
    echo "======================="
    
    # 运行PLINK2命令
    ./plink2/plink2 \
    --bfile /mnt/project/Jzhang_data/exome_data/qc_step3/ukb_wes_step3_chr${i} \
    --keep ./ukb_wes_sample_qc_final_EUR_4sd.sample_list.txt \
    --make-bed \
    --threads 30 \
    --out ./QCStep4/ukb_wes_eur_chr${i}
    
    # 检查PLINK2是否成功完成
    if [ $? -eq 0 ]; then
        echo "Chromosome ${i} processing completed successfully."
        
        # 检查生成的文件是否存在
        if [ -f "./QCStep4/ukb_wes_eur_chr${i}.bed" ] && [ -f "./QCStep4/ukb_wes_eur_chr${i}.bim" ] && [ -f "./QCStep4/ukb_wes_eur_chr${i}.fam" ]; then
            echo "Output files generated successfully."
            
            # 上传结果
            echo "Uploading chromosome ${i} files..."
            dx upload -r ./QCStep4/ukb_wes_eur_chr${i}.* --destination "project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/qc_step4/"
            
            # 检查上传是否成功
            if [ $? -eq 0 ]; then
                echo "Upload for chromosome ${i} completed successfully."
                
                # 删除文件以节省空间
                echo "Removing chromosome ${i} files to save space..."
                rm ./QCStep4/ukb_wes_eur_chr${i}.*
                echo "Files removed."
            else
                echo "Warning: Upload for chromosome ${i} failed. Files will be kept."
            fi
        else
            echo "Error: One or more output files for chromosome ${i} not found."
        fi
    else
        echo "Error: PLINK2 processing for chromosome ${i} failed."
    fi
    
    echo "Chromosome ${i} processing complete."
    echo ""
done

echo "All chromosomes processed and uploaded."

# 检查是否有任何文件剩余
remaining_files=$(ls ./QCStep4/ | wc -l)
if [ $remaining_files -eq 0 ]; then
    echo "All files have been processed, uploaded, and removed successfully."
else
    echo "Warning: $remaining_files files remain in the QCStep4 directory."
    ls -lh ./QCStep4/
fi

echo "Script completed at $(date)."