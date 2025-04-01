#!/bin/bash
# 设置工作目录
cd /opt/notebooks/

# 创建输出目录
mkdir -p ./QCStep2

# 下载并安装PLINK 2.0
mkdir -p ./plink2
wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip -O ./plink2/plink2.zip
cd ./plink2
unzip ./plink2.zip
chmod +x ./plink2
cd /opt/notebooks/


# 获取合并文件列表（只获取基本文件名，去掉扩展名）
ls /mnt/project/Jzhang_data/exome_data/merged_chr/*.bed | sed 's/\.bed$//' | sort -u > merged_full_paths.txt
sed 's|.*/||' merged_full_paths.txt > merged.txt

# 下载低复杂度区域文件
wget https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs38.bed.gz
gunzip -f LCR-hs38.bed.gz  # 解压文件
awk '{print $0,NR}' LCR-hs38.bed > LCR-hs38.txt

# 一步完成所有质控处理
while read basename; do
  # 从文件名中提取染色体号
  chr=$(echo $basename | grep -o 'chr[0-9]\+' | sed 's/chr//')
  
  echo "处理染色体 $chr 的质量控制..."
  
  # 合并步骤: 同时过滤基因型缺失率>0.1、HWE p-value<1E-15和低复杂度区域
  ./plink2/plink2 --bfile /mnt/project/Jzhang_data/exome_data/merged_chr/${basename} \
         --geno 0.1 \
         --hwe 1e-15 \
         --exclude range LCR-hs38.txt \
         --make-bed \
         --out ./QCStep2/ukb_wes_step2_chr${chr}
  
  echo "染色体 $chr 质量控制完成"
  
  # 上传此染色体的质控后文件
  echo "上传染色体 $chr 的质控后文件"
  if [ -f ./QCStep2/ukb_wes_step2_chr${chr}.bed ]; then
    dx upload -r ./QCStep2/ukb_wes_step2_chr${chr}* --destination "project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/qc_step2"
    
    # 上传完成后删除本地文件以节省空间
    echo "删除本地质控后文件以节省空间"
    rm -f ./QCStep2/ukb_wes_step2_chr${chr}*
  else
    echo "警告：染色体 $chr 的质控后文件不存在，跳过上传"
  fi
  
  echo "完成处理染色体 $chr"
done < merged.txt

# 清理临时文件
rm -f merged.txt merged_full_paths.txt LCR-hs38.bed LCR-hs38.txt

echo "所有染色体的质量控制处理完成"