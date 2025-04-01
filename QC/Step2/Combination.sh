#!/bin/bash

# 设置工作目录
cd /opt/notebooks/

# 创建输出目录
mkdir -p /opt/notebooks/merged_chr

# 下载并安装PLINK 2.0
mkdir -p ./plink2
wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip -O ./plink2/plink2.zip
cd ./plink2
unzip ./plink2.zip
chmod +x ./plink2
cd /opt/notebooks/

# 获取唯一的染色体编号
chromosomes=$(ls /mnt/project/Jzhang_data/exome_data/qc_step1/ukb23157_c*_b*_v1.bed | grep -o 'c[0-9XY]\+_' | sed 's/c//;s/_//' | sort -n | uniq)

echo "发现的染色体: $chromosomes"

# 处理每条染色体
for chr in $chromosomes; do
  echo "处理染色体 $chr"
  
  # 创建临时合并列表
  merge_list=$(mktemp merge_list_chr${chr}.XXXXXX)
  
  # 查找此染色体的所有批次文件并添加到合并列表中（不带扩展名）
  for bed_file in $(ls /mnt/project/Jzhang_data/exome_data/qc_step1/ukb23157_c${chr}_b*_v1.bed | sort -V); do
    base_name=$(basename $bed_file .bed)
    echo "/mnt/project/Jzhang_data/exome_data/qc_step1/$base_name" >> $merge_list
  done
  
  # 计算要合并的文件数
  file_count=$(wc -l < $merge_list)
  echo "发现染色体 $chr 的 $file_count 个文件"
  
  if [ $file_count -eq 0 ]; then
    echo "染色体 $chr 没有找到文件，跳过"
    rm $merge_list
    continue
  fi
  
  if [ $file_count -eq 1 ]; then
    # 只有一个文件 - 直接复制
    first_file=$(head -n 1 $merge_list)
    echo "染色体 $chr 只有一个文件，直接复制"
    ./plink2/plink2 \
      --bfile $first_file \
      --make-bed \
      --out /opt/notebooks/merged_chr/ukb_wes_chr${chr}_merged
  else
    # 多个文件 - 使用PLINK 2.0一次合并所有批次
    echo "使用PLINK 2.0一次性合并染色体 $chr 的所有 $file_count 个批次文件"
    
    # 创建PLINK 2.0格式的合并列表
    plink2_merge_list=$(mktemp plink2_merge_list_chr${chr}.XXXXXX)
    cat $merge_list > $plink2_merge_list
    
    # 使用pmerge-list命令一次性合并所有文件
    # 添加 --set-all-var-ids 参数解决多等位基因变体问题
    # 添加 --new-id-max-allele-len 200 允许处理长等位基因代码
    ./plink2/plink2 \
      --pmerge-list $plink2_merge_list "bfile" \
      --set-all-var-ids "@:#:\$r:\$a" \
      --new-id-max-allele-len 500 \
      --make-bed \
      --out /opt/notebooks/merged_chr/ukb_wes_chr${chr}_merged
    
    # 检查合并是否成功
    if [ -f /opt/notebooks/merged_chr/ukb_wes_chr${chr}_merged.bed ]; then
      echo "染色体 $chr 的所有批次合并成功"
    else
      echo "染色体 $chr 合并失败，查看日志了解详情"
      
      # 如果仍然失败，我们可以考虑使用替代策略
      # 例如，尝试使用更简单的ID模式
      echo "尝试使用更简单的ID模式重新合并"
      ./plink2/plink2 \
        --pmerge-list $plink2_merge_list "bfile" \
        --set-all-var-ids "@:#" \
        --make-bed \
        --out /opt/notebooks/merged_chr/ukb_wes_chr${chr}_merged
        
      # 如果还是失败，可以考虑使用替代方法
      if [ ! -f /opt/notebooks/merged_chr/ukb_wes_chr${chr}_merged.bed ]; then
        echo "替代ID模式也失败了，染色体 $chr 可能需要手动处理"
      fi
    fi
    
    # 清理临时文件
    rm $plink2_merge_list
  fi
  
  # 上传此染色体的合并文件
  echo "上传染色体 $chr 的合并文件"
  if [ -f /opt/notebooks/merged_chr/ukb_wes_chr${chr}_merged.bed ]; then
    dx upload -r /opt/notebooks/merged_chr/ukb_wes_chr${chr}_merged* --destination "project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/merged_chr"
    
    # 上传完成后删除本地合并文件以节省空间
    echo "删除本地合并文件以节省空间"
    rm -f /opt/notebooks/merged_chr/ukb_wes_chr${chr}_merged*
  else
    echo "警告：染色体 $chr 的合并文件不存在，跳过上传"
  fi
  
  # 清理临时文件
  rm $merge_list
  echo "完成处理染色体 $chr"
done

echo "所有染色体处理完成"