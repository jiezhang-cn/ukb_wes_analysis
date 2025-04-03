#!/bin/bash
# UKB全外显子组数据ANNOVAR注释脚本

# 进入工作目录
cd /opt/notebooks/

# 下载并设置ANNOVAR（如果尚未安装）
if [ ! -d "annovar" ]; then
  echo "下载ANNOVAR..."
  wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
  tar -zxvf annovar.latest.tar.gz
  rm annovar.latest.tar.gz
  echo "ANNOVAR下载并解压完成"
fi

# 创建humandb链接（如果需要）
if [ ! -d "humandb" ] && [ -d "annovar/humandb" ]; then
  ln -s annovar/humandb/ humandb
elif [ ! -d "humandb" ]; then
  mkdir -p annovar/humandb
  ln -s annovar/humandb/ humandb
fi

# 下载必要的数据库（如果缺失）
cd annovar
for DB in refGene ensGene gnomad211_exome gnomad30_genome dbnsfp42a; do
  if [ ! -e "humandb/${DB}_hg38.txt" ] && [ ! -e "humandb/${DB}_hg38.txt.gz" ]; then
    ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar $DB humandb/
  fi
done
cd /opt/notebooks/

# 第一部分：处理bim文件并生成ANNOVAR输入文件
for i in {1..22}; do
  INPUT_BIM="/mnt/project/Jzhang_data/exome_data/qc_step4/ukb_wes_eur_chr${i}.bim"
  if [ -f "$INPUT_BIM" ]; then
    awk -v OFS="\t" '{print $1,$4,$4,$6,$5,$2}' "$INPUT_BIM" > ./input/ukb_wes_eur_chr${i}.avinput
    dx upload -r ./input/ukb_wes_eur_chr${i}.avinput --destination "project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/Annotation/Annovar/input/"
  fi
done

# 第二部分：运行ANNOVAR注释
for i in {1..22}; do
  if [ -f "./input/ukb_wes_eur_chr${i}.avinput" ]; then
    perl annovar/table_annovar.pl ./input/ukb_wes_eur_chr${i}.avinput humandb/ \
    -buildver hg38 \
    -out ./output/ukb_wes_eur_chr${i} \
    -remove \
    -protocol refGene,gnomad211_exome,gnomad30_genome,dbnsfp42a,ensGene \
    -operation g,f,f,f,g \
    -nastring . \
    -csvout
    
    dx upload -r ./output/ukb_wes_eur_chr${i}.* --destination "project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/Annotation/Annovar/output/"
  fi
done

echo "ANNOVAR注释流程完成"