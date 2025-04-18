#!/bin/bash
# 设置工作目录
cd /opt/notebooks/

# 创建输出目录
mkdir -p ./input
mkdir -p ./output

for chr in {1..22}; do
  (echo -e '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'; 
   awk '{print $1"\t"$4"\t"$2"\t"$6"\t"$5"\t.\tPASS\t."}' /mnt/project/Jzhang_data/exome_data/qc_step4/ukb_wes_eur_chr${chr}.bim) > ./input/ukb_wes_eur_chr${chr}_snpeff.vcf
  echo "Processed chromosome ${chr}"
  dx upload -r ./input/ukb_wes_eur_chr${chr}_snpeff.vcf --destination "project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/Annotation/SnpEff/input/"
done


mkdir -p ./snpEff
cd ./snpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip


# 下载 GRCh38.105 数据库
apt-get update
apt-get install -y openjdk-21-jdk

cd snpEff/
echo -e "\n# GRCh38.105\nGRCh38.105.genome : Homo_sapiens\nGRCh38.105.reference : ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/" >> /opt/notebooks/snpEff/snpEff/snpEff.config

java -jar snpEff.jar download GRCh38.105



### Annotating mutation types to select putative loss-of-function variants
# snpEff GRCh38.105
for i in {1..22};
do
java -jar snpEff.jar  GRCh38.105 /opt/notebooks/input/ukb_wes_eur_chr${i}_snpeff.vcf \
> /opt/notebooks/output/ukb_wes_eur_chr${i}_SnpEff.vcf

dx upload -r /opt/notebooks/output/ukb_wes_eur_chr${i}_SnpEff.vcf --destination "project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/Annotation/SnpEff/output/"
done


### Annotating using 5 methods to define putative damaging missense variants
# snpEff dbnsfp
for i in $(seq 1 22); do
  echo "Processing chromosome $i..."
  
  # 步骤1：筛选missense变异
  echo "  Filtering missense variants..."
  java -jar /opt/notebooks/snpEff/snpEff/SnpSift.jar filter \
  "ANN[*].EFFECT has 'missense_variant'" \
  /mnt/project/Jzhang_data/exome_data/Annotation/SnpEff/output/ukb_wes_eur_chr${i}_SnpEff.vcf \
  > /opt/notebooks/input/ukb_wes_eur_chr${i}_missense.vcf
  
  # 步骤2：应用dbNSFP注释
  echo "  Annotating with dbNSFP..."
  java -jar /opt/notebooks/snpEff/snpEff/SnpSift.jar dbnsfp \
  -f genename,Ensembl_geneid,Uniprot_acc,SIFT_pred,Polyphen2_HDIV_pred,MutationTaster_pred,Polyphen2_HVAR_pred,CADD_phred,MetaSVM_pred,REVEL_score,AlphaMissense_pred \
  -db /opt/notebooks/dbNSFP5.1a_grch38.gz \
  -g hg38 \
  /opt/notebooks/input/ukb_wes_eur_chr${i}_missense.vcf \
  > /opt/notebooks/output/ukb_wes_eur_chr${i}_missense_SnpSift.vcf
  
  echo "Completed processing chromosome $i"
done

dx upload -r /opt/notebooks/output/ukb_wes_eur_chr*_missense_SnpSift.vcf --destination "project-GYXj5k0JzYzPbyBYx8V7BFg3:/Jzhang_data/exome_data/Annotation/SnpEff/output/"

echo "All chromosomes processed!"
