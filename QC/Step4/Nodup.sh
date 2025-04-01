## exclude duplicated (result: king.con)

# 设置工作目录
cd /opt/notebooks/

# 下载并安装PLINK 2.0
mkdir -p ./plink2
wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip -O ./plink2/plink2.zip
cd ./plink2
unzip ./plink2.zip
chmod +x ./plink2
cd /opt/notebooks/


# 创建并进入下载目录
mkdir -p ./king
cd ./king

# 下载KING (Linux版本)
wget https://www.kingrelatedness.com/Linux-king.tar.gz

# 解压文件
tar -xzvf Linux-king.tar.gz


cd /opt/notebooks/


# 创建合并列表文件
ls -1 /mnt/project/Jzhang_data/exome_data/qc_step4_HQ/HQ_chr*_v2.bed | sed 's/.bed$//' > merge_list.txt

# 检查列表文件内容
cat merge_list.txt

# 合并step4所有的染色体highquality_variants
    ./plink2/plink2 \
      --pmerge-list merge_list.txt "bfile" \
      --make-bed \
      --out /opt/notebooks/ukb_wes_chr_all_king_sample_qc
      
      
      
#去除重复个体
./king/king -b /opt/notebooks/ukb_wes_chr_all_king_sample_qc.bed --duplicate --cpus 12