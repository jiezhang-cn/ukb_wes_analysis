# Create king matrix of 2nd degree (ukb_wes_chr_all_king_sample_qc.kin0)

# 设置工作目录
cd /opt/notebooks/

# 下载并安装PLINK 2.0
mkdir -p ./plink2
wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip -O ./plink2/plink2.zip
cd ./plink2
unzip ./plink2.zip
chmod +x ./plink2
cd /opt/notebooks/

/plink2/plink2 \
--bfile /opt/notebooks/ukb_wes_chr_all_king_sample_qc \
--make-king \
--king-table-filter 0.0884 \
--make-king-table \
--threads 60 \
--out /opt/notebooks/ukb_wes_chr_all_king_sample_qc