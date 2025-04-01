#### Step1 : Genotype QC ####

We firstly conducted the Genotype QC on UKB Population level exome OQFE variants, pVCF format - Final exome release (UKB data-filed 23157), as following:

1) Split multi-allelic sites to represent separate bi-allelic sites
2) Genotype quality control:
  • For homozygous reference calls: Genotype Quality < 20; Genotype Depth < 10; Genotype Depth > 200  
  • For heterozygous calls: (A1 Depth + A2 Depth)/Total Depth < 0.9; A2 Depth/Total Depth < 0.2; Genotype likelihood[ref/ref] < 20; Genotype Depth < 10; Genotype Depth > 200  
  • For homozygous alternative calls: (A1 Depth + A2 Depth)/Total Depth < 0.9; A2 Depth/Total Depth < 0.9; Genotype likelihood[ref/ref] < 20; Genotype Depth < 10; Genotype Depth > 20; Genotype Depth > 200

`bcftools_process.sh` is set up to run in the UKB RAP platform, due to the UKB exome data are not allowed to be downloaded locally. This bash script will use bcftools to complete the initial genotype QC on exome data of pVCF format, 
and export the after-QC exome data to plink format (.bed, .bim, .fam).

Please note: the entire UKB exome data (exported plink files) are huge ~2.4 TB, `bcftools_process.sh` will run >5 days using the instance of mem3_ssd1_v2_x48. You may balance the time and costs, and choose a disfferent 
instance in UKB RAP.

Reference:
1) Jurgens, S. J. et al. Analysis of rare genetic variation underlying cardiometabolic diseases and traits among 200,000 individuals in the UK Biobank. Nat. Genet. 54, 240–250 (2022).
2) Yang, L., Ou, YN., Wu, BS. et al. Large-scale whole-exome sequencing analyses identified protein-coding variants associated with immune-mediated diseases in 350,770 adults. Nat Commun 15, 5924 (2024). 
