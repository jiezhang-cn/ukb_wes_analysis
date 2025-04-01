#### Step2: Variant QC ####

We then conducted the Variant QC on the expored plink files of exome data in the Step 1 QC, as following:

1) We merged the plink files by chromosome using the `Combination.sh` script.
2) Variant quality control, we removed the variants using the `Variant_QC.sh` script:
  • a call rate of less than 90% 
  • those deviating from Hardy-Weinberg equilibrium (P<1×10−15). 
  • those deviating present in regions of low complexity

 The ref file for the regions of low complexity was download from https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs38.bed.gz

Reference:
1) Li H. Toward better understanding of artifacts in variant calling from high-coverage samples. Bioinformatics. 2014;30(20):2843-2851.
