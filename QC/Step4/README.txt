#### Step4 : Sample QC part-2 ####

In this part, further QC were performed:
· kinship calculation, 3nd degree (>0.0884) and above relatives were excluded
· remove duplicated samples between instances / centers
· calculate ancestry principal components
· predicted genetic ancestries


1) We firstly derived a high-quality SNPs subset from the ukb exome sequence data (after Steps 1-3 QC) by using `HighQuality.sh` which filtered:
· a minor allele frequency (MAF) > 0.1%
· call rate > 99%
· HWE P-value > 1x10^-6. 
· two rounds of pruning (Plink parameters: -indep-pairwise 200 1000.1 and -indep-pairwise 200 100 0.05). 

2) We then identified duplicated individuals were using KING software in `Nodup.sh`.

3) We sought to remove related individuals to mitigate a possible increase of variance estimates from our analyses. 
Using exome sequence-derived high-quality SNPs subset (as the aforementioned in `HighQuality.sh`) as input to the kinship algorithm included in KING software, we generated pairwise kinship coefficients for all remaining samples in `Relative.sh`.
  We used the ukb_gen_samples_to_remove() function from the R package ukbtools to choose a subset of individuals within which no pair had a kinship coefficient exceeding 0.0884, 
equivalent of up to third-degree relatives. For each related pair, this function removes whichever member has the highest number or relatives above the provided threshold, resulting in a maximal set. 

4) Next, we removed duplicated and third-degree related samples from the high-quality SNPs subset. 
Principal components of ancestry (PCA) was then calculated using high-quality SNPs subset by running --pca approx function in Plink (here not provided the script), 
and the first 10 PCA were used as covariates in association tests for common and rare variants.

5) Finally, we predicted genetic ancestries from the high-quality SNPs subset using peddy v0.4.2 with the ancestry labelled 1,000 Genomes Project as reference, as in `Predicted_ancestry.sh`.  
We further restricted the European ancestry cohort to those (410,831) had a Pr(European) ancestry prediction of more than 0.9 and within ±4 s.d. across the top four principal component means;
Through this step (`Ancestry_ukb.r`), we also identified and used 8,344 (Pr(African) > 0.6), 1,172 (Pr(Hispanic or Latin American) > 0.6), 
703 (Pr(East Asian) > 0.6) and 8,911 (Pr(South Asian) > 0.6) ukb participants for ancestry-independent analyses. 

Reference:
1) Jurgens, S. J. et al. Analysis of rare genetic variation underlying cardiometabolic diseases and traits among 200,000 individuals in the UK Biobank. Nat. Genet. 54, 240–250 (2022).
2) Yang, L., Ou, YN., Wu, BS. et al. Large-scale whole-exome sequencing analyses identified protein-coding variants associated with immune-mediated diseases in 350,770 adults. Nat Commun 15, 5924 (2024). 
3) Wang, Q., Dhindsa, R.S., Carss, K. et al. Rare variant contribution to human disease in 281,104 UK Biobank exomes. Nature 597, 527–532 (2021).
4) Backman, J.D., Li, A.H., Marcketta, A. et al. Exome sequencing and analysis of 454,787 UK Biobank participants. Nature 599, 628–634 (2021). 
5) Pedersen and Quinlan, Who’s Who? Detecting and Resolving Sample Anomalies in Human DNA Sequencing Studies with Peddy, The American Journal of Human Genetics (2017).
