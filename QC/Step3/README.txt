#### Step3 : Sample QC part-1 ####

We then conducted the Sample QC to remove individuals with:
· Ti/Tv ratio, heterozygote/homozygote ratio, and number of singletons beyond 8 standard deviations
· Sample call rate < 90%. 
· discordance genetic sex and self-reported sex
· participants who had rescinded their consent

---

First, the Outlier.sh script uses plink to calculate the Ti, Tv, heterozygote, homozygote, and singleton counts, along with the number of sequenced sites for each individual.
Then, using the result from Outlier.sh, Calculation.R computes the ratios, and call rate by dividing each individual's variant count by the total variant count, retaining the IDs of individuals who pass QC. 
Finally, the Sample_QC.sh script keeps only those individuals who passed the sample QC.

---

Please Note: this Step3 QC also removes individuals who had rescinded their consent (374 ukb participants withdrawal until to 2024/12/17), 
and those mismatched genetic sex (UKB data-filed 22001) and self-reported sex (UKB data-filed 31) according to UKB's official calculations. These processes are not provided in the above scripts!

Reference:
1) Yang, L., Ou, YN., Wu, BS. et al. Large-scale whole-exome sequencing analyses identified protein-coding variants associated with immune-mediated diseases in 350,770 adults. Nat Commun 15, 5924 (2024).
