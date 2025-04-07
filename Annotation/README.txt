####  Annotation  ####
---

ANNOVAR：
Common variants (MAF>1%) were mapped to genes by using ANNOVAR. 

---

SnpEff:
In the SAIGE analysis, we used SnpEff annotation to classify mutation types: 
1) Predicted loss-of-function (plof) variants were defined as  Frameshift insertion/deletion, splice acceptor, splice donor, stop gain, start loss, and stop loss mutations.
2) Predicted deleterious missense (pmis) variants were defined as those consistently predicted to be deleterious by five in silico tools: SIFT, LRT, PolyPhen2 HDIV, PolyPhen2 HVAR, and MutationTaster. 

We then organized the results using plof_pmis_filter.py. GroupFileGenerate.py was used to transfer the annotation results into the format of SAIGE analyses.


---

Reference:
1) Van Hout, C. V. et al. Exome sequencing and characterization of 49,960 individuals in the UK Biobank. Nature 586, 749–756 (2020).
2) Backman, J.D., Li, A.H., Marcketta, A. et al. Exome sequencing and analysis of 454,787 UK Biobank participants. Nature 599, 628–634 (2021). 
3) Cingolani, P. et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly 6, 80–92 (2012). 

