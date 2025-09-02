####  Annotation  ####
---

Cardiovascular-kidney-metabolic syndrome：
Cardiovascular-kidney-metabolic (CKM) syndrome is defined as a health disorder attributable to connections among obesity, diabetes, chronic kidney disease (CKD), and cardiovascular disease (CVD), including heart failure, atrial fibrillation, coronary heart disease, stroke, and peripheral artery disease. CKM syndrome includes those at risk for CVD and those with existing CVD.

A Presidential Advisory From the American Heart Association provides a CKM staging construct that reflects the pathophysiology, spectrum of risk, and opportunities for prevention and care optimization within CKM syndrome: 
stage 0, no CKM risk factors; 
stage 1, excess or dysfunctional adiposity; 
stage 2, metabolic risk factors (hypertriglyceridemia, hypertension, diabetes, metabolic syndrome) or moderate- to high-risk chronic kidney disease; 
stage 3, subclinical CVD in CKM syndrome or risk equivalents (high predicted CVD risk or very high-risk CKD); 
stage 4, clinical CVD in CKM syndrome.

---

ExWAS using POLMM-GENE:
1) Make Sparse GRM to adjust for sample relatedness (run_sparseGRM.sh)
2) Fit NULL model using Sparse GRM (fitNULLPOLMM.sh)
3) Single-variant tests using POLMM (single_variantPOLMM.sh)
4) Set-based tests using POLMM-GENE (gene_setPOLMM.sh)

---

Reference:
1) Ndumele CE, et al. Cardiovascular-Kidney-Metabolic Health: A Presidential Advisory From the American Heart Association. Circulation. 2023;148(20):1606-1635. 
2) Bi W, et al. Scalable mixed model methods for set-based association studies on large-scale categorical data analysis and its application to exome-sequencing data in UK Biobank. Am J Hum Genet. 2023;110(5):762-773.
3) Bi W, et al. Efficient mixed model approach for large-scale genome-wide association studies of ordinal categorical phenotypes. Am J Hum Genet. 2021;108(5): 825-839.

