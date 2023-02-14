# Technical validation of imputation  🛠️
Author: Sonia Sarmiento

### Description: This repository is a collection of resources and tools to use to assess the validity of imputation of low-coverage samples. 
#### Coding: The imputation is done in a HPC cluster in bash, but the validation and all the codes in this repository are done in R. 

## Introduction 
This repository tests the accuracy of imputation by comparing the results obtained from high-coverage and low-coverage sequencing of the same individuals. The incongruences between the two sets of data are analyzed to determine the reliability of the imputation process. This type of validation is important to ensure that the imputed data is accurate and can be used for meaningful analysis. 

### Section 1: Comparison with PCA
PCA comparison can evaluate the accuracy of imputation in low-coverage samples by comparing the results of high-coverage and low-coverage sequencing of the same individuals. The imputed data is visualized in a lower-dimensional space along with high-coverage data to see how well it aligns. Close alignment suggests accurate imputation, while significant differences suggest potential issues with the process. 

* Input: VCF containing samples in both high and low-coverage for one Chromosome (with different names, e.g. 'HC' and 'TC' suffix)
* Script: ```1_PCA.R```
* Output: Reprojected PCA

In the results, high-coverage individuals are depicted as big squares and the low-coverage in small circles.


### Section 2: Comparison of genotypes along the chromosome
This code is designed to compare genotypes between two VCF files. Specifically, the goal is to evaluate the accuracy imputation by comparing imputed low-coverage sequences to high-coverage sequences of the same individual. 

* Input: High-Coverage and Low-Coverage GDS* for replicate individuals
* Script: ```2_Genotype_mismatches.md```
* Output: Distribution of mismatches and missing values and Genotype mismatches along the genome. 

**IMPORTANT**: input must include only SNPs existing in both VCFs. For this purpose, we can use the following code in bash:
```bash
# Keep intersect snps of unrelated
bcftools isec -n=2 -w1 -O z -o highcoverage_intersect.vcf.gz highcoverage.vcf.gz lowcoverage.vcf.gz

# Keep intersect snps of related
bcftools isec -n=2 -w1 -O z -o lowcoverage_intersect.vcf.gz lowcoverage.vcf.gz highcoverage.vcf.gz

```

* To obtain GDS, you can load gcc and r, open R, and run the following command:
```r
library(SNPRelate) 
vcf <- 'name.vcf.gz'
snpgdsVCF2GDS(vcf,'name.gds')
```
