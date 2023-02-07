# Technical validation imputation  🛠️
Author: Sonia Sarmiento

### This repository is a collection of resources and tools to use to assess the validity of imputation of low-coverage samples. 

## Introduction 
This repository tests the accuracy of imputation by comparing the results obtained from high-coverage and low-coverage sequencing of the same individuals. The incongruences between the two sets of data are analyzed to determine the reliability of the imputation process. This type of validation is important to ensure that the imputed data is accurate and can be used for meaningful analysis. 

### Section 1: Comparison with PCA
PCA comparison can evaluate the accuracy of imputation in low-coverage samples by comparing the results of high-coverage and low-coverage sequencing of the same individuals. The imputed data is visualized in a lower-dimensional space along with high-coverage data to see how well it aligns. Close alignment suggests accurate imputation, while significant differences suggest potential issues with the process. 

* Input: VCF containing samples in both high and low-coverage for one Chromosome (with different names, e.g. 'HC' and 'TC' suffix)
* Script: ```1_PCA.R```
* Output: Reprojected PCA ```1_Reprojected_PCA_Super-Scaffold_14.png```

In the results, high-coverage individuals are depicted as big squares and the low-coverage in small circles (in our data, some samples have 3 low-coverage sequences (triplicates) and other 1 (duplicates)). 


### Section 2: Comparison of genotypes along the chromosome

* Input: High-Coverage and Low-Coverage GDS for replicate individuals
* Script: ```2_Genotype_mismatches.md```
* Output: Distribution of mismatches and missing values (```2_Mismatches and missing genotypes distribution.jpg```) and Genotype mismatches along the Super-Scaffold_14 (biggest scaffold in barn owls) for individuals with high and low proportion of mismatches separately (```2_Mismatches between high- and imputed low-coverage - colors - separated indvs.jpg```).

