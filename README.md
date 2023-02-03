# Technical validation imputation
Author: Sonia Sarmiento

### This repository is a collection of resources and tools to use to assess the validity of imputation of low-coverage samples. 

## Introduction 
This repository tests the accuracy of imputation by comparing the results obtained from high-coverage and low-coverage sequencing of the same individuals. The incongruences between the two sets of data are analyzed to determine the reliability of the imputation process. This type of validation is important to ensure that the imputed data is accurate and can be used for meaningful analysis. 

### Section 1: Comparison with PCA

R code to compute a PCA from individuals sequenced at both high and low-coverage as a first view to look at imputation success. We plot high-coverage individuals (big squares) to compute the real distance between the samples, and reproject them at low-coverage (small circles) to see how far they fall from the high-coverage versions.

* Input: VCF containing samples in both high and low-coverage for one Chromosome (with different names, e.g. 'HC' and 'TC' suffix)
* Script: ```1_PCA.R```
* Output: Reprojected PCA ```1_Reprojected_PCA_Super-Scaffold_14.png```

