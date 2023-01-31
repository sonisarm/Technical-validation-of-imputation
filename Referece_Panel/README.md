Using a reference panel when performing low-coverage genotype calling and imputation increases accuracy.
A reference panel consists of high-coverage sequences that can be phased using read-base, pedigree and population information. 
To assess the best type of phasing, we compared phasing outputs with Mendelian Inheritance obtained from trios (1 offspring - 2 parents). 
For this, we used 6 confirmed offspring by our pedigree from our dataset, which were phased with other unrelated samples. 

Input:
1) Unphased VCFs (parents + offspring)
2) Phased VCFs (offspring only)
3) PED file - columns: offspring | dad id | mom id
