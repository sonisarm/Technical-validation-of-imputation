##########################################################
### Author: Sarmiento Cabello, Sonia                   ###
### Version: 1.0.                                      ###
### Objective: PCA to compare same individuals         ###
###            sequenced at  high & low-coverage.      ###
##########################################################

# Load libraries
library("SNPRelate")

# Load data (VCF with all individuals that you want to compare in both high and low coverage (with different names, e.g. IND1_HC, IND1_LC))
vcf.fn<-"HighCov_and_LowCov_Super-Scaffold_14.vcf.gz" 

#Produce GDS
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
snpgdsSummary("ccm.gds")
genofile <- openfn.gds("ccm.gds")

# get sample names
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
sample.id

#create pop.txt by reading the samples vertically and removing the TC or HC
cat(paste( sample.id, collapse='\n' ) )

#open pop.txt and check that they match
pop_code <- scan("pop.txt", what=character())
head(cbind(sample.id, pop_code))

# PCA
ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE)
# variance proportion (%)
pc.percent <- ccm_pca$varprop*100
head(round(pc.percent, 2))

#### First Plot (ALL) - no reprojection ####
tab <- data.frame(sample.id = ccm_pca$sample.id, pop = factor(pop_code)[match(ccm_pca$sample.id, sample.id)], 
        EV1 = ccm_pca$eigenvect[,1],    # the first eigenvector
        EV2 = ccm_pca$eigenvect[,2],    # the second eigenvector
        stringsAsFactors = FALSE)

#Generate colors
library(RColorBrewer)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
pie(rep(1,32), col=sample(color, 32))
col.list=sample(color, 32)
#Plot
# 1) change size of HC
temp = unlist(lapply(strsplit(sample.id, "_"), function(x){return(x[length(x)])}))
Size = rep(1, length(temp))
Size[temp=="HC"]=2
# 2) plot
plot(tab$EV2, tab$EV1, col=col.list[as.integer(tab$pop)],  xlab="eigenvector 2", ylab="eigenvector 1", pch = (as.integer(tab$pop)%%25), cex = Size)


#### Second Plot (ALL) - with reprojection ####
indv <- unique(pop_code)
sample.HC <- paste0(indv, "_HC")

#PCA only for HC
ccm_pca_HC <- snpgdsPCA(genofile, sample.id=sample.HC, autosome.only=FALSE, remove.monosnp=FALSE)

tab.HC <- data.frame(sample.HC = (ccm_pca_HC$sample.id %in% sample.HC), pop = factor(sample.HC), 
                  EV1 = ccm_pca_HC$eigenvect[,1],    # the first eigenvector
                  EV2 = ccm_pca_HC$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
# 1) Get only HC
tab.HC <- tab.HC %>% filter(tab.HC$sample.HC=='TRUE') 

# 2) Plot only HC
# set a color palette
library(RColorBrewer)
col.list=sample(color, 32)
#plot HC only 
PCA_HC <- plot(tab.HC$EV2, tab.HC$EV1, col=col.list[as.integer(tab.HC$pop)],  xlab="eigenvector 2", ylab="eigenvector 1", pch = (as.integer(tab.HC$pop)%%25), cex = 2)

# Get SNP loads from HC 
SNPsLoad = snpgdsPCASNPLoading(ccm_pca_HC, genofile)

dim(SNPsLoad$snploading) #gets you how many samples and SNPs you have

# Reproject the first PCA values into all individuals
ReprojectionAll = snpgdsPCASampLoading(SNPsLoad, genofile)

plot(ReprojectionAll$eigenvect[,1], ReprojectionAll$eigenvect[,2], pch = 20)
points(test$eigenvect[,1], test$eigenvect[,2], col = "red")


tab <- data.frame(sample.id = ReprojectionAll$sample.id, pop = factor(pop_code)[match(ReprojectionAll$sample.id, sample.id)], 
                  EV1 = ReprojectionAll$eigenvect[,1],    # the first eigenvector
                  EV2 = ReprojectionAll$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)


#Generate colors
library(RColorBrewer)
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
col.list=sample(color, 32)
pie(rep(1,32), col=col.list)
#Plot
# 1) change size of HC
temp = unlist(lapply(strsplit(ReprojectionAll$sample.id, "_"), function(x){return(x[length(x)])}))
Type = rep(1, length(temp))
Type[temp=="HC"]=2
Size = rep(.5, length(temp))
Size[temp=="HC"]=2
# 2) plot
pdf("1_Reprojected_PCA_Super-Scaffold_14_v1.pdf")
plot(ReprojectionAll$eigenvect[,1], ReprojectionAll$eigenvect[,2], col=col.list[as.integer(tab$pop)],  xlab="eigenvector 2", ylab="eigenvector 1", pch = (as.integer(tab$pop)%%25), cex = Size)
#abline(v = 0)
#abline(h = 0)
dev.off()


library(ade4)
pdf("1_Reprojected_PCA_Super-Scaffold_14_v2.pdf")
s.chull(ReprojectionAll$eigenvect[,c(1,2)], fac = tab$pop, label = NA)
points(ReprojectionAll$eigenvect[,c(1,2)], pch = c(20, 0)[Type], cex = Size/2)
dev.off()
