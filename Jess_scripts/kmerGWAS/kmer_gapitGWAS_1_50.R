# GAPIT ALL
if(!require(gplots)) install.packages("gplots", repos='http://cran.us.r-project.org')
if(!require(LDheatmap)) install.packages("LDheatmap", repos='http://cran.us.r-project.org')
if(!require(genetics)) install.packages("genetics", repos='http://cran.us.r-project.org')
if(!require(ape)) install.packages("ape", repos='http://cran.us.r-project.org')
if(!require(compiler)) install.packages("compiler", repos='http://cran.us.r-project.org')
if(!require(grid)) install.packages("grid", repos='http://cran.us.r-project.org')
if(!require(bigmemory)) install.packages("bigmemory", repos='http://cran.us.r-project.org')
if(!require(EMMREML)) install.packages("EMMREML", repos='http://cran.us.r-project.org')
if(!require(scatterplot3d)) install.packages("scatterplot3d", repos='http://cran.us.r-project.org')
if(!require(lme4)) install.packages("lme4", repos='http://cran.us.r-project.org')

if(!'multtest'%in% installed.packages()[,"Package"]){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("multtest", force=TRUE)
  BiocManager::install("snpStats")
}

source("/cluster/scratch/stapleyj/gwas_kmer/gapit_functions.txt")  
source("/cluster/scratch/stapleyj/gwas_kmer/emma.txt")    
setwd("/cluster/scratch/stapleyj/gwas_kmer/gapit_out/res1")

file.list <- scan("/cluster/scratch/stapleyj/gwas_kmer/hmp_all.list", what="character")

myY <- read.table("/cluster/scratch/stapleyj/gwas_kmer/halo_pheno24_all.txt", head = TRUE)

for (i in 2:50) {
  nm <- file.list[i]
  setwd(paste0("../res",i))
  myG <- read.table(paste0("../../hmp_files_all/hmp_chrN_files/",nm), head =FALSE)
  gwas.halo <- GAPIT(Y=myY[,c(1,2)], G=myG, SNP.MAF=0.05, PCA.total = 2,
                     Geno.View.output = FALSE, model="MLM")
}

q()
