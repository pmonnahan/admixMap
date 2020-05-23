.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(wrapr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(testit))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(GWASTools))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(SNPRelate))

######################### Read in arguments #########################

option_list <- list( 
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false", 
              dest="verbose", help="Print little output"),
  make_option(c("-p", "--plink_prefix"),
              help="plink prefix for bed, bim, and fam file"),
  make_option(c("-c", "--pca_file"), 
              help="pca file as produced by plink"),
  make_option(c("-n", "--pc_number"), 
              help="number of pc's to include as covariates"),
  make_option(c("-O", "--outDir"), help = "output directory"),
  make_option(c("-o", "--outPrefix"), help = "output prefix"),
  make_option(c("-s", "--samplesFile"), default="all", 
              help="path to file containing samples to include"),
  make_option(c("-C", "--phenoFile"), metavar="additional covariates", default = "none",
              help="file containing additional covariates to include in the model.  Must have column labelled with ID specifying sample IDs in a way that will match what is seen in the plink fam file"),
  make_option(c("-c", "--cores"), type="integer", default=1,
              help="Number of cores to use for parallel processing")
)

opt <- parse_args(OptionParser(option_list=option_list))
registerDoParallel(cores = as.numeric(opt$cores))

#Load/convert plink data
snpgdsBED2GDS(bed.fn = paste(opt$p, ".bed", sep = ""), 
              bim.fn = paste(opt$p, ".bim", sep = ""), 
              fam.fn = paste(opt$p, ".fam", sep = ""), 
              out.gdsfn = paste(opt$O, opt$o, ".gds", sep = ""))


geno <- GdsGenotypeReader(filename = paste(opt$O, opt$o, ".gds", sep = ""))
genoData <- GenotypeData(geno)



#Load phenotype and covariate data
mydat = read.table(paste(opt$p, ".fam", sep = ""))
mydat %<>% mutate(ID = V2, sex=V5,pheno=V6) %>% select(ID,sex,pheno)
# pca = read.table(pca_file, comment.char = "")
# pca = pca[,1:(as.numeric(opt$n)+2)]
# pca %<>% mutate(ID = V2) %>% select(-V2)
# mydat %<>% left_join(pca, by = c("ID"))

# Not tested
if (opt$C != "none"){
  add_covar = read.table(opt$C, comment.char = "")
  mydat %<>% left_join(add_covar, by ="ID")
  covar_names = add_covar %>% select(-ID) %>% colnames()
}

#Convert pheno/covariate data to ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(mydat)

#Format covariate string
Covars = c("sex") + covar_names
for (i in 3:ncol(pca)){
  varName = paste("PC",i - 2, sep="")
  scanAnnot[[varName]] <- pca[,i]
  Covars = c(Covars, paste("PC",i - 2, sep=""))
}

# Generate kinship matrix STOPPED HERE
# mypcrel contains Kinship Estimates from a previous PC-Relate analysis
myGRM <- pcrelateToMatrix(myp crel)
myGRM[1:5,1:5]

# fit the null mixed model
nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = Covars, 
                        cov.mat = myGRM, family = binomial)

# Perform GWAS...this needs to be parallelized by chromosome
genoIterator <- GenotypeBlockIterator(HapMap_genoData, snpBlock=5000)
assoc <- assocTestSingle(genoIterator, null.model = nullmod)
