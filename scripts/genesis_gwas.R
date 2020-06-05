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
  make_option(c("-k", "--kin_file"), 
              help="pairwise divergence values from .kin0 from KING-robust kinship estimation"),
  make_option(c("-n", "--pc_number"), 
              help="number of pc's to include as covariates"),
  make_option(c("-O", "--outDir"), default = getwd(), help = "output directory"),
  make_option(c("-o", "--outPrefix"), help = "output prefix"),
  make_option(c("-t", "--kin_threshold"), help = "threshold of kinship above which individuals are considered related"),
  make_option(c("-l", "--ld_threshold"), default = 0.1, help = "LD (correlation) threshold"),
  make_option(c("-s", "--samplesFile"), default="all", 
              help="path to file containing samples to include"),
  make_option(c("-C", "--covarFile"), metavar="additional covariates", default = "none",
              help="file containing additional covariates to include in the model.  Must have column labelled with ID specifying sample IDs in a way that will match what is seen in the plink fam file"),
  make_option(c("-c", "--cores"), type="integer", default=1,
              help="Number of cores to use for parallel processing")
)

opt <- parse_args(OptionParser(option_list=option_list))
registerDoParallel(cores = as.numeric(opt$cores))

gdsFile = paste(opt$outDir, "/", opt$outPrefix, ".gds", sep = "")

outName = paste0(opt$outDir, "/", opt$outPrefix)

snpBlockSize = 5000

#Load/convert plink data.  Cannot contain '#' in samnple names
snpgdsBED2GDS(bed.fn = paste(opt$plink_prefix, ".bed", sep = ""), 
              bim.fn = paste(opt$plink_prefix, ".bim", sep = ""), 
              fam.fn = paste(opt$plink_prefix, ".fam", sep = ""), 
              out.gdsfn = gdsFile)

geno <- GdsGenotypeReader(filename = gdsFile)
genoData <- GenotypeData(geno)

#Load phenotype and covariate data
mydat = read.table(paste(opt$plink_prefix, ".fam", sep = ""))
mydat %<>% mutate(scanID = V2, sex=case_when(V5==1 ~ 'M', V5==2 ~ 'F', TRUE ~ NA_character_),pheno=V6 - 1) %>% select(scanID,sex,pheno)

# Generate kinship matrix 
# mypcrel contains Kinship Estimates from a previous PC-Relate analysis

# It is important that the order of individuals in the matrices kinMat and divMat match the order of individuals in the genoData

KINGmat <- kingToMatrix(opt$kin_file, estimator = "Kinship", 
                        thresh = as.numeric(opt$kin_threshold))

#LD prune the genoData
gds <- snpgdsOpen(gdsFile, allow.duplicate = TRUE)
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6,
                          ld.threshold=sqrt(as.numeric(opt$ld_threshold)), verbose=TRUE, num.thread = as.numeric(opt$cores))
pruned <- unlist(snpset, use.names=FALSE)
length(pruned)

# pruned = read.table(opt$ld_prune_snps)

pcair = pcair(genoData, kinobj = KINGmat, divobj = KINGmat, snp.include = pruned)

genoIterator1 <- GenotypeBlockIterator(genoData, snpBlock = snpBlockSize, snpInclude = pruned)

mypcrelate = pcrelate(genoIterator1, pcair$vectors[,1,drop=FALSE])
myGRM <- pcrelateToMatrix(mypcrelate)

# Not tested
Covars = c("sex")
if (opt$covarFile != "none"){
  add_covar = read.table(opt$covarFile, comment.char = "")
  mydat %<>% left_join(add_covar, by ="scanID")
  covar_names = add_covar %>% select(-scanID) %>% colnames()
  Covars = Covars + covar_names
}

#Convert pheno/covariate data to ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(mydat)
# 
# #Format covariate string
# 
# for (i in 3:ncol(pca)){
#   varName = paste("PC",i - 2, sep="")
#   scanAnnot[[varName]] <- pca[,i]
#   Covars = c(Covars, paste("PC",i - 2, sep=""))
# }

# fit the null mixed model
nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = Covars, 
                        cov.mat = myGRM, family = binomial)

# Perform GWAS...this needs to be parallelized by chromosome
genoIterator <- GenotypeBlockIterator(genoData, snpBlock=5000)
assoc <- assocTestSingle(genoIterator, null.model = nullmod)
write.table(assoc, outName, quote=F, row.names = F)





#OLD PCA stuff
# pca = read.table(pca_file, comment.char = "")
# pca = pca[,1:(as.numeric(opt$n)+2)]
# pca %<>% mutate(ID = V2) %>% select(-V2)
# mydat %<>% left_join(pca, by = c("ID"))

