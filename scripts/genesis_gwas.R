.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(wrapr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(testit))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(GWASTools))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(SPAtest))
suppressPackageStartupMessages(library(Biobase))

######################### Read in arguments #########################

option_list <- list( 
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false", 
              dest="verbose", help="Print little output"),
  make_option(c("-p", "--plink_prefix"),
              help="plink prefix for bed, bim, and fam file"),
  # make_option(c("-k", "--kin_file"), 
  #             help="pairwise divergence values from .kin0 from KING-robust kinship estimation"),
  make_option(c("-k", "--pruned_file"), 
              help="pruned SNPs for inclusion in PC calcs"),
  make_option(c("-n", "--pc_number"), 
              help="number of pc's to include as covariates"),
  make_option(c("-O", "--outDir"), default = getwd(), help = "output directory"),
  make_option(c("-o", "--outPrefix"), help = "output prefix"),
  make_option(c("-t", "--kin_threshold"), help = "threshold of kinship above which individuals are considered related"),
  make_option(c("-l", "--ld_threshold"), default = 0.1, help = "LD (correlation) threshold"),
  make_option(c("-b", "--snp_block_size"), default = 5000, help = "Number of SNPs to read in to memory it at a time when doing the GWAS"),
  make_option(c("-f", "--distribution_family"), default = 'binomial', help = "Error distribution family for response variable.  For case/control, set to 'binomial'"),
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

snpBlockSize = as.numeric(opt$snp_block_size)

pbim = read.table(opt$pruned_file)

if (opt$distribution_family == "binomial"){
  Test = 'Score.SPA'
} else {
  Test = 'Score'
}
########################## END Read in Arguments ############################

########################## Define functions #################################
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  print(n)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

########################## END Define functions ################################

########################## Main Body ######################################
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
# 
# KINGmat <- kingToMatrix(opt$kin_file, estimator = "Kinship", 
#                         thresh = as.numeric(opt$kin_threshold))

#LD prune the genoData
gds <- snpgdsOpen(gdsFile, allow.duplicate = TRUE)
# snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6,
#                           ld.threshold=sqrt(as.numeric(opt$ld_threshold)), verbose=TRUE, num.thread = as.numeric(opt$cores))
# pruned <- unlist(snpset, use.names=FALSE)
# length(pruned)

snpIDs = getSnpID(genoData)
pruned = which(snpIDs %in% pbim$V2)

# pruned = read.table(opt$ld_prune_snps)
KINGmat <- snpgdsIBDKING(gds, snp.id=pruned, verbose=FALSE)

kingMat <- KINGmat$kinship
dimnames(kingMat) <- list(KINGmat$sample.id, KINGmat$sample.id)

pcair = pcair(genoData, kinobj = kingMat, divobj = kingMat, snp.include = pruned)


genoIterator1 <- GenotypeBlockIterator(genoData, snpBlock = snpBlockSize, snpInclude = pruned)

mypcrelate = pcrelate(genoIterator1, pcair$vectors[,1:as.numeric(opt$pc_number),drop=FALSE])
myGRM <- pcrelateToMatrix(mypcrelate, thresh=opt$kin_threshold, scaleKin=2) # Using a threshold actually has a minimal effect here.  

#output plot of kinship
kinship <- mypcrelate$kinBtwn

png(paste0(outName,'.kinplot'))
ggplot(kinship, aes(k0, kin)) +
  geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
  geom_point(alpha=0.5) +
  ylab("kinship estimate") +
  theme_bw()
dev.off()

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

#Add PCs to annot data
pc.df <- as.data.frame(pcair$vectors)
names(pc.df) <- paste0("PC", 1:ncol(pcair$vectors))
pc.df$sample.id <- row.names(pcair$vectors)
pc.df <- pc.df %>% mutate(scanID = sample.id) %>% select(-sample.id) %>% left_join(pc.df, pData(scanAnnot), by="scanID")

# add PCs to sample annotation in SeqVarData object
annot <- AnnotatedDataFrame(pc.df)
sampleData(seqData) <- annot

# fit the null mixed model
if (opt$distribution_family == "binomial"){
  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = c(Covars, paste0("PC", 1:as.numeric(opt$pc_number))), cov.mat = myGRM, family = binomial)
} else{ 
  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = c(Covars, paste0("PC", 1:as.numeric(opt$pc_number))), cov.mat = myGRM)
  }

# Perform GWAS
chroms = getChromosome(genoData)
Chroms = unique(chroms)
results = foreach(i = 1:length(Chroms), .combine = rbind) %doPar% {
  snpIDs = getSnpID(genoData, chroms==Chroms[i])
  genoIterator <- GenotypeBlockIterator(genoData, snpBlock=snpBlockSize, snpInclude = snpIDs)
  assocTestSingle(genoIterator, null.model = nullmod, test = Test)
}

 # Score.SPA is specifically used for binomial models
write.table(results, outName, quote=F, row.names = F)

pc.df %>% select(c(scanID, paste0("PC", 1:as.numeric(opt$pc_number)))) %>% write.table(str_replace(outName, "genesis", "pcrelate"), quote=F, row.names = F)

png(paste0(outName,'.qqplot'))
if (opt$distribution_family == "binomial"){
  results %<>% filter(!is.na(SPA.pval))
  gg_qqplot(results$SPA.pval)
} else {
  results %<>% filter(!is.na(Score.pval))
  gg_qqplot(results$Score.pval)
}
dev.off()


