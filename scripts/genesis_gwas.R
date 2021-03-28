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
  make_option(c("-p", "--pdat"),
              help="pdat file from genesis_prep.R"),
  # make_option(c("-k", "--kin_file"), 
  #             help="pairwise divergence values from .kin0 from KING-robust kinship estimation"),
  make_option(c("-k", "--grm"), 
              help="genetic relatedness matrix output by genesis_prep.R"),
  make_option(c("-g", "--gds"), 
              help="GDS file containing genotype info, created by genesis_prep.R"),
  make_option(c("-n", "--pc_number"), 
              help="number of pc's to include as covariates"),
  make_option(c("-O", "--outDir"), default = getwd(), help = "output directory"),
  make_option(c("-o", "--outPrefix"), help = "output prefix"),
    make_option(c("-b", "--snp_block_size"), default = 5000, help = "Number of SNPs to read in to memory it at a time when doing the GWAS"),
  make_option(c("-f", "--distribution_family"), default = 'binomial', help = "Error distribution family for response variable.  For case/control, set to 'binomial'"),
  make_option(c("-s", "--samplesFile"), default="all", 
              help="path to file containing samples to include"),
  make_option(c("-C", "--covarString"), metavar="additional covariates", default = "none",
              help="comma delimited (no spaces) list of covariates contained in the pdat file"),
  make_option(c("-c", "--cores"), type="integer", default=2,
              help="Number of cores to use for parallel processing")
)

opt <- parse_args(OptionParser(option_list=option_list))
registerDoParallel(cores = as.numeric(opt$cores) - 1)

pdat = opt$pdat

grm = openfn.gds(opt$grm)

gdsFile = opt$gds

Covars = c(paste0("PC", 1:as.numeric(opt$pc_number)))

if (opt$covarString != 'none'){
  Covars = c(Covars, strsplit(opt$covarString, ",")[[1]])
}

outName = paste0(opt$outDir, "/", opt$outPrefix)

snpBlockSize = as.numeric(opt$snp_block_size)

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

geno <- GdsGenotypeReader(filename = gdsFile, allow.fork = T)
genoData <- GenotypeData(geno)

pc.df = read.table(pdat, header = T)
myGRM = read.gdsn(index.gdsn(grm, "kinship"))
rownames(myGRM) = read.gdsn(index.gdsn(grm, "sample.id"))
colnames(myGRM) = read.gdsn(index.gdsn(grm, "sample.id"))


# add PCs to sample annotation in SeqVarData object
scanAnnot <- ScanAnnotationDataFrame(pc.df)

# fit the null mixed model
if (opt$distribution_family == "binomial"){
  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = Covars, cov.mat = myGRM, family = binomial)
} else{ 
  nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = Covars, cov.mat = myGRM)
  }

# Perform GWAS
chroms = getChromosome(genoData)
Chroms = unique(chroms)
results = foreach(i = 1:length(Chroms), .combine = rbind) %dopar% {
  snpIDs = getSnpID(genoData, chroms==Chroms[i])
  genoIterator <- GenotypeBlockIterator(genoData, snpBlock=snpBlockSize, snpInclude = snpIDs)
  rslt = assocTestSingle(genoIterator, null.model = nullmod, test = Test)
  write.table(rslt, str_replace(outName, '.txt',paste0('.',i,'.txt')), quote=F, row.names = F)
  rslt
}

# if (opt$distribution_family == "binomial"){
#   results %<>% mutate(Score.pval=SPA.pval) %>% select(-SPA.pval)
# }

 # Score.SPA is specifically used for binomial models
write.table(results, outName, quote=F, row.names = F)


png(paste0(outName,'.qqplot'))
if (opt$distribution_family == "binomial"){
  results %<>% filter(!is.na(SPA.pval))
  gg_qqplot(results$SPA.pval)
} else {
  results %<>% filter(!is.na(Score.pval))
  gg_qqplot(results$Score.pval)
}
dev.off()


