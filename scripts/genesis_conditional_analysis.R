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
  make_option(c("-a", "--chromosome"), help = "chromosome of the specific marker"),
  make_option(c("-b", "--bp"), help = "base-pair position of specific marker"),
  make_option(c("-d", "--distance"), help = "distance on either side of the position in which to perform conditional analysis"),
  make_option(c("-m", "--marker"), help = "genotype data for specific marker to be used as covariate in conditional analysis.  Expected to be output via --recodeAD option in PLINK"),
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
  make_option(c("-c", "--cores"), type="integer", default=1,
              help="Number of cores to use for parallel processing")
)

opt <- parse_args(OptionParser(option_list=option_list))
registerDoParallel(cores = as.numeric(opt$cores))

mdat = read.table(opt$markers, header = T)

chrom = as.numeric(opt$chromosome)
pos = as.numeric(opt$bp)
dist = as.numeric(opt$distance)

pdat = opt$pdat

grm = openfn.gds(opt$grm)

gdsFile = opt$gds

Covars = c(paste0("PC", 1:as.numeric(opt$pc_number)), 'marker')
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

########################## Main Body ######################################

geno <- GdsGenotypeReader(filename = gdsFile, allow.fork = T)

chroms = getChromosome(geno)
positions = getPosition(geno)
snpIDs = getSnpID(geno, (chroms==chrom & positions > (pos - dist) & positions < (pos + dist) & positions != pos))
gdsSubset(gdsFile, paste0(outName, ".genesis.gds"), snp.include=snpIDs)
gdsSubsetCheck(gdsFile, paste0(outName, ".genesis.gds"), snp.include=snpIDs)

close(geno)

geno <- GdsGenotypeReader(filename = paste0(outName, ".genesis.gds"), allow.fork = T)
genoData <- GenotypeData(geno)

pc.df = read.table(pdat, header = T)

print(head(mdat))
print(head(pc.df))

mdat$A = mdat[,7]
mdat$D = mdat[,8]

mdat %<>% mutate(scanID = V2) %>% select(scanID, A, D)

pd.df %<>% left_join(mdat, by = "scanID") 

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
genoIterator <- GenotypeBlockIterator(genoData, snpBlock=snpBlockSize)
results = assocTestSingle(genoIterator, null.model = nullmod, test = Test)

# Output results
write.table(results, outName, quote=F, row.names = F)



