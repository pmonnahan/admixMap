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

######################### Read in arguments #########################

option_list <- list( 
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false", 
              dest="verbose", help="Print little output"),
  make_option(c("-f", "--famFile"), metavar="plink .fam file",
              help="plink .fam file containing case/control or phenotype values and sex"),
  make_option(c("-c", "--covariates"), default = "none",
              help="file containing other covariate data"),
  make_option(c("-o", "--outPrefix"), default="BLINK/", 
              help = "output prefix"),
  make_option(c("-d", "--PCAfile"), 
              # help=".pcrelate.txt from genesis_gwas.R script"),
              help=".eigenvec file from plink"),
  make_option(c("-n", "--PCAnumber"), default = 4,
              help="number of PCs you wish to retain")
)

opt <- parse_args(OptionParser(option_list=option_list))

######################### END: Read in arguments #########################
fam = read.table(opt$famFile, comment.char="", na.strings = "-9")

fam %<>% mutate(taxa = V2, phenotype = V6, sex = V5)

fam %>% select(taxa, phenotype) %>% write.table(paste0(opt$outPrefix, ".txt"), row.names = F, quote = F)

## For using PLINK PCs
pca = read.table(opt$PCAfile, comment.char="")
pca_num = as.numeric(opt$PCAnumber) + 2
pca %<>% select(2:all_of(pca_num)) %>% mutate(taxa=str_replace(V2,"#","")) %>% select(-V2)

## For when using PCrelate PCs
# pca = read.table(opt$PCAfile, comment.char="", header=T)
# pca %<>% mutate(taxa=str_replace(scanID,"#","")) %>% select(-scanID)

if (opt$c != "none"){
  cov = read.table(opt$covariates, comment.char="")
  cov %<>% mutate(taxa=str_replace(V2,"#",""))
  fam %<>% select(taxa,sex) %>% left_join(cov, by = "taxa") %>% left_join(pca, by = "taxa")
  fam %>% write.table(paste0(opt$outPrefix, ".cov"), row.names = F, quote = F)
  
} else {
  fam %<>% select(taxa,sex) %>% left_join(pca, by = "taxa")
  fam %>% write.table(paste0(opt$outPrefix, ".cov"), row.names = F, quote = F)
}