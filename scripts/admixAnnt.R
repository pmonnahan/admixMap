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

######################### Read in arguments #########################

option_list <- list( 
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false", 
              dest="verbose", help="Print little output"),
  make_option(c("-i", "--input_directory"), metavar="plink .fam file",
              help="directory containing the .dat files generated via admixMunge.R"),
  make_option(c("-p", "--preanntFile"), 
              help="result of bedtools intersect b/n significant admixture mapping regions and GFF file"),
  make_option(c("-b", "--buffer"), 
              help="amount of base pairs on either side that the regions were extended"),
  make_option(c("-o", "--outFile"), help = "output file name")
)

opt <- parse_args(OptionParser(option_list=option_list))


rslts_dir = "~/Documents/Research/osteoCNV/data/"
fam_file = ""

files = list.files(opt$input_directory, pattern = '\\.dat')

buffer = as.numeric(opt$buffer)

preannt = read.table(opt$preanntFile, header=F)
preannt %<>% mutate(chm = V1, spos = V2 + buffer, epos = V3 - buffer, ANC.est = V4, ANC.se = V5, ANC.stat = V6, ANC.pval = V7, gene.start = V11, gene.stop = V12, gene = V16) %>% select(c(chm,spos,epos,ANC.est,ANC.se,ANC.stat,ANC.pval,gene.start,gene.stop,gene))

mungeDat = function(files){
  Dat = foreach(i = 1:length(files), .combine=rbind) %do% {
    file = files[i]
    dat = read.table(paste(opt$input_directory,file,sep="/"), header=T)
    dat %<>% group_by(chm,spos,epos,pheno) %>% summarize_all(mean) %>% mutate(Pheno=case_when(pheno==1~'case', pheno==0~'ctrl')) %>% select(-c(sample, snps, sgpos,egpos,sex,pheno)) %>% pivot_longer(-c(chm,spos,epos,Pheno)) %>% filter(!grepl('glob',name)) %>% mutate(avg.anc = value / 2) %>% select(-value) %>% pivot_wider(id_cols=c(chm,spos,epos), names_from=c('name','Pheno'),values_from='avg.anc')
    dat %<>% inner_join(preannt, by = c("chm","spos","epos"))
    dat
  }
  return(Dat)
}

aDat = mungeDat(files)
write.table(aDat, opt$outFile, row.names=F, quote=F)
