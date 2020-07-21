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
  make_option(c("-i", "--inputFile"), 
              help="original gwas results"),
  make_option(c("-s", "--snpEff"), help = "parsed snpEff annotation"),
  make_option(c("-r", "--rsIDs"), help="text file with marker name and rsID"),
  make_option(c("-f", "--freqs"), help="case-control frequencies from plink"),
  make_option(c("-o", "--outPre"), help="output prefix")
)

opt <- parse_args(OptionParser(option_list=option_list))

dat = read.table(opt$inputFile, header = T)

freq = read.table(opt$freqs, header = T)
freq %<>% separate(SNP,into=c("chr2","pos2")) %>% mutate(chr = as.numeric(chr2), pos = as.numeric(pos2)) %>% select(-c("CHR","chr2","pos2"))

rs = read.table(opt$rsIDs)
rs %<>% mutate(chr = V1, pos = V2, rsID = V3, ref = V4, alt = V5) %>% select(-c(V1,V2,V3,V4,V5))
print(head(rs))

annt = read.table(opt$snpEff, fill=T)
annt %<>% mutate(chr = V1, pos = V2, ref = V3, alt = V4, gene = V5, severity = V6, type = V7, ann.allele = V8) %>% select(-c(V1,V2,V3,V4,V5,V6,V7,V8))
print(head(annt))

soft = str_split(opt$inputFile,"[.]")[[1]]
soft = soft[length(soft) - 1]
print(soft)

if (soft=="blink"){
  print(opt$outPre)
  annt %>% left_join(dat, by = c("chr","pos")) %>% left_join(rs, by = c("chr","pos"), suffix = c(".annt","rs")) %>% left_join(freq, by = c("chr","pos")) %>% arrange(p_value) %>% select(chr, pos, rsID, gene, type, p_value, maf, MAF_A, MAF_U, NCHROBS_A, NCHROBS_U, ref.annt, alt.annt, refrs, altrs, A1, A2) %>% write.table(opt$outPre, row.names = F, quote = F)
  # annt %<>% left_join(dat, by = c("chr","pos")) 
  # annt %<>% left_join(rs, by = c("chr","pos"), suffix = c(".annt","rs")) 
  # annt %<>% left_join(freq, by = c("chr","pos")) 
  # annt %>% arrange(p_value) %>% select(chr, pos, rsID, gene, type, p_value, maf, MAF_A, MAF_U, NCHROBS_A, NCHROBS_U, ref.annt, alt.annt, refrs, altrs, A1, A2) %>% write.table(opt$outPre, row.names = F, quote = F)
} 

if (soft=="genesis"){
  annt %>% left_join(dat, by = c("chr","pos")) %>% left_join(rs, by = c("chr","pos"), suffix = c(".annt","rs")) %>% mutate(OR = exp(Est)) %>% left_join(freq, by = c("chr","pos")) %>% arrange(SPA.pval) %>% select(chr, pos, rsID, gene, type, SPA.pval, OR, Est, freq, MAC, MAF_A, MAF_U, NCHROBS_A, NCHROBS_U, A1, A2) %>% write.table(opt$outPre, row.names = F, quote = F)
}