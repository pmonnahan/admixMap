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
  make_option(c("-d", "--inputDirectory"), 
              help="Directory containing output of RFMix"),
  make_option(c("-o", "--outPrefix"), help = "output prefix"),
  make_option(c("-s", "--samplesFile"), default="all", 
              help="path to file containing samples to include")
)

opt <- parse_args(OptionParser(option_list=option_list))

######################### Define functions #########################
getGlobalAnc = function(rfmix_dir, samples){
  files <- list.files(rfmix_dir, pattern = "\\.Q$")
  DF = data.frame()
  for (i in 1:length(files)){
    if (i==1){
      header <- read.table(paste(rfmix_dir, files[i], sep="/"), nrows = 1, skip = 1, header = FALSE, comment.char = "", stringsAsFactors = FALSE)
      header[,1] = 'sample'
    }
    dat = read.table(paste(rfmix_dir, files[i], sep="/"), skip=2, header=F, comment.char = "")
    DF = rbind(DF, dat)
  }
  DF %<>% mutate(samp = str_replace(V1, "#", "")) %>% select(-V1)
  if (samples != "all"){
    samples = read.table(samples, comment.char = "")
    DF %<>% filter(samp %in% samples$V1)
  }
  DF %<>% pivot_longer(-samp, names_to="pop", values_to="ancestry") %>% group_by(samp, pop) %>% summarize(ancestry = mean(ancestry)) %>% spread(pop,ancestry)
  colnames( DF ) <- unlist(header)
  return(DF)
}

Q = getGlobalAnc(opt$inputDirectory, opt$samplesFile)
write.table(Q, opt$outPrefix, row.names=F, quote = F)