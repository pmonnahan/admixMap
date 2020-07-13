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
  make_option(c("-f", "--inputFile"), 
              help=".dat output file from admixMunge.R"),
  make_option(c("-o", "--outPrefix"), default="AdmixMap", 
              help = "output prefix"),
  make_option(c("-P", "--AncestryPredictors"), default = "AFR,EUR,AMR",
              help=""),
  make_option(c("-p", "--OtherPredictors"), default = "sex",
              help=""),
  make_option(c("-c", "--cores"), type="integer", default=1,
              help="Number of cores to use for parallel processing")
)

opt <- parse_args(OptionParser(option_list=option_list))
opts <- list(chunkSize=2)
registerDoParallel(cores = as.numeric(opt$cores) - 1)


ancpre = strsplit(opt$AncestryPredictors, ",")[[1]]
othpre = strsplit(opt$OtherPredictors, ",")[[1]]

pre.list = list()
for (a in 1:length(ancpre)){
  pre.list[[length(pre.list)+1]] = c(ancpre[a], paste0(ancpre[a],".glob"), othpre)
}

######################### END: Read in arguments #########################

######################### Define functions #########################
calcGLM = function(ndat, predictors, family = "binomial", random = -9){
  if (random==-9){
    mod = glm(reformulate(paste(predictors, collapse='+'), response = "pheno"), data = ndat, family = family)
  }
  else{
    mod = glmer(reformulate(paste(predictors, collapse='+'), response = "pheno"), data = ndat, family = family, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
  }
  coeffs = tidy(mod) %>% mutate(Term = case_when(term=="(Intercept)" ~ "intercept", TRUE ~ term)) %>% select(-term) %>% pivot_longer(c(estimate, std.error,statistic, p.value)) %>% mutate(term = paste(Term, name, sep = ".")) %>% select(-c(name,Term)) %>% pivot_wider(names_from = term)
  rslt = tibble(chm = ndat[1,1], spos = ndat[1,2], epos = ndat[1,3], sgpos = ndat[1,4], egpos = ndat[1,5], snps = ndat[1,6]) %>% bind_cols(glance(mod),coeffs)
  return(rslt)
}

wrapGLM = function(file, predictors, family = "binomial", random = -9){
    dat = read.table(file, header = T, comment.char = "")
    pos = dat %>% distinct(spos)
    res = foreach(i = 1:nrow(pos), .combine='rbind', .errorhandling="remove") %:%
      foreach(j = 1:length(predictors), .combine = 'rbind', .errorhandling="remove") %dopar% {
      ndat = dat[dat$spos==pos[i,],]
      res1 = calcGLM(ndat, predictors[[j]], family, random)
      anc = predictors[[j]][1]
      exp_cols = (length(predictors[[j]]) * 4) + 18
      res2 = res1 %>% rename_with( ~ str_replace(.,anc,"ANC"))
      res2 %<>% mutate(anc = anc)
      if (ncol(res2) == exp_cols){
        res2
      }
    }
  return(res)
}

######################### END: Define Functions #########################


######################### Run main components ###########################


res = wrapGLM(opt$inputFile, pre.list)
write.table(res, paste(opt$outPrefix, '.admixmap.txt', sep = ""), quote = F, row.names = F)



