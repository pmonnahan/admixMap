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

#Parse arguments
opt <- parse_args(OptionParser(option_list=option_list))
opts <- list(chunkSize=2)

#Set up multiple cores
registerDoParallel(cores = as.numeric(opt$cores) - 1)


ancpre = strsplit(opt$AncestryPredictors, ",")[[1]]
othpre = strsplit(opt$OtherPredictors, ",")[[1]]


######################### END: Read in arguments #########################

######################### Define functions #########################

#Run General(ized) Linear model for a specific window
calcGLM = function(ndat, predictors, family = "binomial", random = -9){
  if (random==-9){ #Run GLM
    mod = glm(reformulate(paste(predictors, collapse='+'), response = "pheno"), data = ndat, family = family)
  }
  else{ #Use REML if we have random factors (not tested)
    mod = glmer(reformulate(paste(predictors, collapse='+'), response = "pheno"), data = ndat, family = family, control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
  }
  coeffs = tidy(mod) %>% mutate(Term = case_when(term=="(Intercept)" ~ "intercept", TRUE ~ term)) %>% select(-term) %>% pivot_longer(c(estimate, std.error,statistic, p.value)) %>% mutate(term = paste(Term, name, sep = ".")) %>% select(-c(name,Term)) %>% pivot_wider(names_from = term) #Parse results of GLM for tidy storage
  rslt = tibble(chm = ndat[1,1], spos = ndat[1,2], epos = ndat[1,3], sgpos = ndat[1,4], egpos = ndat[1,5], snps = ndat[1,6]) %>% bind_cols(glance(mod),coeffs) #Combine GLM results with window info
  return(rslt)
}

#Paralellize GLM calculations across windows and ancestries
wrapGLM = function(file, predictors, family = "binomial", random = -9){
    dat = read.table(file, header = T, comment.char = "")
    pos = dat %>% distinct(spos)
    glob = tibble(name = colnames(dat)) %>% filter(grepl(".glob", name)) 
    glob = glob[1:(nrow(glob) - 1),] #Retrieve all but one global ancestry values as per Grinde 2019
    pre.list = list() # A separate GLM is fit for each ancestry, so we make a list to hold sets of predictor variables
    for (a in 1:length(predictors)){ 
      pre.list[[length(pre.list)+1]] = c(predictors[a], as.vector(glob$name), othpre) #Predictors are local ancestry for current ancestry, all global ancestries (minus the last), and any other covariates.
    }
    # Loop over windows and loop over ancestries specified in 'predictors'
    res = foreach(i = 1:nrow(pos), .combine='rbind', .errorhandling="remove") %:%
      foreach(j = 1:length(pre.list), .combine = 'rbind', .errorhandling="remove") %dopar% {
      ndat = dat[dat$spos==pos[i,],] #Extract just the data for this window
      res1 = calcGLM(ndat, pre.list[[j]], family, random) #Calc GLM for this window and ancestry
      anc = pre.list[[j]][1]
      # res2 = res1 %>% rename_with( ~ str_replace(.,anc,"ANC"))
      names(res1)[19:22] = c("ANC.estimate","ANC.std.error","ANC.statistic","ANC.p.value") #Predictor estimates are labelled generically for long format
      res1 %<>% mutate(anc = anc) #Label these results with the focal ancestry
      #If the GLM failed for whatever reason, the results will not have these expected number of columns
      exp_cols = (length(pre.list[[j]]) * 4) + 18
      if (ncol(res1) == exp_cols + 1){ 
        res1 #Only store results if GLM was successfull in order to avoid conflicting rows upon binding.
      }
    }
  return(res)
}

######################### END: Define Functions #########################


######################### Run main components ###########################


res = wrapGLM(opt$inputFile, ancpre)
write.table(res, paste(opt$outPrefix, '.admixmap.txt', sep = ""), quote = F, row.names = F)



