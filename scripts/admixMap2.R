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
  make_option(c("-d", "--inputDirectory"), 
              help="Directory containing output of RFMix"),
  make_option(c("-o", "--outPrefix"), default="AdmixMap", 
              help = "output prefix"),
  make_option(c("-s", "--samplesFile"), default="all", 
              help="path to file containing samples to include"),
  make_option(c("-p", "--phenoFile"), metavar="additional covariates", default = "none",
              help="file containing additional covariates to include in the model.  Must have a column named 'sample' that corresponds to samples in -f"),
  make_option(c("-P", "--Predictors"), default = "sex,AFR,AFR.glob",
              help=""),
  make_option(c("-c", "--cores"), type="integer", default=1,
              help="Number of cores to use for parallel processing"),
  make_option(c("-m", "--mode"), default="both",
              help="modes are: both, stats, or munge.  if 'munge' then data will be prepared and written, but no stats will be calculated.  'stats' can be used if 'munge' has been previously run (e.g. if munging multiple datasets and then running stats on combined data")  
)

opt <- parse_args(OptionParser(option_list=option_list))
registerDoParallel(cores = as.numeric(opt$cores))
predictors = strsplit(opt$Predictors, ",")[[1]]
######################### END: Read in arguments #########################

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
  nrow(dat)
  DF = rbind(DF, dat)
}
if (samples != "all"){
  samples = read.table(samples)
  DF %<>% filter(V1 %in% samples$V1)
}
DF %<>% pivot_longer(-V1, names_to="pop", values_to="ancestry") %>% group_by(V1, pop) %>% summarize(ancestry = mean(ancestry)) %>% spread(pop,ancestry)
colnames( DF ) <- unlist(header)
return(DF)
}

readMSP = function(filename, cases=NA, filter_case = FALSE, label_case = TRUE){
  q1h = read.table(filename, header=F, comment.char = "", nrows=2, fill = TRUE, stringsAsFactors = FALSE)
  q1p = q1h[1,] %>% pivot_longer(colnames(q1h)) %>% filter(value != "") %>% separate(value, c("Pop","anc")) %>% filter(row_number()>2)
  q1h = q1h[2,-6]
  q1 = read.table(filename, header=F, comment.char = "", skip = 2)
  colnames( q1 ) <- unlist(q1h)
  q1 %<>% pivot_longer(-c(`#chm`,spos,epos,sgpos,egpos,snps), names_to="sample", values_to="anc")
  q1 %<>% separate(sample, c("sample","hap"),"[.]") %>% mutate(chm = `#chm`) %>% select(-`#chm`)
  q1 %<>% group_by(chm,spos,epos,sgpos,egpos,snps, sample) %>% count(anc) %>% ungroup() %>% mutate(anc = as.character(anc)) %>% left_join(q1p, by = "anc") %>% select(-c(anc,name)) %>% pivot_wider(id_cols = c(chm,spos,epos,sgpos, egpos,snps,sample), names_from=Pop, values_from=n, values_fill = list(n=0))

  return(q1)
}

preMunge = function(fileDir, samples, fam_pheno, other_pheno){
  Q = getGlobalAnc(fileDir, samples)
  pdat = read.table(fam_pheno, comment.char = "")
  pdat %<>% mutate(sample = V2, sex = V5, pheno = V6 - 1) %>% filter(pheno >= 0) %>% left_join(Q, by="sample")
  if (other_pheno != "none"){
    odat = read.table(other_pheno, head = T)
    pdat %<>% left_join(odat, by="sample")
  }
  
  if (samples != "all"){
    samples = read.table(samples, comment.char = "")
  }
  
  files = list.files(fileDir, pattern = "\\.msp.tsv$")
  dat = tibble()
  for (f in 1:length(files)){
    in_name = paste(fileDir, files[f], sep = "")
    anc_dat = readMSP(in_name)
    anc_dat %<>% left_join(pdat, by="sample", suffix = c("",".glob"))
    if (samples != "all"){
      anc_dat %<>% filter(sample %in% samples$V1)
    }
    write.table(anc_dat, str_replace(in_name, ".msp.tsv", ".dat"), quote=F, row.names = F, col.names = T)
  }  
}

calcGLM = function(fileDir, predictors, family = "binomial", random = -9){
  files = list.files(fileDir, pattern = "\\.dat$")
  Out = tibble()
  for (f in 1:length(files)){
    dat = read.table(paste(fileDir, files[f], sep = ""), header = T, comment.char = "")
    pos = dat %>% distinct(spos)
    res = foreach(i = 1:nrow(pos), .combine=rbind) %dopar% {
      chrom = dat %>% distinct(chm)
      ndat = dat[dat$spos==pos[i,],]
      if (random==-9){
        mod = glm(reformulate(paste(predictors, collapse='+'), response = "pheno"), data = ndat, family = family)
      }
      else{
        mod = glmer(reformulate(paste(predictors, collapse='+'), response = "pheno"), data = ndat, family = family, control = glmerControl(optimizer = "bobyqa"),
            nAGQ = 10)
      }
      coeffs = tidy(mod) %>% mutate(Term = case_when(term=="(Intercept)" ~ "intercept", TRUE ~ term)) %>% select(-term) %>% pivot_longer(c(estimate, std.error,statistic, p.value)) %>% mutate(term = paste(Term, name, sep = ".")) %>% select(-c(name,Term)) %>% pivot_wider(names_from = term)
      tibble(chm = ndat[1,1], spos = pos[i,], epos = ndat[1,3], sgpos = ndat[1,4], egpos = ndat[1,5], snps = ndat[1,6]) %>% bind_cols(glance(mod),coeffs)
    }
    Out = bind_rows(Out, res)
  }
  return(Out)
}

######################### END: Define Functions #########################


######################### Run main components ###########################
if (opt$mode != 'stats'){
  preMunge(opt$inputDirectory, opt$samplesFile, opt$famFile, opt$phenoFile)
  # dat = bind_rows(dat, q1)
}
if (opt$mode != 'munge'){
res = calcGLM(opt$inputDirectory, predictors)
write.table(res, paste(opt$outPrefix, 'admixmap.txt', sep = ""), quote = F, row.names = F)
}


