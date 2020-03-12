.libPaths(c("/usr/local/lib/R/site-library", .libPaths())) # This is necessary so that singularity doesn't look in bound home directory for libraries.  The pre-pended directory ought to be the one WITHIN the singularity image.

suppressPackageStartupMessages(library(lfa))
suppressPackageStartupMessages(library(optparse))

######################### Read in arguments #########################

option_list <- list( 
  make_option(c("-p", "--plink_prefix"),
              help="Full path of plink prefix"),
  make_option(c("-o", "--output"), help="output_file"),
  make_option(c("-d", "--num_latent_factors"),
              help="Number of latent factors to infer, including the intercept"))

opt <- parse_args(OptionParser(option_list=option_list))

######################### END: Read in arguments #########################

######################### Main Body #########################

dat = read.bed(opt$p)
dat = dat[complete.cases(dat),]
LF = lfa(dat, opt$d)
p = sHWE(dat, LF, 1)
write.table(p, opt$o, quotes=F,row.names=F,col.names = F)
