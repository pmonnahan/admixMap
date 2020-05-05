.libPaths(c("/usr/local/lib/R/site-library", .libPaths())) # This is necessary so that singularity doesn't look in bound home directory for libraries.  The pre-pended directory ought to be the one WITHIN the singularity image.
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(optparse))

######################### Read in arguments #########################

option_list <- list( 
  make_option(c("-p", "--plink_prefix"),
              help="Plink binary prefix"),
  make_option(c("-o", "--output"), help="output_file"),
  make_option(c("-i", "--input_sHWE_runs"),
              help="comma-separated list (no spaces) with the outputs of sHWE"),
  make_option(c("-t", "--threshold"), default = 0.001,
              help="p-value threshold for calling a SNP as not in HWE"))

opt <- parse_args(OptionParser(option_list=option_list))

inputs = strsplit(opt$i, ",")[[1]]
bins = as.numeric(opt$bins)
bim = paste0(opt$p,".bim")
threshold = as.numeric(opt$t)

######################### END: Read in arguments #########################
print(inputs)
res = foreach(i = 1:length(inputs), .combine = rbind) %do% {
  snps = read.table(winner$file)
  IDs = read.table(bim)
}


print(winner$file)
IDs = read.table(bim)

IDs = IDs[snps$V1 < threshold, ]

write.table(IDs$V2, opt$o, quote=F,row.names = F, col.names = F)

