.libPaths(c("/usr/local/lib/R/site-library", .libPaths())) # This is necessary so that singularity doesn't look in bound home directory for libraries.  The pre-pended directory ought to be the one WITHIN the singularity image.
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(optparse))

######################### Read in arguments #########################

option_list <- list( 
  make_option(c("-f", "--fam_file"),
              help="PLINK fam file containing all of the SNP IDs"),
  make_option(c("-o", "--output"), help="output_file"),
  
  #Make option so that it reads in trailing arguments...or just do comma separated list
  make_option(c("-d", "--num_latent_factors"),
              help="Number of latent factors to infer, including the intercept"))

opt <- parse_args(OptionParser(option_list=option_list))

######################### END: Read in arguments #########################

res = foreach(i = 1:length(inputs), .combine = rbind) %do% {
p = read.table(input)
counts = hist(p, bins)$counts[-1] # Get counts in bins and remove 0 bin
props = counts/sum(counts)
tibble("file"=inputs[i], "entropy" = Entropy(props))
}

winner = res %>% filter(entropy = min(res$entropy)) %>% select(file)

snps = read.table(winner$file)

IDs = read.table(bimfile)