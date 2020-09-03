.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(wrapr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(broom))

######################### Read in arguments #########################

option_list <- list( 
  make_option(c("-g", "--geno"),
              help=".raw file from plink's --recodeAD"),
  make_option(c("-b", "--betas"), 
              help="GENESIS results with beta values"),
  make_option(c("-o", "--outPrefix"), default="PRS", 
              help = "output prefix")
)

opt <- parse_args(OptionParser(option_list=option_list))
g <- read.table(opt$geno, header=T)
b <- read.table(opt$betas, header=T)
out <- opt$outPrefix

b %<>% mutate(iv = 1 / (Est.SE^2))

s = sum(b$iv)

g %<>% select(-c(FID,PAT,MAT))

dat = g %>% pivot_longer(cols = -c(IID,SEX,PHENOTYPE)) %>% filter(!grepl('_HET', name)) %>% separate(name,c("id","allele"),'_') %>% mutate(variant.id = str_replace(str_replace(id,'X',''),'[.]',':')) %>% left_join(b, by='variant.id')

prs = dat %>% mutate(lr2 = (value * Est) / (Est.SE^2), lr = value * Est) %>% group_by(IID, SEX, PHENOTYPE) %>% summarize(ivPRS = sum(lr2,na.rm=T) / s, PRS = sum(lr, na.rm=T)) 

png(paste0(out, '_PRS_density.png'))
prs %>% ggplot(aes(x=PRS,fill=as.factor(PHENOTYPE))) + geom_density(alpha=0.4) + scale_fill_discrete(name='', breaks = c("1","2"), labels=c('Control', 'Case'))
dev.off()

png(paste0(out, '_ivPRS_density.png'))
prs %>% ggplot(aes(x=ivPRS,fill=as.factor(PHENOTYPE))) + geom_density(alpha=0.4) + scale_fill_discrete(name='', breaks = c("1","2"), labels=c('Control', 'Case')) + xlab("Inverse Variance Weighted PRS")
dev.off()

write.table(prs, paste0(out, '.prs.txt'), quote=F, row.names = F)
