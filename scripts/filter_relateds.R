######################### Load packages #########################

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))

######################### Read in arguments #########################

option_list <- list( 
  make_option(c("-r", "--relateds"),
              help="Output of PLINK IBS calculation; i.e. XXX.genome. Filtered for related individuals"),
  make_option(c("-o", "--output"), help="output_file"))

opt <- parse_args(OptionParser(option_list=option_list))

######################### Main Body #########################

rel = read.table(opt$relateds, header = T, comment.char = "") #Read in IBD results from PLINK

rel %<>% mutate(ID1=paste(FID1,IID1,sep="]"), ID2=paste(FID2,IID2,sep="]")) #Concatenate IDs for ease of comparison
multi_rel = rel %>% select(ID1,ID2) %>% pivot_longer(c(ID1,ID2)) %>% count(value) %>% filter(n>1) #Identify individuals with more than one related individual
sing_rel = rel %>% filter(!(ID1 %in% multi_rel$value) & !(ID2 %in% multi_rel$value)) #Identify individuals with just a single related individual
sing_exclude = sing_rel %>% mutate(choice = rbinom(n(),1,0.5)) %>% mutate(Choice=case_when(choice==1~as.character(ID1), TRUE~as.character(ID2))) %>% select(Choice) #Randomly select one of the two individuals from the single relative pairs.  
multi_exclude = multi_rel %>% mutate(Choice=value) %>% select(Choice) #Remove all individuals that are related to multiple others

exclude = rbind(sing_exclude,multi_exclude)
exclude %<>% separate(Choice, c("FID","IID"), "]") %>% select(FID,IID)

write.table(exclude, opt$output, quote = F, row.names = F, col.names = F)

