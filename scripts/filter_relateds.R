suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))

option_list <- list( 
  make_option(c("-r", "--relateds"),
              help="Output of PLINK IBS calculation; i.e. XXX.genome. Filtered for related individuals"),
  make_option(c("-o", "--output"), help="output_file"))

opt <- parse_args(OptionParser(option_list=option_list))

rel = read.table(opt$relateds, head = T)

rel %<>% mutate(ID1=paste(FID1,IID1,sep="|"), ID2=paste(FID2,IID2,sep="|"))
multi_rel = rel %>% select(ID1,ID2) %>% gather() %>% count(value) %>% filter(n>1)
sing_rel = rel %>% filter(!(ID1 %in% multi_rel$value) & !(ID2 %in% multi_rel$value))
sing_exclude = sing_rel %>% mutate(choice = rbinom(n(),1,0.5)) %>% mutate(Choice=case_when(choice==1~as.character(ID1), TRUE~as.character(ID2))) %>% select(Choice)
multi_exclude = multi_rel %>% mutate(Choice=value) %>% select(Choice)

exclude = rbind(sing_exclude,multi_exclude)
exclude %<>% separate(Choice, c("FID","IID"), "|") %>% select(FID,IID)

write.table(exclude, opt$output, quote = F, row.names = F, col.names = F)

