adm = read.table("~/Documents/Research/ALL/data/admiral_impt2_AA.admixmap.txt", head = T)
cp = read.table("~/Documents/Research/ALL/data/AncInf.admixmap.txt", head=T)

dat = foreach(i=1:22,.combine="rbind") %do% {
  dat1 = adm %>% filter(chm==i) %>% distinct(spos,epos)
  dat2 = cp %>% filter(chm==i) %>% distinct(spos,epos)
  i1 = IRanges(dat1$spos,dat1$epos)
  i2 = IRanges(dat2$spos,dat2$epos)
  p = findOverlapPairs(i1,i2)
  t = pintersect(p)
  int = cbind(data.frame(first=p@first,second=p@second),data.frame(start = t@start,end = t@start + t@width))
  int1 = int %>% mutate(spos = first.start, epos = first.end, int.start = start, int.end = end) %>% select(spos, epos, int.start, int.end)
  int2 = int %>% mutate(spos = second.start, epos = second.end, int.start = start, int.end = end) %>% select(spos, epos, int.start, int.end)
  ndat1 = adm %>% filter(chm==i) %>% right_join(int1,by=c("spos","epos"))
  ndat2 = cp %>% filter(chm==i) %>% right_join(int2,by=c("spos","epos"))
  ndat1 %>% inner_join(ndat2,by=c("chm","int.start","int.end","anc"), suffix=c(".umn",".yale")) %>% mutate(pos = ceiling((int.start + int.end)/2)) %>% mutate(SNP = paste(chm,pos,sep=":")) %>% filter(int.end - int.start > 2)
}


dat %>% filter(anc=="EUR") %>% mutate(REF_ALLELE="A",ALT_ALLELE="G",BETA=ANC.estimate.yale,PVALUE=ANC.p.value.yale,N=df.null.yale) %>% select(SNP,REF_ALLELE,ALT_ALLELE,BETA,PVALUE,N) %>% write.table("~/Documents/Research/ALL/data/Yale_admixmap_EUR_metal.txt", row.names = F,quote=F)

dat %>% filter(anc=="AFR") %>% mutate(REF_ALLELE="A",ALT_ALLELE="G",BETA=ANC.estimate.yale,PVALUE=ANC.p.value.yale,N=df.null.yale) %>% select(SNP,REF_ALLELE,ALT_ALLELE,BETA,PVALUE,N) %>% write.table("~/Documents/Research/ALL/data/Yale_admixmap_AFR_metal.txt", row.names = F,quote=F)

dat %>% filter(anc=="EUR") %>% mutate(REF_ALLELE="A",ALT_ALLELE="G",BETA=ANC.estimate.umn,PVALUE=ANC.p.value.umn,N=df.null.umn) %>% select(SNP,REF_ALLELE,ALT_ALLELE,BETA,PVALUE,N) %>% write.table("~/Documents/Research/ALL/data/UMN_admixmap_EUR_metal.txt", row.names = F,quote=F)

dat %>% filter(anc=="AFR") %>% mutate(REF_ALLELE="A",ALT_ALLELE="G",BETA=ANC.estimate.umn,PVALUE=ANC.p.value.umn,N=df.null.umn) %>% select(SNP,REF_ALLELE,ALT_ALLELE,BETA,PVALUE,N) %>% write.table("~/Documents/Research/ALL/data/UMN_admixmap_AFR_metal.txt", row.names = F,quote=F)