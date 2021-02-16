library(tidyverse)

sumstats <- data.frame()
for (i in 1:22){
  formula <- paste0("data/gwas_results/BMI.chr",i,".BMI.glm.assoc")
  add <- read.csv(formula,sep="\t")
  sumstats <- bind_rows(sumstats,add)
  print(i)
  
}
colnames(sumstats) <- c("chr","rsid","pos","a1","test","nmiss","beta","stat","p","beta_se")
sumstats$n_eff <- 80000

sumstats %>% write_tsv("data/LDpred2/BMI_merged.glm.assoc")

###################### Format merged genetic map

bim <- data.frame()
for(i in 2:22){
  formula <- paste0("data/ukb_merged/EUR_LD_train_i",i,".bim")
  add <- read_csv(formula,sep="\t")
  bim <- bind_rows(bim,add)
  
}
colnames(bim) <- c("chr","rsid","gen_dist","pos","a0","a1")
bim %>% write_tsv("data/LDpred2/BMI_merged.glm.bim")

####################### Combined GWAS sumstats with genetic map

sumstats <- sumstats %>% inner_join(bim,by=c("chr","rsid","pos")) %>%
  select(chr,rsid,pos,a0,a1,beta,beta_se,p,n_eff)
rm(bim)
sumstats %>% write_tsv("data/LDpred2/BMI_merged_sumstats_ldpred.tsv")
