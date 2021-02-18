.libPaths("/rigel/mfplab/users/jm4454/rpackages/")

library(tidyverse)
library(varhandle)
library(bigsnpr)

# BMI done already
phenotypes <- c("Height","RBC", "Platelet", "MCV", "Monocyte", "WBC", "MCH", "Eosinophil", "Lymphocyte")

for (pheno in phenotypes){

print(pheno)

sumstats <- data.frame()
for (i in 1:22){
  formula <- paste0("data/gwas_results/GWAS_UKBB_genotype/",pheno,".chr",i,".",pheno,".glm.assoc")
  add <- read.csv(formula,sep="\t")
  sumstats <- bind_rows(sumstats,add)
  print(i)

}
colnames(sumstats) <- c("chr","rsid","pos","a1","test","nmiss","beta","stat","p","beta_se")

#sumstats %>% write_tsv("data/LDpred2/BMI_merged.glm.assoc")

###################### Format merged genetic map

bim <- data.frame()
for(i in 1:22){
  print(i)
  rds_file <- paste0("data/LDpred2/LD_EUR_train_",i,".rds")
  ukb <- snp_attach(rds_file)
  add <- ukb$map
  rm(ukb)
#  formula <- paste0("data/LDpred2/LD_EUR_train_",i,".bim")
#  add <- read.csv(formula,sep="\t")
  colnames(add) <- c("chr","rsid","gen_dist","pos","a0","a1")
  bim <- bind_rows(bim,add)
}

#bim %>% write_tsv("data/LDpred2/BMI_merged.glm.bim")

####################### Combined GWAS sumstats with genetic map

sumstats1 <- sumstats %>% right_join(bim,by=c("chr","rsid","pos")) %>%
  mutate(a0 = unfactor(a0)) %>%
  mutate(a1.x = unfactor(a1.x)) %>%
  mutate(a1.y = unfactor(a1.y))


sumstats1$a1 <- "A"
for (i in 1:nrow(sumstats1)){
  sumstats1$a1[i] <- ifelse(sumstats1$a1.x[i]==sumstats1$a0[i],
                            sumstats1$a1.y[i],
                            sumstats1$a1.x[i])
}


sumstats2 <- sumstats1 %>%
  mutate(n_eff = 80000-nmiss) %>%
  select(chr,rsid,pos,a0,a1,beta,beta_se,p,n_eff)

rm(sumstats,sumstats1,bim)

file_name <- paste0("data/LDpred2/",pheno,"_merged_sumstats_ldpred.tsv")
sumstats2 %>% write_tsv(file_name)

}

