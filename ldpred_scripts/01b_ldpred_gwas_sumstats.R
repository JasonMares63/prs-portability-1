print("02_LDpred_GWAS_sumstats")

.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
library(tidyverse)
library(bigsnpr)

phenotypes <- c("BMI","Height","RBC", "Platelet", "MCV", "Monocyte", "WBC", "MCH", "Eosinophil", "Lymphocyte")

for (pheno in phenotypes){
  
  print(pheno)
  
  sumstats1 <- data.frame()
  for (i in 1:22){
    formula <- paste0("data/gwas_results/",pheno,".chr",i,".",pheno,".glm.linear")
    add <- read.csv(formula,sep="\t")
    sumstats1 <- bind_rows(sumstats1,add)
    print(i)
    
  }
  colnames(sumstats1) <- c("chr","pos","rsid","a0","alt","a1","test","n_eff",
                           "beta","beta_se","stat","p","err_code")
  
  #sumstats %>% write_tsv("data/LDpred2/BMI_merged.glm.assoc")
  
  ###################### Format merged genetic map
  
  bim <- data.frame()
  for(i in 1:22){
    print(i)
    rds_file <- paste0("data/LDpred/LD_EUR_train_",i,".rds")
    ukb <- snp_attach(rds_file)
    add <- ukb$map
    rm(ukb)
    add <- data.frame(matrix(unlist(add), ncol=length(add), byrow=F),
                      stringsAsFactors = F)
    colnames(add) <- c("chr","rsid","gen_dist","pos","a0","a1")
    bim <- bind_rows(bim,add)
  }
  bim$chr <- as.integer(bim$chr)
  bim$pos <- as.integer(bim$pos)
  #bim %>% write_tsv("data/LDpred/BMI_merged.glm.bim")
  

  ####################### Combined GWAS sumstats with genetic map
  
  sumstats2 <- sumstats1 %>% inner_join(bim,by=c("chr","rsid","pos")) %>%
    mutate(a1.x = as.character(a1.x))
  
  sumstats2$beta <- ifelse(sumstats2$a1.x!=sumstats2$a1.y,
                           sumstats2$beta*-1,
                           sumstats2$beta)
  
  sumstats2$stat <- ifelse(sumstats2$a1.x!=sumstats2$a1.y,
                           sumstats2$stat*-1,
                           sumstats2$stat)

  
  
  sumstats2 <- sumstats2 %>%
    select(chr,rsid,pos,a0.y,a1.y,beta,beta_se,p,n_eff) %>%
    rename(a0 = a0.y) %>%
    rename(a1 = a1.y)
  
  rm(sumstats1,bim)
  
  file_name <- paste0("data/LDpred/",pheno,"_merged_sumstats_ldpred.tsv")
  sumstats2 %>% write_tsv(file_name)
  
}

