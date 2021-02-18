library(tidyverse)

phenotypes <- c("Height", "BMI" ,"RBC", "Platelet",
                 "MCV", "Monocyte", "WBC", "MCH",
                 "Eosinophil", "Lymphocyte")
for (pheno in phenotypes) {

  results_name <- paste0("data/LDpred2/",pheno,"_ldpred_auto_results.tsv")
  results <- read_tsv(results_name)

  old_results_name <- paste0("data/LDpred2/",pheno,"_merged_sumstats_ldpred.tsv")
  old_results <- read_tsv(old_results_name)


  new_results <- old_results %>% inner_join(results,by=c("chr","pos")) %>%
    select(-beta_se,-p,-rsid.y,-beta.y) %>%
    rename(rsid = rsid.x) %>%
    rename(beta = beta.x)

  new_results_name <- paste0("data/LDpred2/prs/",pheno,"_LDpred_sumstats.tsv")
  new_results %>% write_tsv(new_results_name)

}
