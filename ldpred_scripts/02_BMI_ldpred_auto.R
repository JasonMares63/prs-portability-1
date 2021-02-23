.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
library(bigsnpr)
library(Matrix)
library(tidyverse)


sumstats <- read_tsv("data/LDpred2/BMI_merged_sumstats_ldpred.tsv")
tmp <- tempfile(tmpdir = "data/LDpred2/tmp-data")
NCORES <- 1

for (chr in 1:22) {
  
  rds_file <- paste0("data/LDpred2/LD_EUR_train_",chr,".rds")
  ukb <- snp_attach(rds_file)
  G <- ukb$genotypes
  CHR <- as.integer(ukb$map$chromosome)
  POS <- ukb$map$physical.pos
  
  
  # Indices in 'sumstats'
  ind.chr <- which(sumstats$chr == chr)
  # indices in 'G' / rsid replacing `NUM_ID`
  ind.chr2 <- sumstats$rsid[ind.chr]
  # Indices in 'corr' / ind.chr replacing which(CHR == chr) 
  ind.chr3 <- match(ind.chr2,sumstats$rsid[ind.chr])
  
  corr0 <- readRDS(paste0("data/LDpred2/corr_chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    df_beta <- sumstats[ind.chr, c("chr","beta", "beta_se", "n_eff", "rsid")]
    ld <- colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    df_beta <- rbind(df_beta, sumstats[ind.chr, c("chr","beta", "beta_se", "n_eff", "rsid")])
    ld <- c(ld, colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
  
}


(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))
h2_est <- ldsc[["h2"]]

print("LDpred2 - auto taking place")
#LDpred2 - auto
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),ncores=NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

beta_auto_df <- as.data.frame(beta_auto)
final_beta_auto <- rowMeans(beta_auto_df)

sumstats1 <- sumstats %>%
  cbind(beta_auto_df) %>%
  mutate(avg_ldpred_beta=final_beta_auto) %>%
  select(-a0,-a1,-beta_se,-p,-n_eff)

sumstats1 %>% to_tsv("data/LDpred2/BMI_ldpred_auto_results.tsv")
