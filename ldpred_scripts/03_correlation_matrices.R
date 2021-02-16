library(bigsnpr)
ukb <- snp_attach("data/UKBB_imp_HM3.rds")
G <- ukb$genotypes
CHR <- as.integer(ukb$map$chromosome)
POS <- ukb$map$physical.pos

(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)

for (chr in 1:22) {
  
  ind.chr <- which(CHR == chr)
  #ind.row=ind.val
  snp_cor(G, ind.col = ind.chr,
          alpha = 1, infos.pos = POS2[ind.chr], size = 3 / 1000,
          ncores = NCORES
    ),
    file = paste0("data/corr/chr", chr, ".rds")
  )
}

for (chr in 1:22) {  
  # Indices in 'sumstats'
  ind.chr <- which(sumstats$chr == chr)
  # indices in 'G' / rsid replacing `NUM_ID`
  ind.chr2 <- sumstats$rsid[ind.chr]
  # Indices in 'corr'
  ind.chr3 <- match(ind.chr2,which(CHR ==chr))
  
  corr0 <- readRDS(paste0("data/corr/chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    df_beta <- sumstats[ind.chr, c("beta", "beta_se", "n_eff", "rsid")]
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    df_beta <- rbind(df_beta, sumstats[ind.chr, c("beta", "beta_se", "n_eff", "rsid")])
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
  
}


(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))
h2_est <- ldsc[["h2"]]

#LDpred2 - auto
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.9, 30),
                               ncores = NCORES)

beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
pred_auto <- big_prodMat(G, beta_auto, ind.row = 1:nrow(G),
                           ind.col = df_beta[["rsid"]])
sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- rowMeans(beta_auto[, keep])

pred_test <- big_prodMat(G, final_beta_auto, ind.row = 1:nrow(G),
                           ind.col = df_beta[["rsid"]])


