.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
library(bigsnpr)
library(runonce)
#library(RhpcBLASctl)


######################### Prepare correlation matrices
#(NCORES <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) - 1L)
NCORES <- 1
n_val <- 1000

for (chr in 1:22) {
  print(chr)
  rds_file <- paste0("data/LDpred2/LD_EUR_train_",chr,".rds")
  ukb <- snp_attach(rds_file)
  G <- ukb$genotypes
  CHR <- as.integer(ukb$map$chromosome)
  POS <- ukb$map$physical.pos
  
  
  POS2 <- bigsnpr::snp_asGeneticPos(CHR, POS, dir = "data/LDpred2/tmp-data",
                           ncores = NCORES)
  
  
  set.seed(1)
  ind.val <- sample(nrow(G), n_val)
  #ind.test <- setdiff(rows_along(G), ind.val)
  
  
  ind.chr <- which(CHR == chr)
  save_run(snp_cor(G, ind.row=ind.val,ind.col = ind.chr,
                   alpha = 1, infos.pos = POS2[ind.chr], size = 3 / 1000,
                   ncores = NCORES),
           file = paste0("data/LDpred2/corr_chr", chr, ".rds")
  )
}

