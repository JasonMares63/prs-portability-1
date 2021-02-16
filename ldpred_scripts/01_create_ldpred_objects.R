library(bigsnpr)
library(tidyverse)

# Creates .bk and .rd files
for(i in 1:22){
  snp_readBed(paste0("data/LDpred/EUR_LD_train_",i,".bed")
}

NCORES <- nb_cores()


