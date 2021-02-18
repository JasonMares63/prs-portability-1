
.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
library(bigsnpr)
library(tidyverse)

# Creates .bk and .rd files
for(i in 14:22){
  snp_readBed(paste0("data/LDpred2/LD_EUR_train_",i,".bed"))
}



