.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
library(tidyverse)


fst_df <- read_tsv("data/fst/fst0.est", col_names="Fst",col_types="c")
fst_df$Fst <- fst_df$Fst %>% parse_number()
fst_df$label <- rep(c("IID","Mean_Fst","Weighted_Fst"),50)

for (i in 1:656){
file <- paste0("data/fst/fst",i,".est")
fst_df_append <- read_tsv(file, col_names="Fst",col_types="c")
fst_df_append$Fst <- fst_df_append$Fst %>% parse_number()
fst_df_append$label <- rep(c("IID","Mean_Fst","Weighted_Fst"),nrow(fst_df_append)/3)

fst_df <- rbind(fst_df,fst_df_append)
}


IID <- fst_df %>% filter(label=="IID") %>% select(-label) 
mean_fst <- fst_df %>% filter(label=="Mean_Fst") %>% select(-label)
w_fst <- fst_df %>% filter(label=="Weighted_Fst") %>% select(-label)

rm(fst_df)

fst <- cbind(IID,mean_fst,w_fst)
colnames(fst) <- c("IID","Mean_Fst","Weighted_Fst")
fst <- fst %>% arrange(IID)

fst %>% write_tsv('data/fst/final_fst.tsv')
fst %>% write_tsv('data/phenotypes/final_fst_no_filter.tsv')


