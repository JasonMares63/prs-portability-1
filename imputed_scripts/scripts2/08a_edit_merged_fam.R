library(tidyverse)

#library(argyle)
#packageVersion("argyle")


fam_file <- read_delim('data/ukb_merged/merged.fam',
                      col_names = c('#FID', 'IID', 'V1', 'V2', 'V3', 'V4'),
                       delim = ' ', trim_ws = T)


#fam_file %>% nrow %>% print

AFR_df <- read_delim('data/ukb_populations/AFR_all.txt', delim = ' ', trim_ws = T)
AMR_df <- read_delim('data/ukb_populations/AMR_all.txt', delim = ' ', trim_ws = T)
EAS_df <- read_delim('data/ukb_populations/EAS_all.txt', delim = ' ', trim_ws = T)
SAS_df <- read_delim('data/ukb_populations/SAS_all.txt', delim = ' ', trim_ws = T)


eur_all <- read_delim('data/ukb_populations/EUR_all.txt', delim = ' ', trim_ws = T)
eur_test <- read_delim('data/ukb_populations/EUR_test.txt', delim = ' ', trim_ws = T)
eur_train <- eur_all %>%
  anti_join(eur_test, by = c('#FID', 'IID'))
eur_train %>% nrow %>% print

combined_labels <- bind_rows(
  eur_train %>% mutate(pop = 'train'),
  eur_test %>% mutate(pop = str_glue('test{IID}') %>% as.character),
  AFR_df %>% mutate(pop = str_glue('test{IID}') %>% as.character),
  AMR_df %>% mutate(pop = str_glue('test{IID}') %>% as.character),
  EAS_df %>% mutate(pop = str_glue('test{IID}') %>% as.character),
  SAS_df %>% mutate(pop = str_glue('test{IID}') %>% as.character),
)

combined <- fam_file %>%
  left_join(combined_labels, by = c('#FID', 'IID'))

# Print a couple QC checks
#combined %>%
#  summarize_all(function(col) sum(is.na(col))) %>%
#  print


#combined %>%
#  group_by(pop) %>%
#  tally %>%
#  arrange(desc(n)) %>%
#  head %>%
#  print

combined %>%
  select(-'#FID') %>%
  select('#FID' = pop, 'IID', starts_with('V')) %>%
  write_tsv('data/ukb_merged/merged_fst.fam')
#  write_tsv('data/ukb_bbj_pgen/merged_fst.fam')

# Compute Fst
# pop_labels <- combined %>% pull(pop)

# Pointer to the genotypes (doesn't load into memory)
# plink_genotypes <- plinkify(c('data/ukb_merged/merged'), where = 'data/fst/')
# plink_genotypes %>% print

# fst_matrix <- argyle::weir.fst.plink(plink_genotypes, flags = '--memory 30000')
# fst_matrix %>% write_tsv('data/ukb_merged/fst.tsv')
