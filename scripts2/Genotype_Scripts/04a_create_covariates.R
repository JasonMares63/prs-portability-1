
library(tidyverse)

#Import Covariates

sex_df = read_table('/rigel/mfplab/projects/ukb_hakhamanesh/phenotype_data/extracted_phenotypes/ukb.sex.txt',
col_names=c("IID","is_male"))
age_df = read_table('/rigel/mfplab/projects/ukb_hakhamanesh/phenotype_data/extracted_phenotypes/ukb.age.txt',
col_names=c('IID', 'age'))
pc_df = read_tsv('data/ukb_merged/projection.sscore')

df <- sex_df %>% full_join(age_df,by="IID") %>% full_join(pc_df,by="IID")
df <- df %>%
	filter(!is.na(age)) %>%
# Code sex as 0 = missing, 1 = male, 2 = female, as in plink .sample files
	mutate(sex_covar=(2-is_male)) %>%
	mutate(sex_covar=replace(sex_covar,is.na(sex_covar),0)) %>%
	mutate(age_sq = age^2) %>%
	mutate(age_sex = age*sex_covar) %>%
	mutate(age_sq_sex = age_sq * sex_covar) %>%
	select(`#FID`,IID,sex_covar,age,age_sq,age_sex,age_sq_sex,contains("PC"))


# Import population files

AFR_df <- read_delim('data/ukb_populations/AFR_all.txt', delim = ' ', trim_ws = T)
AMR_df <- read_delim('data/ukb_populations/AMR_all.txt', delim = ' ', trim_ws = T)
EAS_df <- read_delim('data/ukb_populations/EAS_all.txt', delim = ' ', trim_ws = T)
SAS_df <- read_delim('data/ukb_populations/SAS_all.txt', delim = ' ', trim_ws = T)


eur_all <- read_delim('data/ukb_populations/EUR_all.txt', delim = ' ', trim_ws = T)
eur_test <- read_delim('data/ukb_populations/EUR_test.txt', delim = ' ', trim_ws = T)
eur_train <- eur_all %>%
  anti_join(eur_test, by = c('#FID', 'IID'))

combined_labels <- bind_rows(
  eur_train %>% mutate(pop = 'EUR_train'),
  eur_test %>% mutate(pop = 'EUR_test'),
  AFR_df %>% mutate(pop = 'AFR'),
  AMR_df %>% mutate(pop = 'AMR'),
  EAS_df %>% mutate(pop = 'EAS'),
  SAS_df %>% mutate(pop = 'SAS'),
)

df_pop <- df %>%
  left_join(combined_labels, by = c('#FID', 'IID'))


# Standard normalize with respect to training EUR population

covar = df_pop %>% select(-`#FID`,-IID,-pop) %>% colnames()

train_eur_df <- df_pop %>% filter(pop=="EUR_train")

for (i in covar){
eur_mean <- train_eur_df %>% pull(i) %>% mean(na.rm=T)
eur_sd <- train_eur_df %>% pull(i) %>% sd(na.rm=T)

df_pop[,i] <- (df_pop[,i]-eur_mean)/eur_sd

}

# Export file

df_pop %>% write_tsv('data/ukb_merged/covar_all_samples.covar')

