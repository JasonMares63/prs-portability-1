library(tidyverse)

# Table of each trait's name, code, field ID, and dtype (including eid)
traits_info <- read_delim('/rigel/mfplab/projects/prs-portability/data/martin_gwas_info.txt', delim = ' ',
                          col_types = cols_only(Trait = col_character(),
                                                UKBBcode = col_character())) %>%
    mutate(
        ukb_field = str_glue('{UKBBcode}-0.0') %>% as.character,
        dtype = 'd'
    ) %>%
    add_row(Trait = 'eid', UKBBcode = 'eid', ukb_field = 'eid', dtype = 'c')

# Load only the 18 columns of interest
field_to_dtype <- traits_info %>%
    select(ukb_field, dtype) %>%
    deframe %>%
    as.list

raw_phenotypes_df <- read_csv('/rigel/mfplab/projects/ukb_martin/data/ukb40732.csv',
                              col_types = do.call(cols_only, field_to_dtype))

# Samples table (because not assured that FID == IID)
psam_df <- read.table(
  'data/ukb_merged/merged.fam'
) %>% mutate(V1=as.character(V1)) %>% mutate(V2=as.character(V2))
colnames(psam_df)[1:2] <- c("#FID","IID")

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
combined_labels$IID <- as.character(combined_labels$IID)
combined_labels$`#FID` <- as.character(combined_labels$`#FID`)

pheno_df <- raw_phenotypes_df %>%
    # Rename fields to their human-readable names
    pivot_longer(-eid, names_to = 'ukb_field') %>%
    inner_join(traits_info %>% select(-dtype, -UKBBcode), by = 'ukb_field') %>%
    pivot_wider(id_cols = eid, names_from = Trait, values_from = value) %>%
    rename(IID = eid) %>%
    left_join(psam_df, by = 'IID')  %>%
    left_join(combined_labels,by=c('IID','#FID'))

traits <- c("BMI","WBC","Height","RBC","MCV","MCH","Lymphocyte","Platelet","Monocyte","Eosinophil")

for (i in traits){
    eur_mean <- pheno_df[,c(i,"pop")] %>% filter(pop=="EUR_train") %>% pull(i) %>% mean(na.rm=T)
    eur_sd <- pheno_df[,c(i,"pop")] %>% filter(pop=="EUR_train") %>% pull(i) %>% sd(na.rm=T)

    pheno_df[,i] <- (pheno_df[,i]-eur_mean)/eur_sd
}

#for (i in traits) {
#plot_df <-pheno_df[,c(i,"pop")] %>% filter(pop!="EUR_train")
#
#plot <- plot_df %>% ggplot() +
#	geom_histogram(aes(plot_df[[1]]),color="black", fill="white",bins=15) + 
#    	facet_grid(cols = vars(pop))+
#	theme_light() + 
#    	xlab(i) 
#
#plot_name <- paste0("img/trait_distributions/",i,"_hist_std_norm.png")
#ggsave(plot_name,plot,,width=10,height=8)
#}

pheno_df %>% select('#FID', IID, all_of(traits)) %>%
	write_delim('data/phenotypes/full_phenotypes.pheno', delim = ' ', na = 'NA')
