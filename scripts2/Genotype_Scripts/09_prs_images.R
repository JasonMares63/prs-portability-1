library(asbio)
library(cowplot)
library(ggrepel)
library(tidyverse)


load_non_prs_df <- function() {
  # Loads all the covariates, phenotypes, and outcomes into a dataframe
  
  # Load covariates (age, sex, age-sex interactions, PC1, ..., PC20)
  covar_df <- read_tsv('data/ukb_merged/covar_all_samples.covar')
  
  # Load phenotypes file (all phenotypes IRNT'd)
  phenotypes_df <- read_delim('data/phenotypes/full_phenotypes.pheno', delim = ' ', trim_ws = T)
  
  # Load the individuals who are in the evaluation set for each population
  population_files <- c('data/ukb_populations/AFR_all.txt', 'data/ukb_populations/AMR_all.txt',
                        'data/ukb_populations/EAS_all.txt', 'data/ukb_populations/EUR_test.txt',
                        'data/ukb_populations/SAS_all.txt')
  populations_df <- data.frame()
  for (file in population_files) {
    pop_df <- read_delim(file, delim = ' ', trim_ws = T,
                         col_types = c('#FID' = col_integer(), 'IID' = col_integer())) %>%
      mutate(population = str_extract(file, '(?<=data/ukb_populations/)[A-Z]{3}'))
    populations_df <- bind_rows(populations_df, pop_df)
  }
  
  # Combine all the above tables into a table that will be joined with PRS information for
  # phenotype-threshold combinations.
  phenotypes_df %>%
    pivot_longer(Basophil:WBC, names_to = 'phenotype', values_to = 'phenotype_value') %>%
    inner_join(covar_df, by = c('#FID', 'IID')) %>%
    inner_join(populations_df, by = c('#FID', 'IID'))
}


get_r2_values <- function(df) {
  # Computes the partial R^2 attributable to the PRS
  df %>%
    group_by(population, phenotype, threshold) %>%
    do(
      # Regression on only covariates
      nested = lm(phenotype_value ~ sex_covar + age + age_sq + age_sex + age_sq_sex +
                    PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG +
                    PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG + PC11_AVG + PC12_AVG +
                    PC13_AVG + PC14_AVG + PC15_AVG + PC16_AVG + PC17_AVG + PC18_AVG +
                    PC19_AVG + PC20_AVG, data = .),
      
      # Regression now including the PRS + covariates
      full = lm(phenotype_value ~ prs + sex_covar + age + age_sq + age_sex + age_sq_sex +
                  PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG +
                  PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG + PC11_AVG + PC12_AVG +
                  PC13_AVG + PC14_AVG + PC15_AVG + PC16_AVG + PC17_AVG + PC18_AVG +
                  PC19_AVG + PC20_AVG, data = .)
    ) %>%
    mutate(
      partial = partial.R2(nested, full),
      nested_r2 = summary(nested)$adj.r.squared,
      full_r2 = summary(full)$adj.r.squared
    ) %>%
    select(-nested, -full)
}

make_prs_evaluation_df <- function(non_prs_df) {
  # Combine partial R^2 values for each phenotype and threshold
  col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  prs_df <- data.frame()
  for (file in list.files(path = 'data/prs', pattern = '*.sscore', full.names = T)) {
    print(file)
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', prs = 'SCORE1_AVG') %>%
      mutate(
        phenotype = str_extract(string = file, pattern = '(?<=data/prs/)[A-Za-z]+(?=_)'),
        threshold = str_extract(string = file, pattern = '[0-4]') %>% as.integer
      ) %>%
      inner_join(non_prs_df, by =  c('#FID', 'IID','phenotype'))%>%
      get_r2_values
    
    prs_df <- bind_rows(prs_df, this_df)
  }
  prs_df
}


non_prs_df <- load_non_prs_df()

prs_df <- make_prs_evaluation_df(non_prs_df)
prs_df %>% write_tsv('data/prs/combined_evaluation1.tsv')

prs_df <- read_tsv('data/prs/combined_evaluation1.tsv')

colors <- c('#000000', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3')

prs_df1 <- prs_df %>%
  group_by(population,phenotype) %>%
  top_n(1, full_r2) %>%
  ungroup() %>%
  arrange(phenotype) %>%
  mutate(
    population = population %>% factor(levels = c('EUR', 'AMR', 'SAS', 'EAS', 'AFR'))
  )

bar1 <- prs_df1 %>%
  ggplot(aes(x=phenotype,y=full_r2,fill=population))+
  geom_bar(position="dodge",stat="identity") +
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors) 
ggsave("img/maximum_r2_phenotypes_pop.png",bar1,width=8,height=4,dpi=300)



bar2 <- prs_df1 %>%
  ggplot(aes(x=phenotype,y=partial,fill=population))+
  geom_bar(position="dodge",stat="identity") +
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors) 
ggsave("img/maximum_incrementalr2_phenotypes_pop.png",bar2,width=8,height=4,dpi=300)


plot <- prs_df1 %>%
  ggplot(aes(x=nested_r2,y=partial,color=population))+
  geom_point() +
  geom_smooth(method="lm",se=F)+
  geom_text(aes(label=phenotype),size=2)+
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors)

ggsave("img/plot_nested_partial_r2_pop.png",plot,width=8,height=4,dpi=300)

prs_df2 <- prs_df %>%
  group_by(population,phenotype) %>%
  top_n(1, nested_r2) %>%
  ungroup() %>%
  arrange(phenotype) %>%
  mutate(
    population = population %>% factor(levels = c('EUR', 'AMR', 'SAS', 'EAS', 'AFR'))
  )

bar3 <- prs_df2 %>%
  ggplot(aes(x=phenotype,y=nested_r2,fill=population))+
  geom_bar(position="dodge",stat="identity") +
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors) 
ggsave("img/nested_r2_phenotypes_pop.png",bar3,width=8,height=4,dpi=300)


