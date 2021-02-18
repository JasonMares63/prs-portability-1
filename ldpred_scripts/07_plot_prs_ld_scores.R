.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
library(asbio)
library(cowplot)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(bigsnpr)
library(plyr)

load_non_prs_df <- function(num_groups) {
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
  
  fst_values <- read_tsv('data/fst/final_fst.tsv')
  
  # Combine all the above tables into a table that will be joined with PRS information for
  # phenotype-threshold combinations.
  phenotypes_df <- phenotypes_df %>%
    pivot_longer(Basophil:WBC, names_to = 'phenotype', values_to = 'phenotype_value') %>%
    inner_join(covar_df, by = c('#FID', 'IID')) %>%
    inner_join(populations_df, by = c('#FID', 'IID')) %>%
    inner_join(fst_values,by=c("IID")) 
  
  n <- num_groups #number of groups
  cutpoints <<- seq(min(phenotypes_df$Weighted_Fst,na.rm=T),
                    max(phenotypes_df$Weighted_Fst,na.rm=T),length=n+1)
  cutpoints1 <<- seq(min(phenotypes_df$PC1_AVG, na.rm=T),
                     max(phenotypes_df$PC1_AVG,na.rm=T),length=n+1)
  cutpoints2 <<- seq(min(phenotypes_df$PC2_AVG,na.rm=T),
                     max(phenotypes_df$PC2_AVG,na.rm=T),length=n+1)
  
  phenotypes_df <- phenotypes_df %>%
    mutate(PC1_groups = PC1_AVG %>% cut(.,breaks=cutpoints1)) %>% #Divides into equally-sized components
    mutate(PC2_groups = PC2_AVG %>% cut(.,breaks=cutpoints2)) %>%
    mutate(mean_fst_groups = Mean_Fst %>% ntile(9)) %>% #Divides into equally-sized components
    mutate(weighted_fst_groups = Weighted_Fst %>% cut(.,breaks=cutpoints)) %>%  #Divides into equal interval components
    add_count(PC1_groups) %>%
    rename(PC1_groups_count=n) %>%
    add_count(PC2_groups) %>%
    rename(PC2_groups_count=n) %>%
    add_count(weighted_fst_groups) %>%
    rename(fst_groups_count=n)
  
  return(phenotypes_df)
}

get_r2_values <- function(df,group_var) {
  # Computes the partial R^2 attributable to the PRS
  #    group_by(population, phenotype, threshold) %>%
  df$group <- df %>% select(any_of(group_var)) %>% as.data.frame()
  df1 <- df %>%
    group_by(group, phenotype, threshold) %>%
    do(
      # Regression on only covariates
      nested = lm(phenotype_value ~ sex_covar + age + age_sq + age_sex + age_sq_sex +
                    PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG +
                    PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG 
                  #+ PC11_AVG + PC12_AVG +
                  #  PC13_AVG + PC14_AVG + PC15_AVG + PC16_AVG + PC17_AVG + PC18_AVG +
                  #  PC19_AVG + PC20_AVG
                  , data = .),
      
      # Regression now including the PRS + covariates
      full = lm(phenotype_value ~ prs + sex_covar + age + age_sq + age_sex + age_sq_sex +
                  PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG +
                  PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG 
                #+ PC11_AVG + PC12_AVG +
                # PC13_AVG + PC14_AVG + PC15_AVG + PC16_AVG + PC17_AVG + PC18_AVG +
                #PC19_AVG + PC20_AVG
                , data = .)
      #      full_fst = lm(phenotype_value ~ prs + Weighted_Fst + sex_covar + age + age_sq + age_sex + age_sq_sex +
      #                  PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG +
      #                  PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG + PC11_AVG + PC12_AVG +
      #                  PC13_AVG + PC14_AVG + PC15_AVG + PC16_AVG + PC17_AVG + PC18_AVG +
      #                  PC19_AVG + PC20_AVG, data = .)
      
    ) %>%
    mutate(
      partial = partial.R2(nested, full),
      nested_r2 = summary(nested)$adj.r.squared,
      full_r2 = summary(full)$adj.r.squared,
      #      full_fst_r2 = summary(full_fst)$adj.r.squared,
      incremental_r2 = full_r2 - nested_r2
      #      incremental_r2_fst = full_fst_r2 - full_r2
    ) %>%
    select(-nested, -full)
  df1$group <- 1:nrow(df1)
  colnames(df1)[1] <- group_var
  return(df1)
}

make_prs_evaluation_df <- function(non_prs_df,group_var_as_string) {
  # Combine partial R^2 values for each phenotype and threshold
  col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  prs_df <- data.frame()
  for (file in list.files(path = 'data/LDpred/prs',
                          pattern = '[a-zA-Z]+_LDpred_scores_[0-9]{1}[0-9]?.sscore', full.names = T)) {
    print(file)
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', prs = 'SCORE1_AVG') %>%
      mutate(
        phenotype = str_extract(string = file, pattern = '(?<=data/LDpred/prs/)[A-Za-z]+(?=_)'),
        threshold = str_extract(string = file, pattern = '[0-9]+') %>% as.integer
      ) %>%
      inner_join(non_prs_df, by =  c('#FID', 'IID','phenotype'))%>%
      get_r2_values(.,group_var=group_var_as_string)
    
    prs_df <- bind_rows(prs_df, this_df)
  }
  prs_df
}

#Line plots
fst_prs <- function(prs_df){
  
  prs_df$threshold <- (prs_df$threshold-6) %>% as.factor()
  #levels(prs_df$threshold) <- c("5e-8","1e-6","1e-4","1e-3","1e-2")
  levels(prs_df$threshold) <- c(1,round(seq_log(1e-4,0.9,20),4),"Average LDpred p")
  
  thresholds <- prs_df %>%
    subset(phenotype!="Basophil") %>%
    subset(threshold!=0) %>%
    subset(group_number==1) %>%
    group_by(phenotype) %>%
    slice(which.max(partial)) %>%
    select(phenotype,threshold)
  
  plot_df <- data.frame()
  for(i in 1:length(unique(thresholds$phenotype))){
    sub_df <- prs_df %>%
      subset(phenotype!="Basophil") %>%
      subset(threshold !=5) %>%
      filter(phenotype==thresholds$phenotype[i],
             threshold==thresholds$threshold[i])
    plot_df <- rbind(plot_df,sub_df)
    
  }
  
  plot_df %>%
    subset(group_number != max(plot_df$group_number)) %>%
    ggplot(aes(x=group_number,y=partial, 
               group=phenotype,color=phenotype))+
    geom_point(aes(shape=threshold, size=3)) +
    geom_line() + 
    scale_shape_manual(values=1:nlevels(plot_df$threshold)) +
    guides(shape = guide_legend(title="Threshold",
                                labels=levels(plot_df$threshold)),
           color = guide_legend(title="Phenotypes"),
           size=F)+
    scale_x_continuous(breaks=1:(max(plot_df$group_number)),
                       labels=legend_items_fst[1:(max(plot_df$group_number))])+
    theme_classic() + 
    xlab('Fst Groups') +
    ylab(expression(Partial~R^2))+
    labs(title="Decay of Prediction Accuracy Across Fst Pools",
         caption=range_fst)+
    theme(legend.position=c(.85,.65))
}


pcs_prs <- function(prs_df,pc_num){
  prs_df$threshold <- (prs_df$threshold-6) %>% as.factor()
  #levels(prs_df$threshold) <- c("5e-8","1e-6","1e-4","1e-3","1e-2")
  levels(prs_df$threshold) <- c(1,round(seq_log(1e-4,0.9,20),4),"Average LDpred p")
  
  if(pc_num==1){
    legend_items = legend_items_pc1
    range = range_pc1
    eur = max(prs_df$group_number)
  }
  if(pc_num==2){
    legend_items = legend_items_pc2
    range = range_pc2
    eur = min(prs_df$group_number)
  }
  
  thresholds <- prs_df %>%
    subset(phenotype!="Basophil") %>%
    subset(threshold!=30) %>%
    subset(group_number==eur) %>%
    group_by(phenotype) %>%
    slice(which.max(partial)) %>%
    select(phenotype,threshold)
  
  plot_df <- data.frame()
  for(i in 1:length(unique(thresholds$phenotype))){
    sub_df <- prs_df %>%
      subset(phenotype!="Basophil") %>%
      subset(threshold !=30) %>%
      filter(phenotype==thresholds$phenotype[i],
             threshold==thresholds$threshold[i])
    plot_df <- rbind(plot_df,sub_df)
    
  }
  
  
  
  plot_df %>%
    ggplot(aes(x=group_number,y=partial, 
               group=phenotype,color=phenotype))+
    geom_point(aes(shape=threshold, size=3)) +
    geom_line() + 
    guides(shape = guide_legend(title="Threshold",
                                labels=levels(plot_df$threshold)),
           color = guide_legend(title="Phenotypes"),
           size=F)+
    scale_x_continuous(breaks=1:(max(plot_df$group_number)),
                       labels=legend_items[1:(max(plot_df$group_number))])+
    theme_classic() + 
    xlab(paste0('PC',pc_num,' Groups')) +
    ylab(expression(Partial~R^2))+
    labs(title=paste0("Decay of Prediction Accuracy Across PC",pc_num," Pools"),
         caption=range)  
  
}


# Comparison

comparison <- function(prs_df,ldpred_prs_df){
  
#  prs_df$threshold <- (prs_df$threshold-6) %>% as.factor()
  #levels(prs_df$threshold) <- c("5e-8","1e-6","1e-4","1e-3","1e-2")
#  levels(prs_df$threshold) <- c(1,round(seq_log(1e-4,0.9,20),4),"Average LDpred p")
  
  thresholds <- prs_df %>%
    subset(weighted_fst_groups==1) %>%
    group_by(phenotype) %>%
    slice(which.max(partial)) %>%
    select(phenotype,threshold)
  
  plot_df <- data.frame()
  for(i in 1:length(unique(thresholds$phenotype))){
    sub_df <- prs_df %>%
      subset(phenotype!="Basophil") %>%
      subset(threshold !=5) %>%
      filter(phenotype==thresholds$phenotype[i],
             threshold==thresholds$threshold[i])
    plot_df <- rbind(plot_df,sub_df)
    
  }
  
  plot_df2 <- ldpred_prs_df
  
  new_df <- plot_df %>% inner_join(plot_df2, by=c("weighted_fst_groups",
                                                       "phenotype")) %>%
    mutate(diff_partial=(partial.y/partial.x))
  new_df <- new_df %>% inner_join(median_p,by="phenotype")
  
  new_df %>%
    ggplot(aes(x = weighted_fst_groups %>% as.factor(), y = diff_partial, label = phenotype)) +
    geom_violin(aes(fill = weighted_fst_groups %>% as.factor()),
                scale = "count", width = 1.3) +
    geom_point(color="black") +
    geom_label_repel(data = new_df %>% filter(group_number != 10),
                     nudge_x = -0.5, nudge_y = 0.01, size=2) +
    theme_classic() + 
    xlab('Fst Groups') +
    ylab(expression(Ratio~of~LDpred~and~C+T~Partial~R^2))+
    scale_y_log10() +
    labs(title="Improvement of Prediction Accuracy Across Fst Pools",
         caption=range_fst) +
    scale_fill_discrete(name="Fst Intervals (n)",
                        labels=gsub("\\n","",legend_items_fst)) +
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6),
      axis.text = element_text(size = 13, color = 'black'),
      axis.title = element_text(size = 13, color = 'black'))
  
}

####################

non_prs_df <- load_non_prs_df(9)
#prs_df_pop <- make_prs_evaluation_df(non_prs_df,"population")
prs_df_fst <- make_prs_evaluation_df(non_prs_df,"weighted_fst_groups")
prs_df_PC1 <- make_prs_evaluation_df(non_prs_df,"PC1_groups")
prs_df_PC2 <- make_prs_evaluation_df(non_prs_df,"PC2_groups")
                                     
                                     
fst_df <- prs_df_fst %>% 
  arrange(phenotype,weighted_fst_groups,threshold) %>%
  mutate(group_number= weighted_fst_groups) %>%
  filter(weighted_fst_groups!=10)

pc1_df <- prs_df_PC1 %>% 
  arrange(phenotype,PC1_groups,threshold) %>%
  mutate(group_number= PC1_groups) %>%
  filter(PC1_groups!=10)

pc2_df <- prs_df_PC2 %>% 
  arrange(phenotype,PC2_groups,threshold) %>%
  mutate(group_number= PC2_groups) %>%
  filter(PC2_groups!=10)


fst_df %>% write_tsv('data/LDpred/prs/UKBB_geno_fst.tsv')
fst_df <- read_tsv('data/LDpred/prs/UKBB_geno_fst.tsv')
pc1_df %>% write_tsv('data/LDpred/prs/UKBB_geno_pc1.tsv')
pc2_df %>% write_tsv('data/LDpred//UKBB_geno_pc2.tsv')


fst_df_ct <- read_tsv('data/prs/UKBB_genotype_1KG_EUR_LD/UKBB_genotype_scores_fst.tsv')

#######################

fst_values <- non_prs_df %>% 
  filter(phenotype=="BMI") %>%
  select(weighted_fst_groups)
fst_values1 <- fst_values[[1]] %>% levels() %>% as.character() %>%
  strsplit(split=",") %>% unlist() %>% parse_number() %>% round(digits=2)

pcs <- non_prs_df %>% filter(phenotype=="BMI") %>%
  select(PC1_groups,PC2_groups)
PC1_values <- pcs[[1]] %>% levels() %>% as.character() %>%
  strsplit(split=",") %>% unlist() %>% parse_number() %>% round(digits=2)
PC2_values <- pcs[[2]] %>% levels() %>% as.character() %>%
  strsplit(split=",") %>% unlist() %>% parse_number() %>% round(digits=2)

legend_items_fst <- c()
#interval <- unlist(fst_values) %>% levels()
count <- (table(fst_values[[1]]))
for(i in 1:9){
  add <- paste0("(",fst_values1[2*i-1],", ",fst_values1[2*i],"]\n",
                " (",count[i],")")
  legend_items_fst <- c(legend_items_fst,add)
}
range_fst <- paste0("Fst Intervals are sized approximately: ",fst_values1[2]-fst_values1[1])


legend_items_pc1 <- c()
#interval <- unlist(fst_values) %>% levels()
count <- (table(pcs[[1]]))
for(i in 1:9){
  add <- paste0("(",PC1_values[2*i-1],", ",PC1_values[2*i],"]\n",
                " (",count[i],")")
  legend_items_pc1 <- c(legend_items_pc1,add)
}
range_pc1 <- paste0("PC1 Intervals are sized approximately: ",PC1_values[2]-PC1_values[1])


legend_items_pc2 <- c()
#interval <- unlist(fst_values) %>% levels()
count <- (table(pcs[[2]]))
for(i in 1:9){
  add <- paste0("(",PC2_values[2*i-1],", ",PC2_values[2*i],"]\n",
                " (",count[i],")")
  legend_items_pc2 <- c(legend_items_pc2,add)
}
range_pc2 <- paste0("PC2 Intervals are sized approximately: ",PC2_values[2]-PC2_values[1])


df <- data.frame(p = seq_log(1e-4,0.9,20))
for (file in list.files(path = 'data/LDpred',
                        pattern = '[a-zA-Z]+_ldpred_auto_p_est.tsv', full.names = T)) {
  phenotype <- str_extract(string = file, pattern = '(?<=data/LDpred/)[A-Za-z]+(?=_)')
  add <- read_tsv(file,col_names=F,skip=1)
  pheno <- function(pheno){
    return(phenotype)
  }
  add <- add %>% rename_with(pheno)
  df <- bind_cols(df,add)
  
}

median_p <- ldply(df,median) 
colnames(median_p) <- c("phenotype","prop")
median_p <- median_p[-1,]
median_p[,2] <- round(median_p[,2],3)

################################
fst_df_avg <- fst_df %>% filter(threshold==28)
pc1_df_avg <- pc1_df %>% filter(threshold==28)
pc2_df_avg <- pc2_df %>% filter(threshold==28)

fig3_continuous_fst <- fst_prs(fst_df)
ggsave('img/LDpred_UKBB_geno_Fst_lines.png', fig3_continuous_fst, width = 6, height = 10, dpi = 300)


fig3_continuous_fst <- fst_prs(fst_df_avg)
ggsave('img/LDpred_avg_UKBB_geno_Fst_lines.png', fig3_continuous_fst, width = 6, height = 10, dpi = 300)

fig3_continuous_pc1_lines <- pcs_prs(pc1_df_avg,1)
ggsave('img/LDpred_avg_UKBB_geno_pc1_lines.png', fig3_continuous_pc1_lines, width = 6*1.4, height = 5*1.25, dpi = 300)


fig3_continuous_pc2_lines <- pcs_prs(pc2_df_avg,2)
ggsave('img/LDpred_avg_UKBB_geno_pc2_lines.png', fig3_continuous_pc2_lines, width = 6*1.4, height = 5*1.25, dpi = 300)


comparison_fst <- comparison(fst_df_ct,fst_df_avg)
ggsave('img/LDpred_avg_CT_UKBB_geno_violin1.png', comparison_fst, width = 6*1.4, height = 5*1.25, dpi = 300)
