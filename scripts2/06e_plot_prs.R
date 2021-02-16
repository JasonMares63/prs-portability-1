
library(asbio)
library(cowplot)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)

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

  fst_values <- read_tsv('data/fst/final_fst.tsv')

    # Combine all the above tables into a table that will be joined with PRS information for
  # phenotype-threshold combinations.
  phenotypes_df <- phenotypes_df %>%
    pivot_longer(Basophil:WBC, names_to = 'phenotype', values_to = 'phenotype_value') %>%
    inner_join(covar_df, by = c('#FID', 'IID')) %>%
    inner_join(populations_df, by = c('#FID', 'IID')) %>%
    inner_join(fst_values,by=c("IID")) 
  
  n <- 9 #number of groups
  cutpoints <<- seq(min(phenotypes_df$Weighted_Fst,na.rm=T),
                    max(phenotypes_df$Weighted_Fst,na.rm=T),length=9+1)
  cutpoints1 <<- seq(min(phenotypes_df$PC1_AVG, na.rm=T),
                     max(phenotypes_df$PC1_AVG,na.rm=T),length=9+1)
  cutpoints2 <<- seq(min(phenotypes_df$PC2_AVG,na.rm=T),
                     max(phenotypes_df$PC2_AVG,na.rm=T),length=9+1)
  
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

#file <- "data/prs/Martin_et_al_approach_scores/WBC_0_scores.sscore"
#group_var_as_string <- "weighted_fst_groups"
#col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
 #             'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())


make_prs_evaluation_df <- function(non_prs_df,group_var_as_string) {
  # Combine partial R^2 values for each phenotype and threshold
  col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  prs_df <- data.frame()
  for (file in list.files(path = 'data/prs',
                          pattern = '[a-zA-Z]+_[0-9]_scores.sscore', full.names = T)) {
    print(file)
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', prs = 'SCORE1_AVG') %>%
      mutate(
        phenotype = str_extract(string = file, pattern = '(?<=data/prs/)[A-Za-z]+(?=_)'),
        threshold = str_extract(string = file, pattern = '[0-4]') %>% as.integer
      ) %>%
      inner_join(non_prs_df, by =  c('#FID', 'IID','phenotype'))%>%
      get_r2_values(.,group_var=group_var_as_string)
    
    prs_df <- bind_rows(prs_df, this_df)
  }
  prs_df
}

#Violin Plots
plot_figure_3_incr <- function(prs_df) {
  # Plot the results as per Martin et al. Figure 3
  colors <- c('#000000', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3')
  plot_df <- prs_df %>%
    subset(phenotype!="Basophil") %>%
    subset(threshold!=5) %>%
    group_by(population, phenotype) %>%
    top_n(1, partial) %>%
    ungroup() %>%
    group_by(phenotype) %>%
    mutate(
      eur_partial_incr = max(incremental_r2 * (population == 'EUR')),
      relative_performance_incr = incremental_r2 / eur_partial_incr,
      population = population %>% factor(levels = c('EUR', 'AMR', 'SAS', 'EAS', 'AFR'))
    )
  
  se_df <- plot_df %>%
    ungroup %>%
    group_by(population) %>%
    mutate(mean_rel_eur_incr = mean(relative_performance_incr),
           sem_rel_eur_incr = sd(relative_performance_incr)/sqrt(length(unique(plot_df$phenotype)))) %>%
    subset(phenotype=='Height') %>%
    arrange(desc(mean_rel_eur_incr))
  
    plot_df %>%
      ggplot(aes(x = population, y = relative_performance_incr, label = phenotype)) +
      geom_violin(aes(fill = population, color = population), scale = "count", width = 1.3) +
      geom_point(color = 'black') +
      geom_errorbar(data=se_df, aes(ymin=mean_rel_eur_incr - sem_rel_eur_incr,
                                    ymax=mean_rel_eur_incr + sem_rel_eur_incr),
                    color='black', width=0.1) +
      geom_label_repel(data = plot_df %>% filter(population != 'EUR'),
                       nudge_x = -0.5, nudge_y = 0.01, size=2) +
      scale_colour_manual(values = colors) +
      scale_fill_manual(values = colors) +
      theme_classic() +
      theme(
        legend.position = 'NONE',
        axis.text = element_text(size = 13, color = 'black'),
        axis.title = element_text(size = 13, color = 'black')
      ) +
      xlab('Super population') +
      ylab(expression(Incremental~R^2~relative~to~Europeans))
  
}

plot_figure_3_partial <- function(prs_df) {
  # Plot the results as per Martin et al. Figure 3
  colors <- c('#000000', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3')
  plot_df <- prs_df %>%
    filter(phenotype != "Basophil") %>%
    filter(threshold!=5) %>%
    group_by(population, phenotype) %>%
    top_n(1, partial) %>%
    ungroup() %>%
    group_by(phenotype) %>%
    mutate(
      eur_partial = max(partial * (population == 'EUR')),
      relative_performance = partial / eur_partial,
      eur_partial_incr = max(incremental_r2 * (population == 'EUR')),
      relative_performance_incr = incremental_r2 / eur_partial_incr,
      population = population %>% factor(levels = c('EUR', 'AMR', 'SAS', 'EAS', 'AFR'))
    )
  
  se_df <- plot_df %>%
    ungroup %>%
    group_by(population) %>%
    dplyr::mutate(mean_rel_eur=mean(relative_performance),
                  sem_rel_eur=sd(relative_performance)/length(unique(plot_df$phenotype))) %>%
    subset(phenotype=='Height') %>%
    arrange(desc(mean_rel_eur))
  
plot_df %>%
  ggplot(aes(x = population, y = relative_performance, label = phenotype)) +
  geom_violin(aes(fill = population, color = population), scale = "count", width = 1.3) +
  geom_point(color = 'black') +
  geom_errorbar(data=se_df, aes(ymin=mean_rel_eur - sem_rel_eur, ymax=mean_rel_eur + sem_rel_eur),
                color='black', width=0.1) +
  geom_label_repel(data = plot_df %>% filter(population != 'EUR'),
                   nudge_x = -0.5, nudge_y = 0.01, size=2) +
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  theme(
    legend.position = 'NONE',
    axis.text = element_text(size = 13, color = 'black'),
    axis.title = element_text(size = 13, color = 'black')
  ) +
  xlab('Super population') +
  ylab(expression(Partial~R^2~relative~to~Europeans))
}

plot_figure_3_partial_fst <- function(prs_df) {
  # Plot the results as per Martin et al. Figure 3

  plot_df <- prs_df %>%
    filter(phenotype != "Basophil") %>%
    filter(threshold!=5) %>%
    group_by(group_number, phenotype) %>%
    top_n(1, partial) %>%
    ungroup() %>%
    group_by(phenotype) %>%
    mutate(
      eur_partial = max(partial * (group_number == 1)),
      relative_performance = partial / eur_partial,
      group_number = group_number %>% as.factor()
   )
  
  se_df <- plot_df %>%
    ungroup %>%
    group_by(group_number) %>%
    dplyr::mutate(mean_rel_eur=mean(relative_performance),
                  sem_rel_eur=sd(relative_performance)/length(unique(plot_df$phenotype))) %>%
    subset(phenotype=='Height') %>%
    arrange(desc(mean_rel_eur))
  
  plot_df %>%
    ggplot(aes(x = group_number, y = relative_performance, label = phenotype)) +
    geom_violin(aes(fill = group_number),
                scale = "count", width = 1.3) +
    geom_point(color = 'black') +
    geom_smooth(data=plot_df %>% filter(group_number!=1) %>% mutate(group_number=as.numeric(group_number)),
                method="lm",se=FALSE,color="black")+
    geom_errorbar(data=se_df, aes(ymin=mean_rel_eur - sem_rel_eur, ymax=mean_rel_eur + sem_rel_eur),
                  color='black', width=0.1) +
    geom_label_repel(data = plot_df %>% filter(group_number != 1),
                     nudge_x = -0.5, nudge_y = 0.01, size=2) +
    theme_classic() + 
    xlab('Fst Groups') +
    ylab(expression(Partial~R^2~relative~to~Europeans))+
    labs(title="Decay of Prediction Accuracy Across Fst Pools",
         caption=range_fst) +
    scale_fill_discrete(name="Fst Intervals (n)",
                     labels=gsub("\\n","",legend_items_fst))+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6),
      axis.text = element_text(size = 13, color = 'black'),
      axis.title = element_text(size = 13, color = 'black')
    )
}


#Line plots
fst_prs <- function(prs_df){
  
  prs_df$threshold <- prs_df$threshold %>% as.factor()
  levels(prs_df$threshold) <- c("5e-8","1e-6","1e-4","1e-3","1e-2")
  
 thresholds <- prs_df %>%
   subset(phenotype!="Basophil") %>%
   subset(threshold!=5) %>%
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
    ggplot(aes(x=group_number,y=partial, 
               group=phenotype,color=phenotype))+
    geom_point(aes(shape=threshold, size=3)) +
    geom_line() + 
    guides(shape = guide_legend(title="Threshold",
                                labels=levels(plot_df$threshold)),
           color = guide_legend(title="Phenotypes"),
           size=F)+
    scale_x_continuous(breaks=1:(max(plot_df$group_number)-1),
                     labels=legend_items_fst[1:(max(plot_df$group_number)-1)])+
    theme_classic() + 
    xlab('Fst Groups') +
    ylab(expression(Partial~R^2))+
    labs(title="Decay of Prediction Accuracy Across Fst Pools",
         caption=range_fst)+
    theme(legend.position=c(.75,.65))
}


pcs_prs <- function(prs_df,pc_num){
  prs_df$threshold <- prs_df$threshold %>% as.factor()
  levels(prs_df$threshold) <- c("5e-8","1e-6","1e-4","1e-3","1e-2")
  
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
    subset(threshold!=5) %>%
    subset(group_number==eur) %>%
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


pc_plot <- function(covar_df){
  covar_df %>%
    ggplot(aes(x=PC1_AVG,y=PC2_AVG,color=population))+
    geom_point()+
    labs(title="PCA1 vs PCA2 of testing UKBB individuals")
  
} 
  

pop_prs <- function(prs_df){
  
  plot_df <- prs_df %>%
    subset(phenotype!="Basophil") %>%
    subset(threshold<5) %>%
    group_by(population, phenotype) %>%
    top_n(1, incremental_r2) %>%
    ungroup()  %>%
  #    group_by(phenotype) %>%
      mutate(
        population = population %>% factor(levels = c('EUR', 'AMR', 'SAS', 'EAS', 'AFR'))
     )
  
  
  
  colors <- c('#000000', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3')
  plot_df$threshold <- plot_df$threshold %>% as.factor()
  levels(plot_df$threshold) <- c("5e-8","1e-6","1e-4","1e-3","1e-2")
  
  
  plot_df %>%
    ggplot(aes(x=population,y=partial, 
               group=phenotype,color=phenotype))+
    geom_point(aes(shape=threshold, size=6)) +
    geom_line() + 
    scale_fill_brewer(palette = "Set1")+
    guides(size = FALSE)+
    scale_fill_discrete(name = "threshold", labels = c("5e-8", "1e-6", "1e-4","1e-3","1e-2"))
  
} #Line graphs across populations
####################

non_prs_df <- load_non_prs_df()
prs_df_fst <- make_prs_evaluation_df(non_prs_df,"weighted_fst_groups")
prs_df_PC1 <- make_prs_evaluation_df(non_prs_df,"PC1_groups")
prs_df_PC2 <- make_prs_evaluation_df(non_prs_df,"PC2_groups")

cutpoints1 #PC1
cutpoints2 #PC2
cutpoints #Fst

fst_df <- prs_df_fst %>% 
  mutate(group_number= weighted_fst_groups) %>%
  filter(weighted_fst_groups!=10) 

pc1_df <- prs_df_PC1 %>%
  mutate(group_number= PC1_groups) %>%
  filter(PC1_groups!=10) 

pc2_df <- prs_df_PC2 %>% 
  mutate(group_number= PC2_groups) %>%
  filter(PC2_groups!=10) 


fst_df %>% write_tsv('data/prs/UKBB_genotype_scores_fst.tsv')
pc1_df %>% write_tsv('data/prs/UKBB_genotype_scores_pc1.tsv')
pc2_df %>% write_tsv('data/prs/UKBB_genotype_scores_pc2.tsv')

#fst_df <- read_tsv('data/prs/Martin_et_al_approach_scores_FST_groups.tsv')
#pc1_df <- read_tsv('data/prs/Martin_et_al_approach_scores_PC1_groups.tsv')
#pc2_df <- read_tsv('data/prs/Martin_et_al_approach_scores_PC2_groups.tsv')

fst_df <- fst_df %>% filter(!phenotype %in%
                       c("Basophil","Hb","DBP",
                         "Ht","MCHC", "Neutrophil","SBP"))
pc1_df <- pc1_df %>% filter(!phenotype %in%
                              c("Basophil","Hb","DBP",
                                "Ht","MCHC", "Neutrophil","SBP"))
pc2_df <- pc2_df %>% filter(!phenotype %in%
                              c("Basophil","Hb","DBP",
                                "Ht","MCHC", "Neutrophil","SBP"))
non_prs_df_sub <- non_prs_df %>% filter(!phenotype %in%
                                          c("Basophil","Hb","DBP",
                                            "Ht","MCHC", "Neutrophil","SBP"))
#############################

fst_values <- non_prs_df_sub %>% select(weighted_fst_groups)
fst_values1 <- fst_values[[1]] %>% levels() %>% as.character() %>%
  strsplit(split=",") %>% unlist() %>% parse_number() %>% round(digits=2)

pcs <- non_prs_df_sub %>% select(PC1_groups,PC2_groups)
PC1_values <- pcs[[1]] %>% levels() %>% as.character() %>%
  strsplit(split=",") %>% unlist() %>% parse_number() %>% round(digits=2)
PC2_values <- pcs[[2]] %>% levels() %>% as.character() %>%
  strsplit(split=",") %>% unlist() %>% parse_number() %>% round(digits=2)

legend_items_fst <- c()
#interval <- unlist(fst_values) %>% levels()
count <- (table(fst_values[[1]]))/10
for(i in 1:9){
  add <- paste0("(",fst_values1[2*i-1],", ",fst_values1[2*i],"]\n",
                " (",count[i],")")
  legend_items_fst <- c(legend_items_fst,add)
}
range_fst <- paste0("Fst Intervals are sized approximately: ",fst_values1[2]-fst_values1[1])


legend_items_pc1 <- c()
#interval <- unlist(fst_values) %>% levels()
count <- (table(pcs[[1]]))/10
for(i in 1:9){
  add <- paste0("(",PC1_values[2*i-1],", ",PC1_values[2*i],"]\n",
                " (",count[i],")")
  legend_items_pc1 <- c(legend_items_pc1,add)
}
range_pc1 <- paste0("PC1 Intervals are sized approximately: ",PC1_values[2]-PC1_values[1])


legend_items_pc2 <- c()
#interval <- unlist(fst_values) %>% levels()
count <- (table(pcs[[2]]))/10
for(i in 1:9){
  add <- paste0("(",PC2_values[2*i-1],", ",PC2_values[2*i],"]\n",
                " (",count[i],")")
  legend_items_pc2 <- c(legend_items_pc2,add)
}
range_pc2 <- paste0("PC2 Intervals are sized approximately: ",PC2_values[2]-PC2_values[1])



fig3_continuous_fst <- plot_figure_3_partial_fst(fst_df)
#ggsave('img/fig3_continuous_fst.png', fig3_continuous_fst, width = 6, height = 10, dpi = 300)
ggsave('img/UKBB_geno_continuous_fst.png', fig3_continuous_fst, width = 6, height = 10, dpi = 300)


fig3_continuous_fst_lines <- fst_prs(fst_df)
#ggsave('img/fig3_continuous_fst_lines.png', fig3_continuous_fst_lines, width = 6*1.4, height = 5*1.25, dpi = 300)
ggsave('img/UKBB_geno_continuous_fst_lines.png', fig3_continuous_fst_lines, width = 6*1.4, height = 5*1.25, dpi = 300)

fig3_continuous_pc1_lines <- pcs_prs(pc1_df,1)
#ggsave('img/fig3_continuous_pc1_lines.png', fig3_continuous_pc1_lines, width = 6*1.4, height = 5*1.25, dpi = 300)
ggsave('img/UKBB_geno_continuous_pc1_lines.png', fig3_continuous_pc1_lines, width = 6*1.4, height = 5*1.25, dpi = 300)

fig3_continuous_pc2_lines <- pcs_prs(pc2_df,2)
#ggsave('img/fig3_continuous_pc2_lines.png', fig3_continuous_pc2_lines, width = 6*1.4, height = 5*1.25, dpi = 300)
ggsave('img/UKBB_geno_continuous_pc2_lines.png', fig3_continuous_pc2_lines, width = 6*1.4, height = 5*1.25, dpi = 300)


pca_plot <- pc_plot(non_prs_df)
ggsave('img/test_UKBB_pcas.png',pca_plot,width=6,height=10,dpi=300)

##############
#Boxplot of FSt distributed across population groups
boxplots <- non_prs_df %>% filter(phenotype=="Height") %>%
  ggplot(aes(x=population, y=Weighted_Fst, color=population)) +
  geom_boxplot()+
  geom_hline(yintercept = fst_avg[,2],linetype="dashed",size=0.3)+
  scale_x_discrete(limits=c("EUR", "AMR", "SAS","EAS","AFR"))
ggsave('img/fst_boxplots.png', boxplots, width = 8, height = 4, dpi = 300)


#Line graphs with fst and allowing optimization across each group
fst_prs_partial <- fst_prs(prs_df1)
ggsave('img/fst_prs_partial_phenotypes.png', fst_prs_partial, width = 6, height = 10, dpi = 300)


prs_df1 %>% write_tsv('data/prs/fst_subset_prs_CT.tsv')

###########

#Replace fst_group with population group in get_r2_values function (using ntiles())
prs_df1 <- read_tsv('data/prs/subset_evaluation_jam.tsv')

pop_prs_partial <- pop_prs(prs_df1)
ggsave('img/pop_prs_partial_phenotypes.png', pop_prs_partial, width = 6, height = 10, dpi = 300)


#######################################
# Violin Plots from populations of different methods


prs_df1 <-  read_tsv('data/prs/subset_evaluation_jam.tsv')
fig3_incr <- plot_figure_3_partial(prs_df1)
ggsave('img/pop_partial_prs_violin_CT.png', fig3_incr, width = 8, height = 4, dpi = 300)


prs_df1 <-  read_tsv('data/prs/fst_subset_prs_CT.tsv')
fig3_incr <- plot_figure_3_partial_fst(prs_df)
ggsave('img/fst_partial_prs_violin_CT.png', fig3_incr, width = 8, height = 4, dpi = 300)
