.libPaths("/rigel/mfplab/users/jm4454/rpackages/")
library(asbio)
library(cowplot)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
#library(bigsnpr)
library(Gmedian)
library(philentropy)

load_non_prs_df <- function(num_pc_groups,num_fst_groups) {
  # Loads all the covariates, phenotypes, and outcomes into a dataframe
  
  # Load covariates (age, sex, age-sex interactions, PC1, ..., PC20)
  covar_df <- read_tsv('data/ukb_merged/covar_all_samples.covar')
  medians <- covar_df %>% filter(pop=="EUR_train") %>% 
    select(starts_with("PC")) %>% Gmedian() %>% as.vector()
  covar_df <- covar_df %>% filter(pop!="EUR_train")
  dist_fun <- function(row){
    if (is.na(row)>0){
      return(NA)
    } 
    else{ 
    return(distance(rbind(row,medians),method="euclidean"))
  
      }
  }
  
  dist <- covar_df %>% select(starts_with("PC")) %>% as.matrix() %>%
    apply(1,dist_fun)
  covar_df$PC_dist <- dist
  
  #gg <- covar_df %>% filter(pop!="EUR_train") %>%
  #    ggplot(aes(x=PC_dist)) + geom_histogram(bins=num_pc_groups)
  #ggsave('img/PC_distance_hist.png',gg)
  
  # Load phenotypes file (all phenotypes not IRNT'd)
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
  fst_values <- read_tsv('/rigel/mfplab/users/jm4454/prs-portability-master/data/fst/final_fst.tsv')  
  #fst_values <- read_tsv('data/fst/final_fst.tsv')
  
  # Combine all the above tables into a table that will be joined with PRS information for
  # phenotype-threshold combinations.
  phenotypes_df <- phenotypes_df %>%
    inner_join(covar_df, by = c('#FID', 'IID')) %>%
    inner_join(populations_df, by = c('#FID', 'IID')) %>%
    left_join(fst_values,by=c("IID")) %>% #fix this once we get all Fsts 
    mutate(Weighted_Fst=ifelse(is.na(Weighted_Fst)*(population=="EUR"),
                             0,
                             Weighted_Fst)) %>%
    pivot_longer(BMI:Eosinophil, names_to = 'phenotype', values_to = 'phenotype_value')
  
 
  cutpoints_fst <<- seq(min(phenotypes_df$Weighted_Fst,na.rm=T),
                    max(phenotypes_df$Weighted_Fst,na.rm=T),
                    length=num_fst_groups+1)
  cutpoints_pcs <<- seq(min(phenotypes_df$PC_dist, na.rm=T),
                     max(phenotypes_df$PC_dist,na.rm=T),
                     length=num_pc_groups+1)
  
  phenotypes_df <- phenotypes_df %>%
    mutate(PC_groups = PC_dist %>% cut(.,breaks=cutpoints_pcs)) %>% #Divides into equally-sized components
#    mutate(mean_fst_groups = Mean_Fst %>% ntile(9)) %>% #Divides into equally-sized components
    mutate(weighted_fst_groups = Weighted_Fst %>% cut(.,breaks=cutpoints_fst)) %>%  #Divides into equal interval components
    add_count(PC_groups) %>%
    rename(PC_groups_count=n) %>%
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
                    PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG + PC11_AVG + PC12_AVG +
                    PC13_AVG + PC14_AVG + PC15_AVG + PC16_AVG + PC17_AVG + PC18_AVG +
                    PC19_AVG + PC20_AVG
                  , data = .),
      
      # Regression now including the PRS + covariates
      full = lm(phenotype_value ~ prs + sex_covar + age + age_sq + age_sex + age_sq_sex +
                  PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG +
                  PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG + PC11_AVG + PC12_AVG +
                  PC13_AVG + PC14_AVG + PC15_AVG + PC16_AVG + PC17_AVG + PC18_AVG +
                  PC19_AVG + PC20_AVG,
                  data = .)
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

make_ldpred_evaluation_df <- function(non_prs_df,group_var_as_string) {
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


####################

non_prs_df <- load_non_prs_df(num_pc_groups = 8,num_fst_groups = 9)
#prs_df_pop <- make_prs_evaluation_df(non_prs_df,"population")
ld_df_fst <- make_ldpred_evaluation_df(non_prs_df,"weighted_fst_groups")
ld_df_PC <- make_ldpred_evaluation_df(non_prs_df,"PC_groups")
prs_df_fst  <- make_prs_evaluation_df(non_prs_df,"weighted_fst_groups")
prs_df_PC <- make_prs_evaluation_df(non_prs_df,"PC_groups")
                                     

                                     
prs_fst_df <- prs_df_fst %>% 
  arrange(phenotype,weighted_fst_groups,threshold) %>%
  mutate(group_number = weighted_fst_groups) %>%
  filter(weighted_fst_groups!=10)

prs_pc_df <- prs_df_PC %>% 
  arrange(phenotype,PC_groups,threshold) %>%
  mutate(group_number = PC_groups) %>%
  filter(PC_groups!=9)

ld_fst_df <- ld_df_fst %>% 
  arrange(phenotype,weighted_fst_groups,-threshold) %>%
  mutate(group_number = weighted_fst_groups) %>%
  filter(weighted_fst_groups!=10)

ld_pc_df <- ld_df_PC %>% 
  arrange(phenotype,PC_groups,-threshold) %>%
  mutate(group_number = PC_groups) %>%
  filter(PC_groups!=9)

non_prs_df %>% write_tsv('data/prs_comparisons/non_prs_df.tsv')
prs_fst_df %>% write_tsv('data/prs_comparisons/UKBB_geno_fst_CT.tsv')
ld_fst_df %>% write_tsv('data/prs_comparisons/UKBB_geno_fst_LD.tsv')
prs_pc_df %>% write_tsv('data/prs_comparisons/UKBB_geno_pc_CT.tsv')
ld_pc_df %>% write_tsv('data/prs_comparisons/UKBB_geno_pc_LD.tsv')

non_prs_df <- read_tsv('data/prs_comparisons/non_prs_df.tsv')
prs_fst_df <- read_tsv('data/prs_comparisons/UKBB_geno_fst_CT.tsv')
ld_fst_df <- read_tsv('data/prs_comparisons/UKBB_geno_fst_LD.tsv')
prs_pc_df <- read_tsv('data/prs_comparisons/UKBB_geno_pc_CT.tsv')
ld_pc_df <- read_tsv('data/prs_comparisons/UKBB_geno_pc_LD.tsv')
#######################

fst_values <- non_prs_df %>% 
  filter(phenotype=="BMI") %>%
  select(weighted_fst_groups)
fst_values1 <- fst_values[[1]] %>% levels() %>% as.character() %>%
  strsplit(split=",") %>% unlist() %>% parse_number() %>% round(digits=2)

pcs <- non_prs_df %>% filter(phenotype=="BMI") %>%
  select(PC_groups)
PC_values <- pcs[[1]] %>% levels() %>% as.character() %>%
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


legend_items_pc <- c()
#interval <- unlist(fst_values) %>% levels()
count <- (table(pcs[[1]]))
for(i in 1:9){
  add <- paste0("(",PC_values[2*i-1],", ",PC_values[2*i],"]\n",
                " (",count[i],")")
  legend_items_pc <- c(legend_items_pc,add)
}
range_pc <- paste0("PC Distance Intervals are sized approximately: ",PC_values[2]-PC_values[1])

mean_fst_values <- non_prs_df %>% 
  filter(phenotype=="Height") %>%
  select(Weighted_Fst,weighted_fst_groups) %>%
  group_by(weighted_fst_groups) %>%
  mutate(mean_fst = mean(Weighted_Fst)) %>%
  ungroup() %>%
  mutate(weighted_fst_groups = as.numeric(weighted_fst_groups)) %>%
  distinct(weighted_fst_groups,mean_fst) %>%
  na.omit() %>% arrange(weighted_fst_groups)

mean_pc_values <- non_prs_df %>% 
  filter(phenotype=="Height") %>%
  select(PC_dist,PC_groups) %>%
  group_by(PC_groups) %>%
  mutate(mean_pc = mean(PC_dist)) %>%
  ungroup() %>%
  mutate(PC_groups = as.numeric(PC_groups)) %>%
  distinct(PC_groups,mean_pc) %>%
  na.omit() %>% arrange(PC_groups) %>%
  mutate(mean_pc = as.numeric(mean_pc))

################### Initial sparsity values in LDpred
{
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
}
################################



rel_partial_CT <- function(prs_df,type) {
  # C+ T
  
  prs_df$threshold <- prs_df$threshold %>% as.factor()
  levels(prs_df$threshold) <- c("5e-8","1e-5","1e-4","1e-3","1e-2")
  
  
  plot_df <- prs_df %>%
    filter(threshold!="5e-8") %>%
    group_by(phenotype,threshold) %>%
    mutate(
      eur_partial = max((partial*(group_number==1)),na.rm=T),
      relative_performance = partial / eur_partial) %>%
    ungroup() 
  
    
  if (type=="FST"){
    plot_df <- plot_df %>% left_join(mean_fst_values,by="weighted_fst_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_fst)
    title1 = "Decay of Relative Prediction Accuracy Across Fst Pools"
    xlabel1 = "Group Average Fst"
    subtitle1 = "Each Fst group is represented along the x-axis by the within-group Fst average."
    caption1 = "Total of 8 groups of interval size 0.02, with lowest and highest Fst values at -0.02 and 0.14, respectively.\n
         Note: The fifth group contained less than 1000 observations; all other groups n>2500."
  }
  if (type=="PC") {
    plot_df <- plot_df %>% left_join(mean_pc_values,by="PC_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_pc)
    title1 = "Decay of Relative Prediction Accuracy Across PC Distance Pools"
    xlabel1 = "Group Average PC Distance"
    subtitle1 = "Each PC distance group is represented along the x-axis by the within-group PC distance from European training population average."
    caption1 = range_pc
  }
  
 se_df <- plot_df %>%
    group_by(group_number,threshold) %>%
   mutate(
      mean_rel_eur=mean(relative_performance),
     sem_rel_eur=sd(relative_performance)/sqrt(10)) %>%
    subset(phenotype=='Height') %>%
    arrange(desc(mean_rel_eur))
  
  
  plot_df <- plot_df %>%
    filter(group_number != max(plot_df$group_number))
  
  
    ggplot(plot_df, aes(x=mean_values,y=as.numeric(relative_performance), 
               color=phenotype))+
    geom_point(size=3,alpha=0.3)+
    geom_smooth(data=plot_df %>% filter(group_number!=1),
               method="lm",se=T,color="black",na.rm=T)+
    guides(color = guide_legend(title="Phenotypes"),
           size=F)+
    theme_light() + 
    scale_fill_brewer(palette="Set3")+
    xlab(xlabel1) +
    ylab(expression(Partial~R^2~relative~to~Europeans))+
    labs(title=title1,
         subtitle=subtitle1,
         caption=caption1)+
    #theme(legend.position=c(.90,.80))+
    ylim(0, 1)+
    facet_grid(cols = vars(threshold))
  
}

fst_CT_graph <- rel_partial_CT(prs_fst_df,"FST")
ggsave('img/CT_FST.png', fst_CT_graph, width = 20, height = 10, dpi = 300)


pc_CT_graph <- rel_partial_CT(prs_pc_df,"PC")
ggsave('img/CT_PC.png', pc_CT_graph, width = 20, height = 10, dpi = 300)



rel_partial_LDpred <- function(prs_df, type) {
  
  plot_df <- prs_df %>%
    filter(threshold==28) %>%
    group_by(phenotype) %>%
    mutate(
      eur_partial = max(partial*(group_number==1),na.rm=T),
      relative_performance = partial / eur_partial) %>%
    ungroup() 
  
  if (type=="FST"){
    plot_df <- plot_df %>% left_join(mean_fst_values,by="weighted_fst_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_fst)
    title1 = "Decay of Relative Prediction Accuracy Across Fst Pools with LDpred"
    xlabel1 = "Group Average Fst"
    subtitle1 = "Each Fst group is represented along the x-axis by the within-group Fst average."
    caption1 = "Total of 8 groups of interval size 0.02, with lowest and highest Fst values at -0.02 and 0.14, respectively.\n
         Note: The fifth group contained less than 1000 observations; all other groups n>2500."
  }
  
  if (type=="PC") {
    plot_df <- plot_df %>% left_join(mean_pc_values,by="PC_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_pc)
    title1 = "Decay of Relative Prediction Accuracy Across PC Distance Pools with LDpred"
    xlabel1 = "Group Average PC Distance"
    subtitle1 = "Each PC distance group is represented along the x-axis by the within-group PC distance from European training population average."
    caption1 = range_pc
  }
  
  plot_df %>%
    subset(group_number != max(plot_df$group_number)) %>%
    ggplot(aes(x=mean_values,y=relative_performance, 
               color=phenotype))+
    geom_point(alpha=0.3, size=3)+ 
    geom_smooth(data=plot_df %>% filter(group_number!=1),
                method="lm",se=T,color="black",na.rm=T)+
    guides(shape = guide_legend(title="Threshold",
                                labels=levels(plot_df$threshold)),
           color = guide_legend(title="Phenotypes"),
           size=F)+
    theme_light() + 
    xlab(xlabel1) +
    ylab(expression(Partial~R^2~relative~to~Europeans))+
    labs(title=title1,
         subtitle=subtitle1,
         caption=caption1)+
    ylim(0, 1)
  
  
}

fst_ld_graph <- rel_partial_LDpred(ld_fst_df,"FST")
ggsave('img/LD_FST.png', fst_ld_graph, width = 6, height = 10, dpi = 300)


pc_ld_graph <- rel_partial_LDpred(ld_pc_df,"PC")
ggsave('img/LD_PC.png', pc_ld_graph, width = 6, height = 10, dpi = 300)



combined <- function(prs_df,ldpred_prs_df, type){
  

  ld_df <- ldpred_prs_df %>%
    filter(threshold==28)
  
  plot_df <- rbind(prs_df,ld_df)
  
  plot_df$threshold <- plot_df$threshold %>% as.factor()
  levels(plot_df$threshold) <- c("5e-8","1e-5","1e-4","1e-3","1e-2", "LDpred2")
  
  
  plot_df <- plot_df %>%
    filter(threshold!="5e-8") %>%
    group_by(phenotype,threshold) %>%
    mutate(
      eur_partial = max((partial*(group_number==1)),na.rm=T),
      relative_performance = partial / eur_partial) %>%
    ungroup() 
  
  
  if (type=="FST"){
    plot_df <- plot_df %>% left_join(mean_fst_values,by="weighted_fst_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_fst)
    title1 = "Decay of Relative Prediction Accuracy Across Fst Pools"
    xlabel1 = "Group Average Fst"
    subtitle1 = "Each Fst group is represented along the x-axis by the within-group Fst average."
    caption1 = "Total of 8 groups of interval size 0.02, with lowest and highest Fst values at -0.02 and 0.14, respectively.\n
         Note: The fifth group contained less than 1000 observations; all other groups n>2500."
  }
  if (type=="PC") {
    plot_df <- plot_df %>% left_join(mean_pc_values,by="PC_groups")
    plot_df <- plot_df %>% rename(mean_values = mean_pc)
    title1 = "Decay of Relative Prediction Accuracy Across PC Distance Pools"
    xlabel1 = "Group Average PC Distance"
    subtitle1 = "Each PC distance group is represented along the x-axis by the within-group PC distance from European training population average."
    caption1 = range_pc
  }

  
  plot_df <- plot_df %>%
    filter(group_number != max(plot_df$group_number))
  
  
  ggplot(plot_df, aes(x=mean_values,y=as.numeric(relative_performance), 
                      color=phenotype))+
    geom_point(size=3,alpha=0.3)+
    geom_smooth(data=plot_df %>% filter(group_number!=1),
                method="lm",se=T,color="black",na.rm=T)+
    guides(color = guide_legend(title="Phenotypes"),
           size=F)+
    theme_light() + 
    scale_fill_brewer(palette="Set3")+
    xlab(xlabel1) +
    ylab(expression(Partial~R^2~relative~to~Europeans))+
    labs(title=title1,
         subtitle=subtitle1,
         caption=caption1)+
    #theme(legend.position=c(.90,.80))+
    ylim(0, 1)+
    facet_grid(cols = vars(threshold))
  
}

comb_fst <- combined(prs_fst_df,ld_fst_df, "FST")
ggsave('img/COMB_FST.png',comb_fst,width=20,height=10)
comb_pc <- combined(prs_pc_df,ld_pc_df, "PC")
ggsave('img/COMB_PC.png',comb_pc,width=20,height=10)

comparison <- function(prs_df,ldpred_prs_df, type){
  
  #  prs_df$threshold <- (prs_df$threshold-6) %>% as.factor()
  #levels(prs_df$threshold) <- c("5e-8","1e-6","1e-4","1e-3","1e-2")
  #  levels(prs_df$threshold) <- c(1,round(seq_log(1e-4,0.9,20),4),"Average LDpred p")
  
  thresholds <- prs_df %>%
    filter(group_number==1) %>%
    group_by(phenotype) %>%
    slice(which.max(partial)) %>%
    select(phenotype,threshold)
  
  plot_df1 <- data.frame()
  for(i in 1:10){
    sub_df <- prs_df %>%
      filter(phenotype==thresholds$phenotype[i],
             threshold==thresholds$threshold[i])
    plot_df1 <- rbind(plot_df1,sub_df)
    
  }
  
  plot_df2 <- ldpred_prs_df %>% filter(threshold==28)

  
  if (type=="FST"){
    plot_df1 <- plot_df1 %>% left_join(mean_fst_values,by="weighted_fst_groups")
    plot_df1 <- plot_df1 %>% rename(mean_values = mean_fst)
    title1 = "Improvement of Prediction Accuracy Across Fst Pools with LDpred"
    xlabel1 = "Group Average Fst"
    subtitle1 = "Each Fst group is represented along the x-axis by the within-group Fst average."
    caption1 = "Total of 8 groups of interval size 0.02, with lowest and highest Fst values at -0.02 and 0.14, respectively.\n
         Note: The fifth group contained less than 1000 observations; all other groups n>2500."
    new_df <- plot_df1 %>% inner_join(plot_df2, by=c("weighted_fst_groups",
                                                     "phenotype")) %>%
      mutate(diff_partial=(partial.y/partial.x))
    }
  
  if (type=="PC") {
    plot_df1 <- plot_df1 %>% left_join(mean_pc_values,by="PC_groups")
    plot_df1 <- plot_df1 %>% rename(mean_values = mean_pc)
    title1 = "Improvement of Prediction Accuracy Across PC Distance Pools with LDpred"
    xlabel1 = "Group Average PC Distance"
    subtitle1 = "Each PC distance group is represented along the x-axis by the within-group PC distance from European training population average."
    caption1 = range_pc
    
    new_df <- plot_df1 %>% inner_join(plot_df2, by=c("PC_groups",
                                                     "phenotype")) %>%
      mutate(diff_partial=(partial.y/partial.x))
  }
  
  
  new_df %>%
    ggplot(aes(x = mean_values,
               y = diff_partial, color = phenotype)) +
    geom_point(size=3,alpha=0.3) +
    geom_smooth(data=new_df,
                method="lm",se=T,color="black",na.rm=T)+
    theme_classic() + 
    xlab(xlabel1) +
    ylab(expression(Ratio~of~LDpred~and~C+T~Partial~R^2))+
    scale_y_log10() +
    labs(title=title1,
         caption=caption1,
         subtitle=subtitle1)+
    ylim(0,10)
  
}

comp_fst <- comparison(prs_fst_df,ld_fst_df, "FST")
ggsave('img/COMP_FST.png',comp_fst,width=6,height=10)
comp_pc <- comparison(prs_pc_df,ld_pc_df, "PC")
ggsave('img/COMP_PC.png',comp_pc,width=6,height=10)

