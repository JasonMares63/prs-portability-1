library(tidyverse)
library(cowplot)
library(colorspace)


# Plot 1000 Genomes PCAs
labels_df <- read_tsv('data/integrated_call_samples_v2.20130502.ALL.ped',
                      col_types = cols_only('Individual ID' = col_character(),
                                            'Population' = col_character()))
super_pop_df <- read_tsv('data/20131219.populations.tsv')
labels_df <- labels_df %>%
    inner_join(super_pop_df, by = c('Population' = 'Population Code'))

pc_df <- read_tsv(
    'data/merged/projection.sscore', 
    col_types = cols_only(IID = col_character(), PC1_AVG = col_double(), PC2_AVG = col_double())
)

pc_plot <- pc_df %>%
    inner_join(labels_df, by = c('IID' = 'Individual ID')) %>%
    ggplot(aes(x = PC1_AVG, y = PC2_AVG, color = !!as.name('Super Population'))) +
    geom_point(alpha = 0.5) +
    xlab('PC1') +
    ylab('PC2')

ggsave(filename = 'img/1000_genomes_pca.png', plot = pc_plot, dpi = 300)


# Plot UK Biobank PCA (computed on 1000 genomes)
ukb_pca_df <- read_tsv('data/ukb_filtered/projection.sscore')
ukb_labels_df <- read_tsv('data/ukb_filtered/population_labels_10PCS.tsv.gz')

ukb_df <- ukb_pca_df %>%
    inner_join(ukb_labels_df, by = c('IID'))

rm(ukb_pca_df)
rm(ukb_labels_df)

ukb_pc_plot <- ukb_df %>%
    filter(IID > 0) %>%
    ggplot(aes(x = PC1_AVG, y = PC2_AVG, color = predicted, shape = inconclusive)) +
    geom_point(alpha = 0.5) +
    xlab('PC1') +
    ylab('PC2')

ggsave(filename = 'img/ukb_pca.png', plot = ukb_pc_plot, dpi = 300)


# Plot UK Biobank with conclusiveness separated
counts_df <- ukb_df %>%
    filter(IID > 0) %>%
    mutate(conclusive = inconclusive %>% as.integer %>% recode_factor('0' = 'Conclusive population label (P ≥ 0.9)', 
                                                                      '1' = 'Inconclusive population label (P < 0.9)')) %>%
    group_by(inconclusive, conclusive) %>%
    summarize(num = n()) %>%
    ungroup %>%
    mutate(
        num_f = format(num, big.mark = ",", scientific = FALSE),
        label = str_glue('N = {num_f}')
    )

ukb_pc_sep_plot <- ukb_df %>%
    filter(IID > 0) %>%
    mutate(conclusive = inconclusive %>% as.integer %>% recode_factor('0' = 'Conclusive population label (P ≥ 0.9)', 
                                                                      '1' = 'Inconclusive population label (P < 0.9)')) %>%
    ggplot(aes(x = PC1_AVG, y = PC2_AVG, color = predicted, shape = inconclusive)) +
    geom_point(alpha = 0.5) +
    xlab('PC1') +
    ylab('PC2') +
    facet_wrap(vars(conclusive)) +
    geom_text(data = counts_df, aes(label = label), x = -0.015, y = 0.075, color = 'black') +
    scale_shape_discrete(guide = "none") +
    scale_color_discrete(name = 'Predicted\npopulation')

ggsave(filename = 'img/ukb_pca_sep.png', plot = ukb_pc_sep_plot, dpi = 300, width = 12, height = 5)


# counts_df <- ukb_df %>%
#     filter(IID > 0) %>%
#     mutate(conclusive = inconclusive %>% as.integer %>% recode_factor('0' = 'P >= 0.9', '1' = 'P < 0.9')) %>% 
#     group_by(inconclusive, conclusive) %>%
#     summarize(num = n()) %>%
#     ungroup %>%
#     mutate(
#         num_f = format(num, big.mark = ",", scientific = FALSE),
#         label = str_glue('N = {num_f}')
#     )
# 
# ukb_pc_sep_plot <- ukb_df %>%
#     filter(IID > 0) %>%
#     mutate(conclusive = inconclusive %>% as.integer %>% recode_factor('0' = 'P >= 0.9', '1' = 'P < 0.9')) %>%
#     ggplot(aes(x = PC1_AVG, y = PC2_AVG, color = predicted, shape = inconclusive)) +
#     geom_point(alpha = 0.5) +
#     xlab('PC1') +
#     ylab('PC2') +
#     facet_wrap(vars(conclusive), labeller = label_parsed) + 
#     geom_text(data = counts_df, aes(label = label), x = -0.01, y = 0.075, color = 'black')
# 
# ggsave(filename = 'img/ukb_pca_sep.png', plot = ukb_pc_sep_plot, dpi = 300, width = 12, height = 5)


# Plot UK Biobank with densities for each super population
ukb_density_plot <- ukb_df %>%
    filter(IID > 0) %>%
    ggplot(aes(x = PC1_AVG, y = PC2_AVG)) +
    stat_density_2d(aes(fill = ..level..), 
                    geom = "polygon", color = darken("#0072B2", .2), size = 0.3, bins = 5) +
    scale_fill_gradient(low = desaturate(lighten("#0072B2", .9), .6),
                        high = desaturate(lighten("#0072B2", .6), .6), 
                        guide = "none") +
    facet_wrap(vars(predicted))
 
 ggsave(filename = 'img/ukb_pca_density.png', plot = ukb_density_plot, dpi = 300, width = 12, height = 5)


# Plot combined PCA plots
combined <- plot_grid(pc_plot + ggtitle('1000 Genomes') + theme(legend.position = 'NA'), 
                      ukb_pc_plot + ggtitle('UK Biobank') + theme(legend.position = 'NA'),
                      get_legend(ukb_pc_plot + scale_color_discrete(name  = 'Super-population') + scale_shape_discrete(name = "Inconclusive")), 
                      nrow = 1, labels = c('A', 'B', ''), rel_widths = c(1, 1, 0.25))

ggsave(filename = 'img/combined_pca.png', plot = combined, dpi = 300, width = 12, height = 5)


