library(tidyverse)

irnt <- function(pheno) {
    # Inverse rank normal transformation. Modified from original in PHEASANT 
    # (https://github.com/MRCIEU/PHESANT/blob/d6508ac224b4e483b6d7f7530668b941f5746cdb
    #  /WAS/testContinuous.r#L243-L249)
    set.seed(1234)
    num_present <- sum(!is.na(pheno))
    ranks <-rank(pheno, na.last="keep", ties.method="random")
    quantilePheno = (ranks - 0.5) / num_present
    phenoIRNT = qnorm(quantilePheno)
    return(phenoIRNT);
}


combine_phenotype_tables <- function(phenotypes_vector) {
    all_pheno_df <- NULL
    for (pheno in phenotypes_vector) {
        pheno_df <- read_delim(
            str_glue('data/phenotypes/ukb.{pheno}.txt'), delim = ' ', col_names = c('IID', pheno), 
            col_types = c('IID' = col_integer(), pheno = col_double())
        )
        pheno_df <- pheno_df %>%
	    mutate(
                !!pheno := irnt(pheno_df %>% pull(pheno))
            )
        
        # Join each new phenotype on the existing data.frame, if it exists
        if (all_pheno_df %>% is.null) {
            all_pheno_df <- pheno_df
        } else {
            all_pheno_df <- all_pheno_df %>%
                left_join(pheno_df, by = 'IID')
        }
    }
    return(all_pheno_df)
}


phenotypes <- c('bmi', 'blood_pressure_dias', 'blood_pressure_sys','height')
all_pheno_df <- combine_phenotype_tables(phenotypes)

psam_df <- read_tsv(
    'data/ukb_filtered/merged.psam', 
    col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'SEX' = col_character())
)

all_pheno_df %>%
    left_join(psam_df, by = 'IID') %>%
    select('#FID', IID, BMI = bmi, DBP = blood_pressure_dias, 
           SBP = blood_pressure_sys, Height = height) %>%
    write_delim('data/phenotypes/full_phenotypes.pheno', delim = ' ', na = 'NA')

