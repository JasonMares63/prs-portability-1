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

# Table of each trait's name, code, field ID, and dtype (including eid)
traits_info <- read_delim('data/martin_gwas_info.txt', delim = ' ',
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

raw_phenotypes_df <- read_csv('/rigel/mfplab/users/mnz2108/ukbiobank/data/ukb40732.csv',
                              col_types = do.call(cols_only, field_to_dtype))

# Samples table (because not assured that FID == IID)
psam_df <- read_tsv(
    'data/ukb_merged/merged.psam',
    col_types = cols_only('#FID' = col_character(), 'IID' = col_character())
)

raw_phenotypes_df %>%
    # Rename fields to their human-readable names
    pivot_longer(-eid, names_to = 'ukb_field') %>%
    inner_join(traits_info %>% select(-dtype, -UKBBcode), by = 'ukb_field') %>%
    pivot_wider(id_cols = eid, names_from = Trait, values_from = value) %>%
    rename(IID = eid) %>%
    # Apply inverse rank normal transformation
    mutate_at(.vars = vars(-IID), .funs = irnt) %>%
    # Add FID column using the sample file
    left_join(psam_df, by = 'IID') %>%
    select('#FID', IID, all_of(traits_info %>% filter(Trait != 'eid') %>% pull(Trait))) %>%
    write_delim('data/phenotypes/full_phenotypes.pheno', delim = ' ', na = 'NA')
