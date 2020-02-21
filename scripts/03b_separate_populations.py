import pandas as pd


def make_population_sample_files(all_populations_df):
    """Make files with IIDs for each of the populations used for GWAS"""
    for population in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        (
            all_populations_df
            .query(f'predicted == "{population}"')
            .loc[:, ['#FID', 'IID']]
            .to_csv(f'data/ukb_populations/{population}_all.txt', header=True,
                    index=False, sep=' ')
        )


def get_trait_to_n_samples():
    """
    Each trait uses different number of EUR for GWAS (to enable comparison
    with BBJ). Sample these numbers from the EUR DataFrame and save them as
    sample files for each trait.
    """
    martin_gwas_info = pd.read_csv('data/martin_gwas_info.txt', sep=r'\s+')

    trait_to_n_samples = (
        martin_gwas_info
        .set_index('Trait')
        .loc[:, 'NGWAS']
        .to_dict()
    )
    return trait_to_n_samples


def make_trait_sample_files(all_samples_df, trait_to_n_samples, seed=100):
    """
    Sample the appropriate number of EUR population individuals for the GWAS
    of each trait. As per Martin et al.'s method, the number of individuals
    used for the GWAS is dependent on the number of individuals used in the
    equivalent-trait GWAS in Biobank Japan.
    """
    for trait, n_samples in trait_to_n_samples.items():
        (
            all_samples_df
            .sample(n=n_samples, random_state=seed)
            .to_csv(f'data/ukb_populations/{trait}.txt', header=True,
                    index=False, sep=' ')
        )


if __name__ == '__main__':
    seed = 100

    # DataFrame of predicted super population labels for UKBB samples
    labels_df = (
        pd.read_csv('data/ukb_merged/population_labels_10PCS.tsv.gz', sep='\t')
        # .query('inconclusive == False')
        .drop_duplicates(subset=['#FID', 'IID'])
    )

    # 5000 random samples used as the "Target" for EUR
    eur_test_df = (
        labels_df
        .query('predicted == "EUR"')
        .sample(n=5000, replace=False, random_state=seed)
    )

    (
        eur_test_df
        .loc[:, ['#FID', 'IID']]
        .to_csv('data/ukb_populations/EUR_test.txt', index=False, header=True,
                sep=' ')
    )

    # All labels for all populations (5000 test EUR removed)
    non_test_df = labels_df.drop(eur_test_df.index, inplace=False)

    # Make files with IIDs for each of the populations used for GWAS
    make_population_sample_files(non_test_df)

    # Use only EUR population for GWAS
    eur_train_df = (
        non_test_df
        .query('predicted == "EUR"')
        .loc[:, ['#FID', 'IID']]
    )
    del labels_df, non_test_df

    # Get number of samples for GWAS on each trait
    trait_to_n_samples = get_trait_to_n_samples()

    # Create a sample file for each trait (correct number of samples, all EUR)
    make_trait_sample_files(eur_train_df, trait_to_n_samples, seed=seed)
