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
    Using this file in order to neatly retrieve trait names
    """
    martin_gwas_info = pd.read_csv('/rigel/mfplab/projects/prs-portability/data/martin_gwas_info.txt', sep=r'\s+')

    trait_to_n_samples = (
        martin_gwas_info
        .set_index('Trait')
        .loc[:, 'NGWAS']
        .to_dict()
    )
    return trait_to_n_samples


def make_trait_sample_files(all_samples_df, trait_to_n_samples, seed=100):
    """
    Sample 80000 of EUR population individuals for the GWAS
    of each trait. 
    """
    for trait, n_samples in trait_to_n_samples.items():
        (
            all_samples_df
            .sample(n=80000, random_state=seed)
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

    # Get trait names
    trait_to_n_samples = get_trait_to_n_samples()

    # Create a sample file for each trait (correct number of samples, all EUR)
    make_trait_sample_files(eur_train_df, trait_to_n_samples, seed=seed)

    # Create a sample file for LD
    ld_samples = eur_train_df.sample(n=40000, random_state=seed)
    ld_samples.to_csv(f'data/ukb_populations/LD_EUR_train.txt', header=True,index=False, sep=' ')

