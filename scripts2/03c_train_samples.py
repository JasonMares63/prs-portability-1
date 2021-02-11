
import pandas as pd

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
    of each trait.
    """
    for trait, n_samples in trait_to_n_samples.items():
        (
            all_samples_df
            .sample(n=80000, random_state=seed)
            .to_csv(f'data/ukb_populations/{trait}1.txt', header=True,
                    index=False, sep=' ')
        )

if __name__ == '__main__':
    seed = 100

    eur_train_df = pd.read_csv("data/ukb_populations/EUR_all.txt",sep=" ")
    trait_to_n_samples = get_trait_to_n_samples()
    make_trait_sample_files(eur_train_df, trait_to_n_samples, seed=seed)

    ld_samples = eur_train_df.sample(n=80000, random_state=seed)
    ld_samples.to_csv(f'data/ukb_populations/LD_EUR_train.txt', header=True,index=False, sep=' ')

