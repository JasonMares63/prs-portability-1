import pandas as pd


def normalize_covariates(df):
    covars = df.loc[:, 'sex_covar':'PC20_AVG']
    return_df = df.copy()
    return_df.loc[:, 'sex_covar':'PC20_AVG'] = (
        covars - covars.mean(axis=0)) / covars.std(axis=0)
    return return_df


if __name__ == '__main__':
    sex_df = pd.read_csv('data/phenotypes/ukb.sex.txt', sep=r'\s+',
                         header=None, names=['IID', 'is_male'])

    age_df = pd.read_csv('data/phenotypes/ukb.age.txt', sep=r'\s+',
                         header=None, names=['IID', 'age'])

    pc_df = pd.read_csv('data/ukb_merged/projection.sscore', sep=r'\s+')

    (
        sex_df
        .merge(age_df, on=['IID'], how='outer')
        .merge(pc_df, on=['IID'], how='outer')
        .drop(columns=['ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM'])
        .dropna(subset=['age'])
        .assign(
            # Code sex as 0 = missing, 1 = male, 2 = female, as in plink .sample files
            sex_covar=lambda df: (2 - df['is_male']).fillna(0),
            age_sq=lambda df: df['age']**2,
            age_sex=lambda df: df['age'] * df['sex_covar'],
            age_sq_sex=lambda df: df['age_sq'] * df['sex_covar'],
        )
        .filter(items=['#FID', 'IID', 'sex_covar', 'age', 'age_sq', 'age_sex',
                       'age_sq_sex', *[f'PC{i}_AVG' for i in range(1, 21)]])
        .pipe(normalize_covariates)
        .to_csv('data/ukb_merged/covar_all_samples.covar', sep='\t',
                na_rep='NA', header=True, index=False, float_format='%.7g')
    )
