import pathlib

import pandas as pd


pop_path = pathlib.Path('data/ukb_populations/')
population_files = [file for file in pop_path.glob('[!fst][!EUR]*_all.txt')
                    if 'EUR' not in file.stem]

# Anti-join EUR individuals to find the EUR train individuals (cluster ID=EUR)
eur_test = pd.read_csv(pop_path.joinpath('EUR_test.txt'), sep=r'\s+')
eur_all = pd.read_csv(pop_path.joinpath('EUR_all.txt'), sep=r'\s+')

clusters_df = (
    eur_all
    .merge(eur_test, how='outer', on=['#FID', 'IID'], indicator=True)
    .query('IID > 0')
    .query('_merge != "both"')
    .loc[:, ['#FID', 'IID']]
    .assign(cluster='EUR')
)

# Add other populations as single-individual clusters (cluster ID=IID)
for population_file in population_files:
    population_df = (
        pd.read_csv(population_file, sep=r'\s+')
        .query('IID > 0')
        .assign(cluster=lambda df: df['IID'])
    )
    clusters_df = pd.concat([clusters_df, population_df], sort=False)

clusters_df.to_csv('data/ukb_populations/fst_clusters.txt', sep=' ',
                   header=False, index=False)
