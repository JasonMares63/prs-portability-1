import pickle

import pandas as pd
import sklearn.ensemble


def load_1000_genomes(pc_cols):
    """
    Loads training data from the 1000 Genomes Project. Individuals are labelled
    only by their population codes, so a map must be used to aggregate
    populations into super populations.
    """
    # 1000 Genomes principal components file
    kgp_pc_df = pd.read_csv(
        'data/kgp_merged/projection.sscore',
        sep='\t',
        usecols=['#FID', 'IID', *pc_cols]
    )
    # 1000 Genomes individual population labels
    kgp_pop_df = pd.read_csv(
        '/rigel/mfplab/projects/prs-portability/data/kgp_meta/integrated_call_samples_v2.20130502.ALL.ped',
        sep='\t',
        usecols=['Family ID', 'Individual ID', 'Gender', 'Population']
    )
    # 1000 Genomes map between population and super population
    kgp_super_pop_df = pd.read_csv(
        '/rigel/mfplab/projects/prs-portability/data/kgp_meta/20131219.populations.tsv',
        sep='\t',
        usecols=['Population Code', 'Super Population']
    )

    # Combine DataFrames to train classifier
    kgp_merged = (
        kgp_pop_df
        .merge(kgp_super_pop_df, left_on='Population',
               right_on='Population Code', how='inner')
        .merge(kgp_pc_df, left_on='Individual ID', right_on='IID', how='inner')
        .filter(items=['Family ID','IID', 'Population', 'Super Population',
                       'Gender', *pc_cols])
    )
    return kgp_merged


def train_classifier(df, pc_cols, n_estimators=100, seed=100):
    # Separate data into X and y, as typical in the sklearn API
    principal_components = df.loc[:, pc_cols].values
    labels = df.loc[:, 'Super_Population'].values

    rf = sklearn.ensemble.RandomForestClassifier(n_estimators=n_estimators,
                                                 random_state=seed)
    rf.fit(principal_components, labels)
    return rf


def classify_ukbb(classifier, pc_cols, min_prob=0.9):
    # Load UK Biobank PCA data and store as two arrays (to save memory)
    ukb_pc_df = pd.read_csv(
        'data/ukb_merged/projection.sscore',
        sep='\t',
        usecols=['#FID', 'IID', *pc_cols]
    )
    sample_ids = ukb_pc_df.loc[:, ['#FID', 'IID']]
    pc_array = ukb_pc_df.loc[:, pc_cols].values
    del ukb_pc_df

    predicted = classifier.predict(pc_array)
    probabilities = classifier.predict_proba(pc_array)
    max_probabilities = probabilities.max(axis=1)
    inconclusive = (max_probabilities < min_prob)

    ukb_labels_df = (
        sample_ids
        .assign(
            predicted=predicted,
            max_probability=max_probabilities,
            inconclusive=inconclusive,
        )
    )
    del predicted, max_probabilities, inconclusive

    for i, label in enumerate(classifier.classes_):
        ukb_labels_df[label] = probabilities[:, i]
    return ukb_labels_df


if __name__ == '__main__':

    use_pcs = [f'PC{i}_AVG' for i in range(1, 11)]

    # Load 1000 Genomes data (only those PCs we want to use)
    kgp_df = load_1000_genomes(use_pcs)

    #Create population files for 1KG data
    kgp_df = kgp_df.rename(columns={'Family ID':'#FID'}) 
    kgp_df = kgp_df.rename(columns={'Super Population':'Super_Population'})

    for population in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        (
            kgp_df
            .query(f'Super_Population == "{population}"')
            .loc[:, ['#FID', 'IID']]
            .to_csv(f'data/kgp_populations/{population}_all.txt', header=True,
                    index=False, sep=' ')
        )
    # Train classifier and then proceed to implement on UKBB
    classifier = train_classifier(kgp_df, use_pcs)
    del kgp_df

    with open('data/models/population_classifier.pkl', 'wb') as f:
        pickle.dump(classifier, f)

    ukb_labels_df = classify_ukbb(classifier, use_pcs)
    ukb_labels_df.to_csv('data/ukb_merged/population_labels_10PCS.tsv.gz',
                         sep='\t', compression='infer', index=False)
