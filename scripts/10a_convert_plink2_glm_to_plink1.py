import argparse

import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Convert a Plink 2 association output (.glm.linear) format to '
        'Plink 1 association output format (.assoc)'))
    parser.add_argument('path', help=('path to a Plink 2 GWAS output file '
                                      '(.glm.linear or .glm.logistic) '
                                      'to be converted')
    parser.add_argument('-o', '--output',
                        help=('path to the Plink 1 output file to be created'))
    args = parser.parse_args()

    (
        pd.read_csv(args.path, sep='\t')
        .rename(columns={
            '#CHROM': '#CHR',
            'POS': 'BP',
            'ID': 'SNP',
            'OBS_CT': 'NMISS',
            'T_STAT': 'STAT',
        })
        .filter(items=['#CHR', 'SNP', 'BP', 'A1', 'TEST', 'NMISS', 'BETA',
                       'STAT', 'P', 'SE'])
        .to_csv(args.output, index=False, sep='\t', float_format='%.9g')
    )
