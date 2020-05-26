import argparse

import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Filter SNPs from a GWAS output file by a p-value threshold.'))
    parser.add_argument('path', help='path to GWAS output file')
    parser.add_argument('-o', '--output',
                        help=('path to the SNPs file to be created. This file'
                              'has one variant identifier per line, only.'))

    args = parser.parse_args()

    (
        pd.read_csv(args.path, sep=r'\t')
        .loc[lambda df: df.duplicated(subset=['ID'], keep='first'), 'ID']
        .to_csv(args.output, index=False, header=False, sep='\t')
    )
