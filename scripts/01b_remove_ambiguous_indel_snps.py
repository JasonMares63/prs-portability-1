import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Remove SNPs from a .snps file'))
    parser.add_argument('path', help='path to .snps file to be filtered')
    parser.add_argument('-r', '--remove', help=('path to .snps file for SNPs '
                                                'to remove from `path` file'))
    parser.add_argument('-o', '--output',
                        help='Output directory for final .snps file')
    args = parser.parse_args()

    with open(args.path, 'r') as f:
        input_snps = set(f.read().splitlines())

    with open(args.remove, 'r') as f:
        snps_to_remove = set(f.read().splitlines())

    remaining_snps = sorted(input_snps.difference(snps_to_remove))

    with open(args.output, 'w') as f:
        f.write('\n'.join(remaining_snps))
