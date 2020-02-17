import argparse
import pathlib


def parse_args(parser):
    args = parser.parse_args()
    input_path = pathlib.Path(args.path)
    if not input_path.exists():
        raise ValueError(f"Passed path {input_path} does not exist")
    remove_path = pathlib.Path(args.remove)
    if not remove_path.exists():
        raise ValueError(f"Passed path {input_path} does not exist")
    output_path = pathlib.Path(args.output)
    return input_path, remove_path, output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=(
        'Remove SNPs from a .snps file'))
    parser.add_argument('path', help='path to .snps file to be filtered')
    parser.add_argument('-r', '--remove', help='path to .snps file for SNPs to remove from `path` file')
    parser.add_argument('-o', '--output', help='Output directory for final .snps file')

    input_path, remove_path, output_path = parse_args(parser)

    with open(input_path, 'r') as f:
        input_snps = set(f.read().splitlines())

    with open(remove_path, 'r') as f:
        snps_to_remove = set(f.read().splitlines())

    remaining_snps = sorted(input_snps.difference(snps_to_remove))

    with open(output_path, 'w') as f:
        f.write('\n'.join(remaining_snps))

