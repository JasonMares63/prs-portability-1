import argparse
import csv
import pathlib


parser = argparse.ArgumentParser(description=('Remove ambiguous variants or indels from a .bin file'))
parser.add_argument('path', help='path to the .bim file')
parser.add_argument('-o', '--output', help='Output directory for filtered .snps file')

args = parser.parse_args()
path = pathlib.Path(args.path)
if not path.exists():
    raise ValueError(f"Passed path {path} does not exist")
output_path = pathlib.Path(args.output)

excluded = {
    ('A', 'T'),
    ('T', 'A'),
    ('C', 'G'),
    ('G', 'C'),
}

unfiltered = open(path, 'r')
unfiltered_reader = csv.reader(unfiltered, delimiter='\t', skipinitialspace=True, quotechar='|')

filtered = open(output_path, 'w')
filtered_writer = csv.writer(filtered, delimiter='\t')

for row in unfiltered_reader:
    if len(row) < 4:
        raise ValueError(f"Wrong size row: {row}")
    # Filter indels
    if len(row[4]) > 1 or len(row[5]) > 1:
        filtered_writer.writerow([row[1]])
    # Filter "ambiguous" variants
    if (row[4], row[5]) in excluded:
        filtered_writer.writerow([row[1]])

unfiltered.close()
filtered.close()

