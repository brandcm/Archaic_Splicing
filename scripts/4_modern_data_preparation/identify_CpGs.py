# Colin M. Brand, University of California San Francisco, 11/02/2022

import argparse
import pysam

def parse_args():
        parser = argparse.ArgumentParser()

        parser.add_argument(
                "--variants", type=str, required=True, help="Path to input text file with variants using 1-based coordinates.")

        parser.add_argument("--fasta", type=str, required=True, help="Path to input sequence in FASTA format.")

        parser.add_argument("--out", type=str, required=True, help="Path to output file.")

        args = parser.parse_args()
        return args

def main():
	args = parse_args()

        with open(f'{args.variants}', 'r') as variants, open(f'{args.out}', 'w') as out:
                variant_lines = [ variant_line.split() for variant_line in variants ]
                for variant_line in variant_lines:
                        if variant_line[2] == 'C':
                                sequence = pysam.FastaFile(args.fasta).fetch(variant_line[0], int(variant_line[1]) - 1, int(variant_line[1]) + 1).upper()
                                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(variant_line[0], variant_line[1], variant_line[2], variant_line[3], variant_line[4], variant_line[5], variant_line[6], sequence))
                        else:
                             	out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(variant_line[0], variant_line[1], variant_line[2], variant_line[3], variant_line[4], variant_line[5], variant_line[6], 'n/a'))

if __name__ == '__main__':
    main()
