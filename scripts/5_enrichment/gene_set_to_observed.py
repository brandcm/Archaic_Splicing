# This script takes a set of genes and designates whether each gene
# occurs (True) or not (False) per GWAS Catalog 2019 and Human Phenotype 
# Ontology term. These are output in two separate files with the first
# row as the list of genes and the subsequent rows indicating GWAS or
# HPO term followed by the boolean per gene.

import argparse

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--set_path", type=str, required=True, help="Path to list of genes")
		
	parser.add_argument(
		"--set_name", type=str, required=True, help="Name of gene set")

	return parser.parse_args()

def main():
	args = parse_args()
	
	f = open(args.set_path, 'r')
	genes = [line.rstrip('\n') for line in f]

	for ont in ['GWAS','HPO']:
		with open(f'../{ont}.txt') as file, open(f'{args.set_name}_{ont}_observed.txt', 'w') as out:
			unlisted_genes='\t'.join(genes)
			out.write("\t{}\n".format(unlisted_genes))
			for terms in file:
				terms = terms.strip().split("\t")
				term_name = terms[0]
				booleans = []
				for gene in genes:
					if gene in terms[2:]:
						booleans.append('True')
					else:
						booleans.append('False')
				unlisted_booleans='\t'.join(booleans)
				out.write("{}: {}\t{}\n".format(ont, term_name, unlisted_booleans))

if __name__ == '__main__':
    main()
