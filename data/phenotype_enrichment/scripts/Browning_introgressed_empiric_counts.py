import argparse
import random
import pandas as pd

def parse_args():
	parser = argparse.ArgumentParser()
		
	parser.add_argument(
		"--iterations", type=int, required=True, default=1000, help="Number of observations to generate")
		
	parser.add_argument(
		"--ontology", type=str, required=True, help="Name of ontology")
	
	parser.add_argument(
		"--set_name", type=str, required=True, help="Name of variant set being considered")
	
	parser.add_argument(
		"--array_id", type=int, required=True, default=0, help="Array ID to use if parallelizing")
					
	args = parser.parse_args()
	
	return args

def main():
	args = parse_args()
	set = args.set_name
	ontology = get_ontology(set)
	ontology_empiric = []
	while len(ontology_empiric) < args.iterations:
	    gene_counts = shuffle_and_count(ontology)
	    if gene_counts is not None:
	    	ontology_empiric.append(gene_counts)
	ontology_empiric_df = pd.DataFrame(ontology_empiric).T
	ontology_empiric_df.to_csv(f'empiric_counts/{args.set_name}_{args.ontology}_empiric_counts_{args.array_id}.txt', header = False, index = True, sep = '\t')

def get_ontology(set):
	args = parse_args()
	observed_ontology_df = pd.read_csv(f'observed/{args.set_name}_{args.ontology}_observed.txt', sep='\t', index_col=0)
	terms = list(observed_ontology_df.index[observed_ontology_df.sum(axis=1) > 0])
	ontology = {}
	file = open(f'ontologies/{args.ontology}.txt', 'r')
	lines = file.readlines()
	for line in lines:
		line = line.strip().split("\t")
		if args.ontology + ": " + line[0] in terms:
			ontology[args.ontology + ": " + line[0]] = line[2:]
	file.close()
	return ontology

def shuffle_and_count(ontology):
	data = pd.read_csv('data/archaic_data_for_enrichment.txt', sep='\t', header=0)
	deltas = data['delta_max'].copy()
	random.shuffle(deltas)
	data['delta_max'] = deltas
	shuffled_genes = data.loc[(data['delta_max'] >= 0.2) & (data['Browning_allele_origin'] == 'introgressed'), 'annotation']
	gene_counts = {}
	for i,r in ontology.items():
		gene_counts[i] = 0
		for g in shuffled_genes:
			if g in r:
				gene_counts[i]+=1
	return gene_counts
	
if __name__ == '__main__':
	main()