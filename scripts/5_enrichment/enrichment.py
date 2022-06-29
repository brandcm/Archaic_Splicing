import argparse
import numpy as np
import pandas as pd

def parse_args():
	parser = argparse.ArgumentParser()
		
	parser.add_argument(
		"--set_name", type=str, required=True, help="Name of gene set")

	return parser.parse_args()

def main():
	args = parse_args()
	
	empiric = {}
	for ont in ['GWAS','HPO']:
		empiric[ont] = []
		print(args.set_name, ont)
		observed = pd.read_csv(f'../observed/{args.set_name}_{ont}_observed.txt', sep='\t', index_col=0)
		observed_with_one_observation = (observed.iloc[1:] > 0).sum(1)
		file = open(f'{args.set_name}_{ont}_empiric_counts.txt')
		Lines = file.readlines()
		for line in Lines:
			line = line.strip().split("\t")
			if line[0] in observed_with_one_observation:
				obs = observed.loc[line[0]].sum()
				counts = [int(x) for x in line[1:]]
				mean = np.mean(counts)
				p = (sum(counts >= obs ) + 1) / (len(counts) + 1)
				empiric[ont].append({'label':line[0],'observed':obs,'mean_expected':mean,'enrichment':obs/mean,'p_value':p})
		empiric[ont] = pd.DataFrame(empiric[ont])
		empiric[ont].sort_values('p_value').to_csv(f'../enrichment/{args.set_name}_{ont}_enrichment.txt',sep="\t",header=True,index=False)
		
if __name__ == '__main__':
    main()
