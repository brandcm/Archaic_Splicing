import argparse
import numpy as np
import pandas as pd

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--input", type=str, required=True, help="Path to input file of empirical counts")
		
	parser.add_argument(
		"--subset", type=int, help="Number of empiric counts to consider")
		
	parser.add_argument(
		"--output", type=str, required=True, help="Path to output")
		
	args = parser.parse_args()
	
	return args

def main():
	args = parse_args()
	f = open(args.input, 'r')
	terms = f.readlines()
	out_df = []
	for term in terms:
		line = np.array(term.strip().split('\t'))
		name = line[0]
		if args.subset is not None:
			line = line[1:args.subset+1]
		else:
			line = line[1:]
		pvals = []
		for i, v in enumerate(line):
			if i % 10 == 0:
				print(i,flush=True)
			p = sum(line[:i] >= v) + sum(line[i+1:] >= v)
			p = (p+1)/len(line)
			pvals.append(p)
		out_df.append({**{'name': name},**dict(zip(range(len(pvals)),pvals))})
	out_df = pd.DataFrame(out_df)
	out_df.to_csv(f'{args.output}.txt', sep = '\t', header = None, index = None)

if __name__ == '__main__':
    main()
