# Colin M. Brand, University of California San Francisco, 10/28/2022

import argparse
import numpy as np
import pysam

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument("--chromosome_lengths", type=str, required=True, help="Path to chromosome lengths file.")
	
	parser.add_argument("--exons", type=str, required=True, help="Path to exon coordinates in BED format.")
	
	parser.add_argument("--fasta", type=str, required=True, help="Path to input sequence in FASTA format.")
		
	parser.add_argument("--start_codons", type=str, required=True, help="Path to start codons file from GFF (in 1-based coordinates).")
	
	parser.add_argument(
		"--variants", type=str, required=True,
		help="Path to input text file with variants using 1-based coordinates. Text file should be tab-delimited and include the chromosome, position, and variant.")
	
	parser.add_argument("--out", type=str, required=True, help="Path to output FASTA file.")
	
	args = parser.parse_args()
	return args

def main():
	args = parse_args()

	with open(f'{args.variants}', 'r') as variants, open(f'{args.exons}', 'r') as exons, open(f'{args.start_codons}', 'r') as codons, open(f'{args.out}', 'w') as out:
		
		# write head of out file
		out.write('chrom\tpos\tref_allele\talt_allele\tgene\tstrand\tconsequence\tdelta_pos\texon_number\texonic_intronic\tsequence_effect\tcanonical_transcript_len\tnovel_transcript_len\tcanonical_len-novel_len\tn_amino_acids\tcanonical_start_codon\tnovel_start_codon\tn_premature_stops\tresult\n')
		
		# analyze variants
		exon_lines = [ exon_line.split() for exon_line in exons ] # read in the exons
		start_codon_lines = [ start_codon_line.split() for start_codon_line in codons ] # read in start codons
		variant_lines = [ variant_line.split() for variant_line in variants ] # read in the variants
		for variant_line in variant_lines:
			chromosome = variant_line[0]
			variant = int(variant_line[1]) - 1 # convert to 0-based
			
			global alternate_allele
			alternate_allele = variant_line[3]
			gene = variant_line[4] # identify the gene
			
			# build canonical_transcript transcript
			print(f'Building canonical transcript for {gene}.')
			#canonical_transcript = {} # create empty dictionary for canonical_transcript sequence
			canonical_transcript_list = [] # create empty list for canonical_transcript sequence
			exon_subset = [ exons for exons in exon_lines if gene == exons[4] ] # subset exons for the target gene
			for exon in exon_subset:
				exon_sequence = pysam.FastaFile(args.fasta).fetch(exon[0], int(exon[1]), int(exon[2])).upper() # retrieve sequence for exons
				#canonical_transcript[gene] = sequence # add sequence to dictionary
				canonical_transcript_list.append(exon_sequence) # add sequence to list
			
			global canonical_transcript
			
			if exon_subset[0][3] == '+':
				print(f'{gene} is on the forward strand.')
				canonical_transcript = ''.join(canonical_transcript_list) # join strings to create single sequence per gene
			elif exon_subset[0][3] == '-':
				print(f'{gene} is on the reverse strand.')
				forward_canonical_transcript = ''.join(canonical_transcript_list)
				forward_canonical_transcript = forward_canonical_transcript.replace("A", "t").replace("C", "g").replace("G", "c").replace("T", "a") # complement
				forward_canonical_transcript = forward_canonical_transcript.upper() # convert to upper
				canonical_transcript = forward_canonical_transcript[::-1] # reorient
			
			# get chromosome lengths
			chromosome_lengths = {}
				
			with open(f'{args.chromosome_lengths}', 'r') as lengths:
				length_lines = [ length_line.split() for length_line in lengths ]
				for length_line in length_lines:
					chromosome_lengths[length_line[0]] = int(length_line[1])
			
			# define chromosome sequence for downstream functions
			global sequence
			
			if exon_subset[0][3] == '+':
				sequence = pysam.FastaFile(args.fasta).fetch(exon[0], 0, (int(chromosome_lengths[chromosome])-1)).upper()
			elif exon_subset[0][3] == '-':
				sequence = pysam.FastaFile(args.fasta).fetch(exon[0], 0, (int(chromosome_lengths[chromosome])-1)).upper()[::-1] # reorient
				sequence = sequence.replace("A", "t").replace("C", "g").replace("G", "c").replace("T", "a") # complement
				sequence = sequence.upper() # conver to upper
			
			# build novel transcript
			print(f'Building novel transcript for {gene}.')
			
			# forward strand genes
			if exon_subset[0][3] == '+':
				
				global exon_starts
				exon_starts = [ int(i[1]) for i in exon_subset ]
				
				global exon_stops
				exon_stops = [ int(i[2]) for i in exon_subset ]
				
				exon_lengths = [ exon_stops[i] - exon_starts[i] for i in range(len(exon_stops)) ]
				
				global exon_cumulative_lengths
				exon_cumulative_lengths = np.asarray(exon_lengths).cumsum().tolist()
				
				exon_stops = [i - 1 for i in exon_stops] # subtract one to convert from 0-based half-open to 0-based fully-closed
			
				
				if float(variant_line[5]) > float(variant_line[6]) and float(variant_line[5]) > float(variant_line[7]) and float(variant_line[5]) > float(variant_line[8]):
					print(f'The SAV at {variant_line[0]}: {variant_line[1]} results in an acceptor gain.')
					consequence = 'acceptor gain'
					delta_pos = int(variant_line[9])
					exon_number, exonic_intronic, sequence_effect, novel_transcript = acceptor_gain(variant, delta_pos)
			
				# acceptor loss	
				elif float(variant_line[6]) > float(variant_line[5]) and float(variant_line[6]) > float(variant_line[7]) and float(variant_line[6]) > float(variant_line[8]):
					print(f'The SAV at {variant_line[0]}: {variant_line[1]} results in an acceptor loss.')
					consequence = 'acceptor loss'
					delta_pos = int(variant_line[10])
					exon_number, exonic_intronic, sequence_effect, novel_transcript = acceptor_loss(variant, delta_pos)			

				# donor gain
				elif float(variant_line[7]) > float(variant_line[5]) and float(variant_line[7]) > float(variant_line[6]) and float(variant_line[7]) > float(variant_line[8]):
					print(f'The SAV at {variant_line[0]}: {variant_line[1]} results in a donor gain.')
					consequence = 'donor gain'
					delta_pos = int(variant_line[11])
					exon_number, exonic_intronic, sequence_effect, novel_transcript = donor_gain(variant, delta_pos)
			
				# donor loss
				elif float(variant_line[8]) > float(variant_line[5]) and float(variant_line[8]) > float(variant_line[6]) and float(variant_line[8]) > float(variant_line[7]):
					print(f'The SAV at {variant_line[0]}: {variant_line[1]} results in a donor loss.')
					consequence = 'donor loss'
					delta_pos = int(variant_line[12])
					exon_number, exonic_intronic, sequence_effect, novel_transcript = donor_loss(variant, delta_pos)			
			
				print(f'\n')
				
			# reverse strand genes
			elif exon_subset[0][3] == '-':
			
				# let's reorient the variant and exon starts/stops so that we can reapply the above code to reverse strand genes			
				exon_starts = [ int(chromosome_lengths[chromosome]) - int(i[2]) for i in exon_subset ]
				exon_starts = exon_starts[::-1]
				
				exon_stops = [ int(chromosome_lengths[chromosome]) - int(i[1]) for i in exon_subset ]
				exon_stops = exon_stops[::-1]
				
				exon_lengths = [ exon_stops[i] - exon_starts[i] for i in range(len(exon_stops)) ]
				exon_cumulative_lengths = np.asarray(exon_lengths).cumsum().tolist()
				
				exon_stops = [i - 1 for i in exon_stops] # subtract one to convert from 0-based half-open to 0-based fully-closed
				
				# reorient the variant
				inverse_variant = int(chromosome_lengths[chromosome]) - variant	- 1 # subtract one to get from 1-based to 0-based

				if float(variant_line[5]) > float(variant_line[6]) and float(variant_line[5]) > float(variant_line[7]) and float(variant_line[5]) > float(variant_line[8]):
					print(f'The SAV at {variant_line[0]}: {variant_line[1]} results in an acceptor gain.')
					consequence = 'acceptor gain'
					delta_pos = int(variant_line[9]) * -1 # reorient delta pos for negative strand
					exon_number, exonic_intronic, sequence_effect, novel_transcript = acceptor_gain(inverse_variant, delta_pos)

				# acceptor loss	
				elif float(variant_line[6]) > float(variant_line[5]) and float(variant_line[6]) > float(variant_line[7]) and float(variant_line[6]) > float(variant_line[8]):
					print(f'The SAV at {variant_line[0]}: {variant_line[1]} results in an acceptor loss.')
					consequence = 'acceptor loss'
					delta_pos = int(variant_line[10]) * -1 # reorient delta pos for negative strand
					exon_number, exonic_intronic, sequence_effect, novel_transcript = acceptor_loss(inverse_variant, delta_pos)			

				# donor gain
				elif float(variant_line[7]) > float(variant_line[5]) and float(variant_line[7]) > float(variant_line[6]) and float(variant_line[7]) > float(variant_line[8]):
					print(f'The SAV at {variant_line[0]}: {variant_line[1]} results in a donor gain.')
					consequence = 'donor gain'
					delta_pos = int(variant_line[11]) * -1 # reorient delta pos for negative strand
					exon_number, exonic_intronic, sequence_effect, novel_transcript = donor_gain(inverse_variant, delta_pos)
			
				# donor loss
				elif float(variant_line[8]) > float(variant_line[5]) and float(variant_line[8]) > float(variant_line[6]) and float(variant_line[8]) > float(variant_line[7]):
					print(f'The SAV at {variant_line[0]}: {variant_line[1]} results in a donor loss.')
					consequence = 'donor loss'
					delta_pos = int(variant_line[12]) * -1 # reorient delta pos for negative strand
					exon_number, exonic_intronic, sequence_effect, novel_transcript = donor_loss(inverse_variant, delta_pos)			
			
				print(f'\n')

			# get start codon position
			codon_line = [ codons for codons in start_codon_lines if gene == codons[5] ]
			
			global start_codon
			global relative_codon_start
			
			if codon_line:
				if exon_subset[0][3] == '+':
					start_codon = int(codon_line[0][2]) - 1
					start_codon_upstream_exon_start = int(next(start for start in reversed(exon_starts) if start_codon >= int(start)))
					upstream_exon_start_index = exon_starts.index(start_codon_upstream_exon_start)
			
				elif exon_subset[0][3] == '-':
					start_codon = int(chromosome_lengths[chromosome]) - int(codon_line[0][3])
					start_codon_upstream_exon_start = int(next(start for start in reversed(exon_starts) if start_codon >= int(start)))
					upstream_exon_start_index = exon_starts.index(start_codon_upstream_exon_start)

				if upstream_exon_start_index == 0:
					relative_codon_start = (start_codon - start_codon_upstream_exon_start)
				elif upstream_exon_start_index > 0:
					relative_codon_start = (start_codon - start_codon_upstream_exon_start) + exon_cumulative_lengths[upstream_exon_start_index - 1]
			
			else:
				start_codon = 'not annotated'
				relative_start_codon = 'not annotated'
				
			
			# estimate transcript differences
			print(len(canonical_transcript))
			print(len(novel_transcript))
			seq_diff, n_amino_acids, canonical_start_codon, novel_start_codon, n_premature_stops, result = compare(canonical_transcript, novel_transcript)
			
			# write out file
			out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chromosome, variant + 1, variant_line[2], alternate_allele, gene, exon_subset[0][3], consequence, delta_pos, exon_number, exonic_intronic, sequence_effect, len(canonical_transcript), len(novel_transcript), seq_diff, n_amino_acids, canonical_start_codon, novel_start_codon, n_premature_stops, result))


### Acceptor Gain ###

def acceptor_gain(variant, variant_effect):

	# exonic variant
	if any(lower <= variant <= upper for (lower, upper) in zip(exon_starts, exon_stops)):
		exonic_intronic = 'exonic'
		upstream_exon_start = int(next(start for start in reversed(exon_starts) if variant >= int(start)))
		upstream_exon_start_index = exon_starts.index(upstream_exon_start)
		exon_number = upstream_exon_start_index + 1

		if (variant + variant_effect) >= variant:
			sequence_effect = 'downstream effect, variant included, sequence deletion, no effect due to upstream acceptor'
			new_acceptor_gain_upstream_exon_start_distance = (variant + variant_effect) - upstream_exon_start
			novel_transcript = canonical_transcript[:(exon_cumulative_lengths[upstream_exon_start_index - 1] + new_acceptor_gain_upstream_exon_start_distance)] + alternate_allele + canonical_transcript[(exon_cumulative_lengths[upstream_exon_start_index - 1] + new_acceptor_gain_upstream_exon_start_distance + 1):]
	
		elif upstream_exon_start <= (variant + variant_effect) < variant:
			sequence_effect = 'upstream effect, variant included, sequence deletion, no effect due to upstream acceptor'
			new_acceptor_gain_upstream_exon_start_distance = (variant + variant_effect) - upstream_exon_start
			novel_transcript = canonical_transcript[:(exon_cumulative_lengths[upstream_exon_start_index - 1] + new_acceptor_gain_upstream_exon_start_distance)] + alternate_allele + canonical_transcript[(exon_cumulative_lengths[upstream_exon_start_index - 1] + new_acceptor_gain_upstream_exon_start_distance + 1):]
			
		elif (variant + variant_effect) < upstream_exon_start:
			sequence_effect = 'upstream effect, variant included, sequence addition'
			new_acceptor_gain_upstream_exon_start_distance = upstream_exon_start - (variant + variant_effect)
			insertion_position = exon_cumulative_lengths[upstream_exon_start_index - 1]
			variant_position = (variant - exon_starts[upstream_exon_start_index]) + exon_cumulative_lengths[upstream_exon_start_index - 1]
			novel_sequence = sequence[(upstream_exon_start - new_acceptor_gain_upstream_exon_start_distance):(upstream_exon_start - 1)]
			novel_transcript = canonical_transcript[:insertion_position] + novel_sequence + canonical_transcript[insertion_position:variant_position] + alternate_allele + canonical_transcript[(variant_position)+1:]
	
	# intronic variant
	else:
		exonic_intronic = 'intronic'
		downstream_exon_start = int(next(start for start in exon_starts if int(start) > variant))
		downstream_exon_start_index = exon_starts.index(downstream_exon_start)
		exon_number = downstream_exon_start_index + 1
		
		if (variant + variant_effect) >= downstream_exon_start:
			sequence_effect = 'downstream effect, variant excluded, sequence deletion, no effect due to upstream acceptor'
			novel_transcript = canonical_transcript
	
		elif variant <= (variant + variant_effect) < downstream_exon_start:
			sequence_effect = 'downstream effect, variant excluded, sequence addition'
			new_acceptor_gain_downstream_exon_start_distance = downstream_exon_start - (variant + variant_effect)
			insertion_position = exon_cumulative_lengths[downstream_exon_start_index - 1]
			novel_sequence = sequence[(downstream_exon_start - new_acceptor_gain_downstream_exon_start_distance):downstream_exon_start]
			novel_transcript = canonical_transcript[:insertion_position] + novel_sequence + canonical_transcript[insertion_position:]
		
		elif (variant + variant_effect) < variant:
			sequence_effect = 'upstream effect, variant included, sequence addition'
			variant_downstream_exon_start_distance = downstream_exon_start - variant
			new_acceptor_gain_downstream_exon_start_distance = downstream_exon_start - (variant + variant_effect)
			insertion_position = exon_cumulative_lengths[downstream_exon_start_index - 1]
			novel_sequence_1 = sequence[(downstream_exon_start - new_acceptor_gain_downstream_exon_start_distance):(downstream_exon_start - variant_downstream_exon_start_distance)]
			novel_sequence_2 = sequence[(downstream_exon_start - variant_downstream_exon_start_distance):(downstream_exon_start - 1)]
			novel_transcript = canonical_transcript[:insertion_position] + novel_sequence_1 + alternate_allele + novel_sequence_2 + canonical_transcript[insertion_position:]

	return exon_number, exonic_intronic, sequence_effect, novel_transcript

### Acceptor Loss ###

def acceptor_loss(variant, variant_effect):

	# exonic variant
	if any(lower <= variant <= upper for (lower, upper) in zip(exon_starts, exon_stops)):
		exonic_intronic = 'exonic'
		upstream_exon_start = int(next(start for start in reversed(exon_starts) if variant >= int(start)))
		upstream_exon_start_index = exon_starts.index(upstream_exon_start)
		exon_number = upstream_exon_start_index + 1
		
		if (variant + variant_effect) == upstream_exon_start:
			sequence_effect = 'exon skipped, variant excluded'
			novel_transcript = canonical_transcript[:exon_cumulative_lengths[upstream_exon_start_index - 1]] + canonical_transcript[exon_cumulative_lengths[upstream_exon_start_index]:]
			
		else:
			sequence_effect = 'exon not skipped, variant included'
			variant_upstream_exon_start_distance = variant - upstream_exon_start
			novel_transcript = canonical_transcript[:(exon_cumulative_lengths[upstream_exon_start_index - 1] + variant_upstream_exon_start_distance)] + alternate_allele + canonical_transcript[(exon_cumulative_lengths[upstream_exon_start_index - 1] + variant_upstream_exon_start_distance + 1):]
		
	# intronic variant
	else:
		exonic_intronic = 'intronic'
		downstream_exon_start = int(next(start for start in exon_starts if int(start) >= variant))
		downstream_exon_start_index = exon_starts.index(downstream_exon_start)
		exon_number = downstream_exon_start_index + 1
	
		if (variant + variant_effect) == downstream_exon_start:
			sequence_effect = 'exon skipped, variant excluded'
			novel_transcript = canonical_transcript[:exon_cumulative_lengths[downstream_exon_start_index - 1]] + canonical_transcript[exon_cumulative_lengths[downstream_exon_start_index]:]
		
		else:
			novel_transcript = canonical_transcript
			sequence_effect = 'exon not skipped, variant excluded'

	return exon_number, exonic_intronic, sequence_effect, novel_transcript

### Donor Gain ###

def donor_gain(variant, variant_effect):

	# exonic variant
	if any(lower <= variant <= upper for (lower, upper) in zip(exon_starts, exon_stops)):
		exonic_intronic = 'exonic'
		downstream_exon_stop = int(next(stop for stop in exon_stops if int(stop) >= variant))
		downstream_exon_stop_index = exon_stops.index(downstream_exon_stop)
		exon_number = downstream_exon_stop_index + 1
		
		if (variant + variant_effect) < variant:
			sequence_effect = 'upstream effect, variant excluded, sequence deletion'
			new_donor_gain_downstream_exon_stop_distance = downstream_exon_stop - (variant + variant_effect)
			deletion_position_start = exon_cumulative_lengths[downstream_exon_stop_index] - new_donor_gain_downstream_exon_stop_distance
			deletion_position_end = exon_cumulative_lengths[downstream_exon_stop_index]
			novel_transcript = canonical_transcript[:deletion_position_start] + canonical_transcript[deletion_position_end:]

			
		elif variant <= (variant + variant_effect) <= downstream_exon_stop:
			sequence_effect = 'downstream effect, variant included, sequence deletion'
			new_donor_gain_downstream_exon_stop_distance = downstream_exon_stop - (variant + variant_effect)
			deletion_position_start = exon_cumulative_lengths[downstream_exon_stop_index] - new_donor_gain_downstream_exon_stop_distance
			deletion_position_end = exon_cumulative_lengths[downstream_exon_stop_index]
			variant_position = exon_cumulative_lengths[downstream_exon_stop_index] - (exon_stops[downstream_exon_stop_index] - variant + 1)			
			novel_transcript = canonical_transcript[:variant_position] + alternate_allele + canonical_transcript[(variant_position + 1):deletion_position_start] + canonical_transcript[deletion_position_end:]
		
		elif (variant + variant_effect) > downstream_exon_stop:
			sequence_effect = 'downstream effect, variant included, sequence addition'
			new_donor_gain_downstream_exon_stop_distance = (variant + variant_effect) -  downstream_exon_stop
			insertion_position = exon_cumulative_lengths[downstream_exon_stop_index]
			variant_position = exon_cumulative_lengths[downstream_exon_stop_index] - (exon_stops[downstream_exon_stop_index] - variant + 1)
			novel_sequence = sequence[(downstream_exon_stop + 1):((downstream_exon_stop + 1) + new_donor_gain_downstream_exon_stop_distance)]
			novel_transcript = canonical_transcript[:variant_position] + alternate_allele + canonical_transcript[(variant_position + 1):insertion_position] + novel_sequence + canonical_transcript[insertion_position:]
	
	# intronic variant		
	else:
		exonic_intronic = 'intronic'
		upstream_exon_stop = int(next(stop for stop in reversed(exon_stops) if variant > int(stop)))
		upstream_exon_stop_index = exon_stops.index(upstream_exon_stop)
		exon_number = upstream_exon_stop_index + 1
		
		if (variant + variant_effect) < upstream_exon_stop:
			sequence_effect = 'upstream effect, variant excluded, sequence deletion'
			new_donor_gain_upstream_exon_stop_distance = upstream_exon_stop - (variant + variant_effect)
			deletion_position_start = exon_cumulative_lengths[upstream_exon_stop_index] - new_donor_gain_upstream_exon_stop_distance
			deletion_position_end = exon_cumulative_lengths[upstream_exon_stop_index]
			novel_transcript = canonical_transcript[:deletion_position_start] + canonical_transcript[deletion_position_end:]
			
		elif upstream_exon_stop <= (variant + variant_effect) <= variant:
			sequence_effect = 'upstream effect, variant excluded, sequence addition'
			new_donor_gain_upstream_exon_stop_distance = (variant + variant_effect) -  upstream_exon_stop
			insertion_position = exon_cumulative_lengths[upstream_exon_stop_index]
			novel_sequence = sequence[(upstream_exon_stop + 1):(upstream_exon_stop + 1 + new_donor_gain_upstream_exon_stop_distance)]
			novel_transcript = canonical_transcript[:insertion_position] + novel_sequence + canonical_transcript[insertion_position:]
			
		elif (variant + variant_effect) > variant:
			sequence_effect = 'downstream effect, variant included, sequence addition'
			new_donor_gain_upstream_exon_stop_distance = (variant + variant_effect) -  upstream_exon_stop
			variant_upstream_exon_stop_distance = variant - upstream_exon_stop
			insertion_position = exon_cumulative_lengths[upstream_exon_stop_index]
			novel_sequence_1 = sequence[(upstream_exon_stop + 1):(upstream_exon_stop + 1 + variant_upstream_exon_stop_distance)]
			novel_sequence_2 = sequence[(upstream_exon_stop + 1 + variant_upstream_exon_stop_distance + 1):(upstream_exon_stop + 1 + new_donor_gain_upstream_exon_stop_distance)]
			novel_transcript = canonical_transcript[:insertion_position] + novel_sequence_1 + alternate_allele + novel_sequence_2 + canonical_transcript[insertion_position:]

	return exon_number, exonic_intronic, sequence_effect, novel_transcript

### Donor Loss ###

def donor_loss(variant, variant_effect):

	# exonic variant
	if any(lower <= variant <= upper for (lower, upper) in zip(exon_starts, exon_stops)):
		exonic_intronic = 'exonic'
		
		downstream_exon_stop = int(next(stop for stop in exon_stops if int(stop) >= variant))
		downstream_exon_stop_index = exon_stops.index(downstream_exon_stop)
		exon_number = downstream_exon_stop_index + 1
		
		if (variant + variant_effect) == downstream_exon_stop and variant == downstream_exon_stop:
			sequence_effect = 'intron retained, variant included, effect at donor site'
			novel_transcript = canonical_transcript[:exon_cumulative_lengths[downstream_exon_stop_index] - 1] + alternate_allele + sequence[(exon_stops[downstream_exon_stop_index] + 1):exon_stops[downstream_exon_stop_index + 1] + 1] + canonical_transcript[exon_cumulative_lengths[downstream_exon_stop_index + 1]:]
		
		elif (variant + variant_effect) == downstream_exon_stop and variant != downstream_exon_stop:
			sequence_effect = 'intron retained, variant included'
			variant_position = exon_cumulative_lengths[downstream_exon_stop_index] - (exon_stops[downstream_exon_stop_index] - variant)
			novel_transcript = canonical_transcript[:variant_position] + alternate_allele + canonical_transcript[(variant_position + 1):exon_cumulative_lengths[downstream_exon_stop_index]] + sequence[(exon_stops[downstream_exon_stop_index] + 1):exon_stops[downstream_exon_stop_index + 1] + 1] + canonical_transcript[exon_cumulative_lengths[downstream_exon_stop_index + 1]:]
			
		else:
			sequence_effect = 'intron not retained, variant included'
			variant_position = exon_cumulative_lengths[downstream_exon_stop_index] - (exon_stops[downstream_exon_stop_index] - variant)
			novel_transcript = canonical_transcript[:variant_position] + alternate_allele + canonical_transcript[(variant_position + 1):exon_cumulative_lengths[downstream_exon_stop_index]] + canonical_transcript[exon_cumulative_lengths[downstream_exon_stop_index]:]

	# intronic variant
	else:
		exonic_intronic = 'intronic'
		
		upstream_exon_stop = int(next(stop for stop in reversed(exon_stops) if variant  >= int(stop)))
		upstream_exon_stop_index = exon_stops.index(upstream_exon_stop)
		exon_number = upstream_exon_stop_index + 1
	
		if (variant + variant_effect) == upstream_exon_stop:
			sequence_effect = 'intron retained, variant included'
			variant_upstream_exon_stop_distance = variant - upstream_exon_stop
			novel_transcript = canonical_transcript[:exon_cumulative_lengths[upstream_exon_stop_index]] + sequence[(upstream_exon_stop + 1):(upstream_exon_stop + variant_upstream_exon_stop_distance)] + alternate_allele + sequence[(upstream_exon_stop + variant_upstream_exon_stop_distance):(exon_starts[upstream_exon_stop_index + 1] - 1)] + canonical_transcript[exon_cumulative_lengths[upstream_exon_stop_index]:]
		
		else:
			novel_transcript = canonical_transcript
			sequence_effect = 'intron not retained, variant excluded'
			
	return exon_number, exonic_intronic, sequence_effect, novel_transcript

def compare(canonical, novel):
	global seq_diff
	global canonical_start_codon
	global novel_start_codon
	global n_premature_stops
	global result	
	
	canonical_codons = [canonical[i:i+3] for i in range(relative_codon_start, len(canonical), 3)]
	novel_codons = [novel[i:i+3] for i in range(relative_codon_start, len(novel), 3)]
	for index, codon in enumerate(canonical_codons):
				if codon == 'TAG' or codon == 'TAA' or codon == 'TGA':
					global first_canonical_stop_codon_index
					first_canonical_stop_codon_index = index
					print(f'The first canonical stop codon index is {first_canonical_stop_codon_index}.')
					break
				else:
					continue
	
	if len(canonical_codons) == 0 or len(novel_codons) == 0:
		seq_diff = len(canonical) - len(novel)
		n_amino_acids = '.'
		canonical_start_codon = '.'
		novel_start_codon = '.'
		n_premature_stops = '.'
		result = 'Translated product unknown'
	
	elif len(canonical) == len(novel):
		seq_diff = 0
		n_amino_acids = first_canonical_stop_codon_index
		canonical_start_codon = canonical_codons[0]
		novel_start_codon = novel_codons[0]
		n_premature_stops = '.'
		
		# Novel transcript same length, same protein
		if canonical_codons[:first_canonical_stop_codon_index] == novel_codons[:first_canonical_stop_codon_index]:
			result = 'No effect'
		
		# Novel transcript same length, different protein
		elif canonical_codons[:first_canonical_stop_codon_index] != novel_codons[:first_canonical_stop_codon_index]:
			result = 'Single missense'
					
	elif len(novel) < len(canonical):
		seq_diff = len(canonical) - len(novel)
		n_amino_acids = first_canonical_stop_codon_index
		canonical_start_codon = canonical_codons[0]	
		novel_start_codon = novel_codons[0]
		n_premature_stops = '.'
		
		# CDS difference
		if canonical[:relative_codon_start] == novel[:relative_codon_start] and canonical_start_codon == 'ATG' and novel_start_codon == 'ATG' and canonical_codons[:first_canonical_stop_codon_index] != novel_codons[:first_canonical_stop_codon_index]:
			result = 'Truncated protein'
		
		# 5' or 3' UTR difference	
		elif canonical[:relative_codon_start] != novel[:relative_codon_start] or canonical[:relative_codon_start] == novel[:relative_codon_start] and canonical_codons[:first_canonical_stop_codon_index] == novel_codons[:first_canonical_stop_codon_index]:
			result = 'Truncated UTR'
				
	elif len(novel) > len(canonical):	
		print(f'The length of canonical codons is {len(canonical_codons)}.')
		print(f'The length of novel codons is {len(novel_codons)}.')
		print(canonical_codons)
		print(novel_codons)

		seq_diff = len(novel) - len(canonical)
		canonical_start_codon = canonical_codons[0]
		novel_start_codon = novel_codons[0]
		n_amino_acids = first_canonical_stop_codon_index
		n_premature_stops = '.'
		
		insertion_start_index = next( ((idx) for idx, (x, y) in enumerate(zip(canonical_codons, novel_codons)) if x!=y), None )
		
		if insertion_start_index is None: # this should be the same as if the diff is in the 3' UTR
			result = 'Longer UTR'
		
		# 5' UTR difference
		elif canonical[:relative_codon_start] != novel[:relative_codon_start]:
			result = 'Longer UTR'
			
		# 3' UTR difference
		elif first_canonical_stop_codon_index <= insertion_start_index:
			result = 'Longer UTR'
		
		# CDS difference
		elif first_canonical_stop_codon_index > insertion_start_index:
			print(f'Insertion start index is {insertion_start_index}.')

			codon_diff = len(novel_codons) - len(canonical_codons)

			insertion_end_index = insertion_start_index + 1 + codon_diff
			print(f'Insertion end index is {insertion_end_index}.')
				
			insertion_codons = novel_codons[insertion_start_index + 1:insertion_end_index]
			print(insertion_codons)

			n_premature_stop_codons = 0

			for codon in insertion_codons:
				if codon == 'TAG' or codon == 'TAA' or codon == 'TGA':
					n_premature_stop_codons += 1
				else:
					continue
			
			print(f'The number of PTCs is {n_premature_stop_codons}.')
			
			if n_premature_stop_codons == 0:
				result = 'Longer protein w/o PTCs'
				n_premature_stops = 0
			elif n_premature_stop_codons > 0:
				result = 'Longer protein w/ PTCs'
				n_premature_stops = n_premature_stop_codons
	
	elif len(canonical) != len(novel) and relative_codon_start == 'not annotated':
		seq_diff = len(novel) - len(canonical)
		n_amino_acids = first_canonical_stop_codon_index
		canonical_start_codon = canonical_codons[0]	
		novel_start_codon = novel_codons[0]
		result = 'Novel transcript is longer but cannot determine start codon'

	return seq_diff, n_amino_acids, canonical_start_codon, novel_start_codon, n_premature_stops, result
	

if __name__ == '__main__':
    main()