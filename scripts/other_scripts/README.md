This directory contains other scripts used in data analysis.

- all_archaic_sites.sh queries the archaic genomes for all variants and subsets sites where all archaics are fixed for the human reference allele. Run this script first before format_fixed_refs.sh, spliceai_fixed_refs.sh, and concat_fixed_refs.sh. Note this script will take quite a long time and eat a lot of memory.

- concat_*_refs.sh concatenates the SpliceAI outputs for the fixed and introgressed reference allele sets. 

- count_1KG_spliceosome_missense_variants.sh takes the outputs from Ensembl's VEP, subsets the missense variants, and counts the number of 1) variants, 2) transcripts affected, 3) damaging variants per PolyPhen, 4) deleterious variants per SIFT, and 5) damaging and deleterious variants per output.

- count_archaic_spliceosome_missense_variants.sh does the same as the above script for each archaic genome.

- format_*_refs.sh formats the data from the fixed and introgressed reference allele sets for input to SpliceAI.

- isoform_builder.py assesses the effect of a SAV on the resulting transcript/protein. It requires six arguments: 1) a chromosome lengths file formatted as a tab-delimited two column file with the chromosome name (e.g., "chr1") in the first field and length in the second field, 2) an exon coordinates file in BED format, 3) a FASTA file for the reference genome, 4) a file with 1-based translation start coordinates from a GFF, 5) a tab-delimited variants file containing the chromosome, 1-based variant position, reference allele, alternate allele, gene name, and the eight SpliceAI outputs, and 6) an output file name.

- isoform_builder.sh launches the isoform_builder.py script with the required input files.

- retrieve_start_codons.sh retrieves the 1-based coordinates for 1 start codon per gene. For genes with multiple ORFs, I prioritized those that were experimentally validated (Ensembl rather than HAVANA) and the longest.

- spliceai_*_refs.sh runs SpliceAI on the fixed and introgressed reference allele sets.
