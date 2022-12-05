#!/bin/bash
#$ -N thousand_genomes_spliceosome_variants
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/isoform_builder.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/other_scripts/isoform_builder.err
#$ -l h_rt=8:00:00
#$ -l mem_free=20G

# set path to hg19 change directories
hg19="../../../../data/hg19_fasta/2022-03-14/hg19.fa"
cd ../../data/annotations

# run
python isoform_builder.py --chromosome_lengths hg19.genome --exons hg19_exons.bed --fasta "$hg19" --start_codons hg19_start_codons.txt --variants ../novel_isoforms/isoform_builder_input.txt --out ../novel_isoforms/novel_isoforms.txt
