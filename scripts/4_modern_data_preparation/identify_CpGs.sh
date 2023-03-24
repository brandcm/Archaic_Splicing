#!/bin/bash
#$ -N identify_CpGs
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/4_modern_data_preparation/identify_CpGs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/4_modern_data_preparation/identify_CpGs.err
#$ -l h_rt=24:00:00
#$ -l mem_free=60G

# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ancestral_allele

# change directories
cd ../../data/archaic_variants_in_humans

# set path
hg19="../../../../data/hg19_fasta/2022-03-14/hg19.fa"

# identify CpGs among the gnomAD variants also present in the archaic hominins
python ../../scripts/4_modern_data_preparation/identify_CpGs.py --variants archaic_SNVs_in_gnomAD_allele_frequencies_hg19.txt --fasta "$hg19" --out archaic_SNVs_in_gnomAD_allele_frequencies_hg19_with_CpG_annotation.txt
