#!/bin/bash
#$ -N gene_set_to_observed
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/gene_set_to_observed.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/gene_set_to_observed.err
#$ -l h_rt=1:00:00
#$ -l mem_free=2G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basics

# change directories
cd ../../data/gene_enrichment/observed

# run enrichment
python3 ../../../scripts/5_enrichment/gene_set_to_observed.py --set_path ../../genes/archaic_specific_altai_genes_2.txt --set_name altai_2
python3 ../../../scripts/5_enrichment/gene_set_to_observed.py --set_path ../../genes/archaic_specific_chagyrskaya_genes_2.txt --set_name chagyrskaya_2
python3 ../../../scripts/5_enrichment/gene_set_to_observed.py --set_path ../../genes/archaic_specific_denisovan_genes_2.txt --set_name denisovan_2
python3 ../../../scripts/5_enrichment/gene_set_to_observed.py --set_path ../../genes/archaic_specific_vindija_genes_2.txt --set_name vindija_2
python3 ../../../scripts/5_enrichment/gene_set_to_observed.py --set_path ../../genes/archaic_specific_neanderthal_genes_2.txt --set_name neanderthal_2
python3 ../../../scripts/5_enrichment/gene_set_to_observed.py --set_path ../../genes/archaic_specific_shared_genes_2.txt --set_name shared_2
python3 ../../../scripts/5_enrichment/gene_set_to_observed.py --set_path ../../genes/introgressed_genes_2.txt --set_name introgressed_2
python3 ../../../scripts/5_enrichment/gene_set_to_observed.py --set_path ../../genes/archaic_specific_genes_2.txt --set_name archaic_specific_2
echo "complete"
