#!/bin/bash
#$ -N chagyrskaya_empiric_counts
#$ -t 1-5
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/6_enrichment/chagyrskaya_empiric_counts.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/6_enrichment/chagyrskaya_empiric_counts.err
#$ -l h_rt=24:00:00
#$ -l mem_free=10G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basics

# change directories
cd ../../data/gene_enrichment/observed

# run enrichment
python3 ../../../scripts/6_enrichment/chagyrskaya_empiric_counts.py --set_name chagyrskaya_2 --iterations 2000 --array_id $SGE_TASK_ID
echo "complete"
