#!/bin/bash
#$ -N altai_empiric_counts
#$ -t 1-5
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/altai_empiric_counts.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/altai_empiric_counts.err
#$ -l h_rt=24:00:00
#$ -l mem_free=10G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basics

# change directories
cd ../../data/gene_enrichment/observed

# run enrichment
#python3 ../../../scripts/5_enrichment/altai_empiric_counts.py --set_name altai_2 --iterations 2000 --array_id $SGE_TASK_ID
#echo "complete"

# join files
cd ../empiric_counts
join altai_2_GWAS_empiric_counts_1.txt altai_2_GWAS_empiric_counts_2.txt > altai.tmp
for f in altai_2_GWAS_empiric_counts_3.txt altai_2_GWAS_empiric_counts_4.txt altai_2_GWAS_empiric_counts_5.txt; do join altai.tmp $f > altai.tmpf && mv altai.tmpf altai.tmp; done
mv altai.tmp altai_2_GWAS_empiric_counts.txt

join altai_2_HPO_empiric_counts_1.txt altai_2_HPO_empiric_counts_2.txt > altai.tmp
for f in altai_2_HPO_empiric_counts_3.txt altai_2_HPO_empiric_counts_4.txt altai_2_HPO_empiric_counts_5.txt; do join altai.tmp $f > altai.tmpf && mv altai.tmpf altai.tmp; done
mv altai.tmp altai_2_HPO_empiric_counts.txt
