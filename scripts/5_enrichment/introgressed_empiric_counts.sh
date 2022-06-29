#!/bin/bash
#$ -N introgressed_empiric_counts
#$ -t 1-5
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/introgressed_empiric_counts.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/introgressed_empiric_counts.err
#$ -l h_rt=24:00:00
#$ -l mem_free=10G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basics

# change directories
cd ../../data/gene_enrichment/observed

# run enrichment
python3 ../../../scripts/5_enrichment/introgressed_empiric_counts.py --set_name introgressed_2 --iterations 2000 --array_id $SGE_TASK_ID
echo "complete"

# join files
cd ../empiric_counts
join introgressed_2_GWAS_empiric_counts_1.txt introgressed_2_GWAS_empiric_counts_2.txt > introgressed.tmp
for f in introgressed_2_GWAS_empiric_counts_3.txt introgressed_2_GWAS_empiric_counts_4.txt introgressed_2_GWAS_empiric_counts_5.txt; do join introgressed.tmp $f > introgressed.tmpf && mv introgressed.tmpf introgressed.tmp; done
mv introgressed.tmp introgressed_2_GWAS_empiric_counts.txt

join introgressed_2_HPO_empiric_counts_1.txt introgressed_2_HPO_empiric_counts_2.txt > introgressed.tmp
for f in introgressed_2_HPO_empiric_counts_3.txt introgressed_2_HPO_empiric_counts_4.txt introgressed_2_HPO_empiric_counts_5.txt; do join introgressed.tmp $f > introgressed.tmpf && mv introgressed.tmpf introgressed.tmp; done
mv introgressed.tmp introgressed_2_HPO_empiric_counts.txt
