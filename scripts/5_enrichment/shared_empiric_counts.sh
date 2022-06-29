#!/bin/bash
#$ -N shared_empiric_counts
#$ -t 3
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/shared_empiric_counts.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/5_enrichment/shared_empiric_counts.err
#$ -l h_rt=30:00:00
#$ -l mem_free=10G

# load modules


# load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate basics

# change directories
cd ../../data/gene_enrichment/observed

# run enrichment
python3 ../../../scripts/5_enrichment/shared_empiric_counts.py --set_name shared_2 --iterations 2000 --array_id $SGE_TASK_ID
echo "complete"

# join files
cd ../empiric_counts
join shared_2_GWAS_empiric_counts_1.txt shared_2_GWAS_empiric_counts_2.txt > shared.tmp
for f in shared_2_GWAS_empiric_counts_3.txt shared_2_GWAS_empiric_counts_4.txt shared_2_GWAS_empiric_counts_5.txt; do join shared.tmp $f > shared.tmpf && mv shared.tmpf shared.tmp; done
mv shared.tmp shared_2_GWAS_empiric_counts.txt

join shared_2_HPO_empiric_counts_1.txt shared_2_HPO_empiric_counts_2.txt > shared.tmp
for f in shared_2_HPO_empiric_counts_3.txt shared_2_HPO_empiric_counts_4.txt shared_2_HPO_empiric_counts_5.txt; do join shared.tmp $f > shared.tmpf && mv shared.tmpf shared.tmp; done
mv shared.tmp shared_2_HPO_empiric_counts.txt