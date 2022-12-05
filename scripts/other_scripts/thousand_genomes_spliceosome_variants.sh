#!/bin/bash
#$ -N thousand_genomes_spliceosome_variants
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/4_modern_data_preparation/thousand_genomes_spliceosome_variants.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/4_modern_data_preparation/thousand_genomes_spliceosome_variants.err
#$ -l h_rt=24:00:00
#$ -l mem_free=60G

# load modules
module load CBI
module load bcftools/1.13

# change directories
cd ../../data/thousand_genomes/spliceosome_variants

# set variable
inds=('HG00137' 'HG00338' 'HG00619' 'HG01198' 'HG01281' 'HG01524' 'HG01802' 'HG02142' 'HG02345' 'HG02629' 'HG03060' 'HG03190' 'HG03708' 'HG03711' 'HG03800' 'HG04014' 'NA11830' 'NA18552' 'NA18868' 'NA19011' 'NA19452' 'NA19741' 'NA20537' 'NA21141' )

# get variants that fall within major spliceosome complex genes
for ind in ${inds[@]}; do bcftools view -R ../../annotations/major_spliceosome_genes_hg38.bed ../thousand_genomes_avg_SAVs/"$ind"_1KG_variants.vcf.gz -O z -o "$ind"_spliceosome_variants.vcf.gz; done
for ind in ${inds[@]}; do bcftools index -t "$ind"_spliceosome_variants.vcf.gz; done

# query VCFs to VEP format
for ind in ${inds[@]}; do bcftools query -f '%CHROM\t%POS\t'.'\t%REF\t%ALT\t'.'\t'.'\t'.'\n' "$ind"_spliceosome_variants.vcf.gz > "$ind"_spliceosome_variants.txt; done
echo "query complete"
