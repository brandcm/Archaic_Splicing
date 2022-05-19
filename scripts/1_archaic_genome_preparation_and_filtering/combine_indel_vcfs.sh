#!/bin/bash
#$ -N combine_indel_vcfs
#$ -M colin.brand@ucsf.edu
#$ -m ae
#$ -cwd
#$ -o ~/../../../group/capra/projects/archaic_splicing/scripts/1_archaic_genome_preparation_and_filtering/combine_indel_vcfs.out
#$ -e ~/../../../group/capra/projects/archaic_splicing/scripts/1_archaic_genome_preparation_and_filtering/combine_indel_vcfs.err
#$ -l h_rt=10:00:00
#$ -l scratch=50G

# load modules
module load CBI
module load bcftools/1.13

# copy needed files into the temporary directory
cd "$TMPDIR"
cp ~/../../../group/capra/data/altai_neanderthal/2022-01-12/*new.vcf.gz . &
cp ~/../../../group/capra/data/altai_neanderthal/2022-01-12/*new.vcf.gz.tbi . &
cp ~/../../../group/capra/data/denisovan/2022-01-12/*new.vcf.gz . &
cp ~/../../../group/capra/data/denisovan/2022-01-12/*new.vcf.gz.tbi . &
cp ~/../../../group/capra/data/vindija_neanderthal/2022-01-12/*new.vcf.gz . &
cp ~/../../../group/capra/data/vindija_neanderthal/2022-01-12/*new.vcf.gz.tbi . &
cp ~/../../../group/capra/projects/archaic_splicing/data/archaic_genomes/new_chr_names.txt . &
cp ~/../../../group/capra/projects/archaic_splicing/data/archaic_genomes/new_altai_name.txt . &
cp ~/../../../group/capra/projects/archaic_splicing/data/archaic_genomes/new_vindija_name.txt .
wait
echo "moving files finished"

# rename sample for Altai because the InDel sample name is different to the SNVs sample name
bcftools reheader -s new_altai_name.txt Altai_chrALL_indels_new.vcf.gz -o Altai_2_chrALL_indels_new.vcf.gz

# rename sample for Vindija because the InDel sample name is different to the SNVs sample name
bcftools reheader -s new_vindija_name.txt Vindija33.19_chrALL_indels_new.vcf.gz -o Vindija_2_chrALL_indels_new.vcf.gz
echo "samples renamed"

# index files
bcftools index -t Altai_2_chrALL_indels_new.vcf.gz
bcftools index -t Vindija_2_chrALL_indels_new.vcf.gz

# combine VCFs
bcftools merge Altai_2_chrALL_indels_new.vcf.gz Denisova_chrALL_indels_new.vcf.gz Vindija_2_chrALL_indels_new.vcf.gz -O z -o unsorted_combined_indels.vcf.gz
bcftools sort unsorted_combined_indels.vcf.gz -O z -o combined_indels.vcf.gz
echo "VCFs merged and sorted"

# rename chromosomes and samples, index new VCF
bcftools annotate --rename-chrs new_chr_names.txt combined_indels.vcf.gz -o combined_archaic_indels.vcf.gz
echo "chromosomes renamed"

bcftools index -t combined_archaic_indels.vcf.gz
echo "combined VCF indexed"

# move new VCFs and indices to archaic_splicing directory
mv combined_archaic_indels.vcf.gz ~/../../../group/capra/projects/archaic_splicing/data/archaic_genomes
mv combined_archaic_indels.vcf.gz.tbi ~/../../../group/capra/projects/archaic_splicing/data/archaic_genomes
echo "script completed"

