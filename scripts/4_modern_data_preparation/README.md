This directory houses scripts that generate additional information for both splice altering and non-splice altering variants among four archaic hominins.

- calculate_non_ASW_AFR_allele_frequencies.sh retrieves information needed to recalculate the allele frequency for the 1KG superpopulation without the ASW (African Ancestry in Southwest US). Individuals in this population have higher proportions of non-African admixture than the other African populations so we recalculate this frequency to capture allele frequencies for populations sampled directly from Africa.

- concat_introgressed_variants.sh concatenates introgressed tag SNPs among the ASN, EUR, PNG, and SAS populations from [here](https://drive.google.com/drive/folders/0B9Pc7_zItMCVM05rUmhDc0hkWmc?resourcekey=0-zwKyJGRuooD9bWPRZ0vBzQ). We used the all_tag_snps.*.merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed files. Duplicate entries that occur from their presence in multiple populations are also removed.

- concat_sQTLs.sh concatenates GTEx sQTL data among 49 tissues from [here](https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_sQTL.tar). 

- get_1KG_sample_names.sh retrieves sample names from the 1KG VCF.

- query_1KG_allele_frequencies.sh retrieves allele counts, numbers, and frequencies for archaic variants that are also present in 1KG.

- query_gnomAD_allele_frequencies.sh retrieves allele counts, numbers, and frequencies for archaic variants that are also present in gnomAD.

- subset_1KG.sh subsets variants in 1KG within GENCODE, Human Release 24 for hg38, gene bodies that also occur in the archaic hominins.

- thousand_genomes_spliceosome_variants.sh subsets variants among 147 spliceosome associated genes in 24 randomly sampled 1KG individuals (1 per population).
