This directory contains scripts used to assemble and analyze archaic hominin genomic data, particulary on alternative splicing. There are several directories grouped by the stage of analysis. Numbers in the directory names reflect the order of analysis.

- [1_archaic_genome_preparation_and_filtering](https://github.com/brandcm/Archaic_Splicing/tree/main/scripts/1_archaic_genome_preparation_and_filtering) contains scripts to assemble archaic InDels and SNVs separately and filter variants to retain high-quality sites and genotypes.

- [2_archaic_data_preparation_for_SpliceAI_analysis](https://github.com/brandcm/Archaic_Splicing/tree/main/scripts/2_archaic_data_preparation_for_SpliceAI_analysis) contains scripts to 1) subset the filtered VCFs above to those sites near/in exons per GENCODE, Human Release 23 for hg19, annotations and 2) split variants into smaller VCFs containing <= 5,000 variants for faster, parallelized SpliceAI analysis.

- [3_SpliceAI_analysis](https://github.com/brandcm/Archaic_Splicing/tree/main/scripts/3_SpliceAI_analysis) contains scripts to run SpliceAI on the split VCFs and reassemble the annotations. Note that in this last step, I split variants with multiple annotations into multiple rows for easier downstream analysis.

- [4_modern_data_preparation](https://github.com/brandcm/Archaic_Splicing/tree/main/scripts/4_modern_data_preparation) contains scripts to prepare additional data that is useful for contextualizing and understanding potential archaic splice variants. Many of these scripts are called in a [Jupyter notebook](https://github.com/brandcm/Archaic_Splicing/blob/main/scripts/notebooks/2_get_ancestral_alleles_frequencies_introgressed_variants.ipynb) in the [notebooks](https://github.com/brandcm/Archaic_Splicing/tree/main/scripts/notebooks) directory.

- [5_enrichment](https://github.com/brandcm/Archaic_Splicing/tree/main/scripts/5_enrichment) contains scripts to assess enrichment for splice altering variants among genes that underlie different human phenotypes from the 2019 GWAS Catalog and Human Phenotype Ontology.

- [notebooks](https://github.com/brandcm/Archaic_Splicing/tree/main/scripts/notebooks) contains Jupyter notebooks used in data analysis and visualization. All results and figures presented in the manuscript can be found in this directory.