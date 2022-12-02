This directory contains the inputs, outputs, and scripts that assess enrichment for splice altering variants in various human phenotypes among different gene sets. This approach emulates that of [McArthur et al. 2022](https://www.biorxiv.org/content/10.1101/2022.02.07.479462v1.full). We considered nine different gene sets using a splice altering probability threshold of 0.2: 1) archaic-specific Altai-specific variants, 2) archaic-specific Chagyrskaya-specific variants, 3) archaic-specific Denisovan-specific variants, 4) archaic-specific Vindija-specific variants, 5) archaic-specific Neanderthal-specific variants, 6) archaic-specific shared variants (i.e., those in all four archaics), 7) archaic-specific variants, 8) the Browning et al. 2018 introgressed variants, and 9) the Vernot et al. 2016 introgressed variants.

Why not use existing gene ontology methods? Genes vary in their 'susceptability' to mutations, particularly those that could affect alternative splicing because genes vary in their physical characteristics (e.g., length, number of exons, etc). Relatedly, using all (or many) genes in the genome underpowers enrichment analyses, especially when some genes may not even be represented by analyzed variants.

The pipeline and associated scripts conduct an enrichment analyis using an empirical null distribution. This distribution is created by shuffling the splice altering probability values among all archaic variants, subsetting the relevent gene set (e.g., variant specific to the Altai neanderthal), and then recording how many genes involved in a particular human phenotype from the 2019 GWAS Catolog or Human Phenotype Ontology are observed. Enrichment can then be calculated using the mean number of genes from the empirical distribution compared to the actual number of genes observed in the real data. P values reflect how often the empirical values were more extreme than the observed. We can also calculate an adjusted p value threshold that controls the false discovery rate using a subset of the empirical null.

This directory contains the following subdirectories and files:

- data contains the gene lists for each of the nine sets. These files are same as those in ../genes/ for delta >= 0.2. The "archaic_data_for_enrichment.txt" file was too large to include here but is generated via the analysis notebook.

- empiric_counts contains files with empirically generated null distributions for each term with >= 1 observation per ontology per gene set.

- empiric_FDR contains files with a subset (10%) of the empirically generated null distributions for each term with >= 1 observation per ontology per gene set used to calculate a FDR-corrected p-value threshold.

- observed contains files indicating which genes per term per ontology are present in each of the nine gene sets.

- ontologies contains ontology data to link genes to phenotypes for the two ontologies (GWAS and HPO).

- scripts contains scripts used in the pipeline.

- phenotype_enrichment_pipeline.sh is a bash script that launches the pipeline.

- Snakefile is a Snakemake pipeline that runs the entire analysis.

To rerun this analysis, you will need a Conda environment with snakemake (named "phenotype_enrichment here, see .sh script). Snakemake determines what jobs need to be executed based on file absence/presence. Therefore, you would remove all directories except for data, ontologies, and scripts and then run the phenotype_enrichment_pipeline.sh script. This will result in a replicated analysis. Note that due to the nature of empirical distributions and data shuffling, any reanalysis will be similar but not identical to the results here. 
