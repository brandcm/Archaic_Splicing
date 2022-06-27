This repository contains data used in the gene enrichment analysis.

- GWAS.txt.gz is the gene set library for 1,737 GWAS Catalog terms (retrieved from https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GWAS_Catalog_2019). GWAS_terms_and_systems.txt.gz lists the term and manually curated system label.

- HPO.txt.gz is the gene set library for 1,779 Human Phenotype Ontology terms (retrieved from https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Human_Phenotype_Ontology). HPO_terms_and_systems.txt.gz lists the term and manually curated system label.

- phenotype_enrichment.txt.gz summarizes the enrichment results for each term per gene set. This file is concatenated version of individual files that were previously housed in an 'enrichment' directory that is referenced in multiple scripts/notebooks. The directory is not included here for brevity. 

- empiric_counts contains files with empirically generated null distributions for each term with >= 1 observation per ontology per gene set.

- empiric_FDR contains files with a subset (10%) of the empirically generated null distributions for each term with >= 1 observation per ontology per gene set used to calculate a FDR-corrected p-value threshold.

- observed contains files indicating which genes per term per ontology are present in each of the eight gene sets.
