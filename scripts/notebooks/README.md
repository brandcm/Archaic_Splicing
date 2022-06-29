This directory contains Jupyter notebooks that directly generate dataframes for analysis or launch scripts that assemble/prepare data. If these analyses are reproduced note that order matters here and is reflected in the numbers in the filenames. 

- 1_generate_archaic_dataframe.ipynb generates the initial dataframe for the analysis.

- 2_get_ancestral_alleles_frequencies_introgressed_variants.ipynb generates additional information for the variants included in the analysis including designating the ancestral allele, gathering 1KG allele frequencies, identifying introgressed variants, and including sQTL data from GTEx.

- 3_generate_archaic_dataframe_with_ancestral_alleles_frequencies_introgressed_variants.ipynb generates a new dataframe with the additional data from the previous notebook except for sQTLs which is handled by the next notebook.

- 4_sQTLs_to_hg19.ipynb maps sQTLs from hg38 to hg19 so that those data can be added to the main dataframe.

- 5_analysis.ipynb generates the final, complete dataframe and conducts all data analysis presented in the manuscript.

- 6_plots.ipynb produces all the plots presented in the manuscript. 
