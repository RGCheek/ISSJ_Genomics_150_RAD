# Habitat-linked genetic variation supports microgeographic adaptive divergence in Island Scrub-jays

Code to study the genetic underpinnings of trait divergence at a microgeographic scale. Manuscript published in Molecular Ecology (DOI: 10.1111/mec.16438). 

These data are also available through Dryad (https://doi.org/10.5061/dryad.8sf7m0cpq), and are liscenced under terms of the Creative Commons Attribution-NonCommercial
License, which permits use, distribution and reproduction in any medium, provided the original work is properly cited and is not used for commercial purposes. 

## Data
### Genomic Data

Neutral and adaptive SNP VCFs also available on Dryad. doi:10.5061/dryad.8sf7m0cpq
ISSJ.ZF.ordered_imputed_BEAGLE -- conveniently in different formats for scripts to run

ISSJ.ZF.ordered_imputed_BEAGLE_neutral -- For PCA 

radiator_data_neutral --For mantel and neutral popstats since you don't need complete data for these to run


### Phenotype & Environment data 

ISSJ_resistance_data.csv -- For MLPE

issj_filtered.csv -- Individuals that were filtered due to missing data. Need this to filter the environmental data 

ISSJ_individual_data.csv -- Phenotype data 

ISSJ_environmental_data_final.csv -- For RDA

## Scripts

001-bioinformatics_scripts_final

radiator_script_150_ISSJ

100-Analyses_Neutral_Pop_Gen_final

200–Identification_of_Population_Structure_Associated_with_Habitat_final

300-Identification_of_loci_Underlying_Variation_in_Morphology_final_bedops 
  Reorder_VCF_issj
  process_satsuma_assembledchroms

## log files
### GWAS results 
GEA_gene_ids_25kb.bed
gemma_gene_ids_25kb.bed
rda_gene_ids_25kb.bed

rda_goids_ensembl.csv
GEA_goids_ensembl.csv
gemma_goids_ensembl.csv

### Gemma
ISSJ.ZF.ordered_imputed_BEAGLE.cXX.txt
ISSJ_centered_relatedness.cXX.txt
ISSJ_standardized_relatedness.sXX.txt

ISSJ.ZFcorr.csv

gene_bill_morph_GO_function.csv


