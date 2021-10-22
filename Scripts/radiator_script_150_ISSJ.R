library(radiator)
library(SeqArray)
library(tidyverse)

##filter testing radiator package 
#Strata file is similar to the pop map fole for stacks. 
#Is a tab deliminated file but requires the collumns "INDIVIDUALS" and	"STRATA"
#make sure the strata file looks correct. We are expecting 150 individuals and one population
radiator::summary_strata("C:/Users/Rebecca/Desktop/r_coding/data/filtered_strata.issj.tsv")

#Radiator only works if run from a local directory
#need to adjust the parallel.core setting depending on your system otherwise the package will error out. use parallel::detectCores() to see how many available cores you have, and subtract however many to have 1.

#Run radiator using the raw vcf file, and take the blacklist of individuals based on misingness and duplicate genomes and remove them from the popmap to filter them out in popualtions. using this command populations -P ./ -M ./popmap_124_indiv -O ./populations_124_indiv
issj_radiator <- radiator::filter_rad(
  data = "C:/Users/Rebecca/Desktop/r_coding/data/data/populations.snps.names.vcf",
  strata = "C:/Users/Rebecca/Desktop/r_coding/data/data/strata.issj.tsv", 
  interactive.filter=T,
  output="vcf",
  parallel.core = parallel::detectCores() - 3)


#filter strata file to just include the whitelist of individuals
strata <- read_tsv("C:/Users/Rebecca/Desktop/r_coding/data/strata.issj.tsv")


issj_filtered <- read.csv("G:/My Drive/ISSJ_Genomics_150_RAD/ISSJ_Genomics_150_RAD/data/issj_filtered.csv") %>% 
  dplyr::rename("INDIVIDUALS"="tissue_number") 

filtered_strata <- anti_join(strata, issj_filtered,  by="INDIVIDUALS")%>% 
  dplyr::select(!"TARGET_ID")

write_tsv(filtered_strata, "C:/Users/Rebecca/Desktop/r_coding/data/filtered_strata.issj.tsv")


#can then re-run the dataset of filtered individuals and then focus on filtering genotypes. Look at ISSJ_radiator-screen-output_123indiv in the log files directory for the 
#specific parameters and resulting number of SNPs
#note this only works if:
# you have a fresh R session with only radiator and SNPArray packages loaded
#PCs have issues with parallel processing so just have 1 or 1L selected for parallel.core
#Only works if you run in a seperate r script (not markdown) and is run from a local directory
#

issj_radiator <- radiator::filter_rad(
  data = "C:/Users/Rebecca/Desktop/r_coding/data/populations.123indiv.snps.vcf",
  strata = "C:/Users/Rebecca/Desktop/r_coding/data/filtered_strata.issj.tsv", 
  interactive.filter=TRUE,
  output="vcf",
  parallel.core = 1)


