#24 Aug 2020
# Genomic landscape of divergence in Island Scrub-jays
# R script used to reorder vcf file based on satsuma output
#Modified from script from CH Bossu
# ------------------------------------------------------------

setwd("G:/My Drive/ISSJ_Genomics_150_RAD/ISSJ_Genomics_150_RAD/data/ZEFI_Annotation/")

# STEP 1: Read in and reformat the satsuma and scaffold data


library(tidyverse)
scaff_ord <- read_csv("ISSJ_scaffold_order_from_ZFinch.csv") %>%
  rename(scaffold = sca)


# Load as a separate file the query scaffold names and sizes
  scaff_lengths <- read.table("final.assembly2.fasta.gz.ann.txt",header=F) %>%
  dplyr::select(c("V2","V5")) %>% 
  rename(scaffold=V2,length=V5)

missing_lengths <- anti_join(scaff_ord, scaff_lengths)

missing_ords <- anti_join(scaff_lengths, scaff_ord)


scaffs <- scaff_ord %>%
  left_join(scaff_lengths) %>%
  rename(ZFCHROM = chr, CHROM = scaffold)

scaffs



vcf <- read.table("C:/Users/Rebecca/Desktop/r_coding/data/data/Filtered_Data/radiator_data.vcf", header=F) %>% 
  rename(CHROM=V1,POS=V2) 

combo <- left_join(vcf, scaffs, by="CHROM")



zf_ified <- combo %>%
  mutate(ZFPOS = {
    ML = floor(mean.loc)  # some temp variables to make it easier to express
    Lo2 = floor(length/2)
    L = length
    ifelse(sca.ori == 1,
           ML - Lo2 + POS,           # forward orientation
           ML - Lo2 + (L - POS))     # reverse orientation
  }) %>%
  dplyr::select(ZFCHROM, ZFPOS, everything()) %>%
  mutate(ZFCHROM = factor(ZFCHROM, levels = unique(ZFCHROM))) %>%   # this is to get them to sort correctly
  arrange(ZFCHROM, ZFPOS)

head(zf_ified)

#Write the coordinate information into a file for use after GWAS 
ZEFIcor <- zf_ified[,1:5]
#write.csv(ZEFIcor, file="C:/Users/Rebecca/Desktop/r_coding/data/data/Filtered_Data/ISSJ.ZFcorr.csv", quote=F, row.names=F) 


#now, remove the data from the original alignment so it's jut the ZF positional information left
zf_ified$CHROM <- NULL
zf_ified$POS <- NULL
zf_ified$mean.loc <- NULL
zf_ified$sca.ori <- NULL
zf_ified$length <- NULL

#change the column name and position name to make the file a usable vcf
colnames(zf_ified)[1] <- "CHROM"
colnames(zf_ified)[2] <- "POS"

#Filter so that negative positions and unmapped scaffolds are removed:
zf_ified <- filter(zf_ified, POS >0)
zf_ified <- na.omit(zf_ified)

#zf_ified <- zf_ified[!grepl("_Un", zf_ified$`CHROM`),]
head(zf_ified)
library(pgirmess)

#write.delim(zf_ified, file="C:/Users/Rebecca/Desktop/r_coding/data/data/Filtered_Data/ISSJ.ZF.ordered.txt", quote=F,  sep = "\t")


# Note: To make into a functional vcf file need to add header. This is done by copy pasting the header and first row of the original vCF
#to the top of the file uing Notepad++.(or some other test editor) Ignore the fist #:
#  ###fileformat=VCFv4.3																																																																																																																																			
#  ##fileDate=20200707@2120																																																																																																																																			
#  ##source=radiator_v.1.1.5																																																																																																																																			
#  "##INFO=<ID=NS,Number=1,Type=Integer,Description=""Number of Samples With Data"">"																																																																																																																																			
#  "##FORMAT=<ID=GT,Number=1,Type=String,Description=""Genotype"">"																																																																																																																																			
#  #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	F09-100	F09-105	F09-108	F09-111	F09-112	F09-118	F09-119	F09-124	F09-132	F09-136	F09-138	F09-139	F09-21	F09-23	F09-26	F09-33	F09-35	F09-36	F09-39	F09-40	F09-46	F09-47	F09-48	F09-53	F09-56	F09-64	F09-69	F09-76	F09-82	F09-86	F09-91	F09-93	F09-96	F10-105	F10-106	F10-108	F10-111	F10-12	F10-120	F10-126	F10-13	F10-138	F10-14	F10-149	F10-15	F10-150	F10-152	F10-156	F10-161	F10-163	F10-165	F10-166	F10-167	F10-169	F10-170	F10-172	F10-175	F10-177	F10-22	F10-25	F10-30	F10-31	F10-33	F10-36	F10-39	F10-42	F10-53	F10-54	F10-57	F10-60	F10-62	F10-64	F10-65	F10-66	F10-76	F10-77	F10-78	F10-82	F10-84	F10-86	F10-87	F10-89	F11-003	F11-005	F11-006	F11-007	F11-009	F11-017	F11-018	F11-021	F11-027	F11-050	F11-058	F11-059	F11-060	F11-061	F11-065	F11-067	F11-068	F11-070	F11-071	F11-072	F11-077	F11-081	F11-083	F11-085	F11-087	F11-088	F11-090	F11-095	F11-097	F11-100	F11-102	F11-105	F11-115	F11-118	F11-141	F11-142	F11-146	F11-147	F11-154	F11-164	F11-169
# 

