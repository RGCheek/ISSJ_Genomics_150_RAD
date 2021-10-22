# Script for processing Satsuma Synteny output for Assembled Chromosomes
# This is based on a script written by Benjamin Van Doren 
# 24 Aug 2020                                                                         
# -----------------------------------------------------------

# Set working directory
setwd("G:/My Drive/ISSJ_Genomics_150_RAD/ISSJ_Genomics_150_RAD/data/ZEFI_Annotation/")

# A necessary function
expand.indices <- function(df) {
  if (nrow(df)>1) {
    df.expanded = apply(df,1,function(X) { return(X[1]:X[length(X)]) })
    if (class(df.expanded)=="matrix") {
      df.expanded = lapply(apply(df.expanded,2,list),unlist)
    }
    df.cat = do.call("c",df.expanded)
  } else {
    df.cat = df[1,1]:df[1,2]
  }
  return(df.cat)
}

# Read in original Satsuma output file 
sat <- readLines("./satsuma_summary_all_chr.chained.out")

# The following shortens and breaks up each line of Satsuma Synteny output from this:
# KB442450.1_Taeniopygia_guttata_chromosome_Z_un_1040_1490_LAII01000182.1_Zosterops_lateralis_melanops_scaffold_181_437837_438271_0.793333_-
# to this:
# KB442450.1	chromosome_Z	1040	1490	LAII01000182.1	437837	438271	0.793333	-

#sat <- gsub(",_whole_genome_shotgun_sequence","",sat)
#sat <-gsub("_Taeniopygia_guttata_isolate_Blue55","\t",sat)
#sat <-gsub("_chromosome","Chr",sat)
#sat <- gsub("_Zosterops_lateralis_melanops_scaffold_","\tSca_",sat)

# Now that we're done with this one-time pre-processing, save the file.
# We can now skip the above steps in the future.
fileConn <- file("./PROCESSED_Assembled.txt")
fileConn <- file("./satsuma_summary_all_chr.chained.out_processed")
writeLines(sat, fileConn)
close(fileConn)

# Read in the pre-processed file, not yet as a table and split into columns by the tab character
#chr.map <- readLines("./PROCESSED_Assembled.txt")

library(data.table)
chr.map<-fread("./satsuma_summary_all_chr.chained.out_processed",header = F,sep="\t")
head(chr.map)

# name the columns
#colnames(chr.map) <- c("chr.code","chr","chr.start","chr.end","sca.code","sca.no","sca.start","sca.end","conserve.strength","sca.ori")
colnames(chr.map) <- c("chr.code","chr","chr.start","chr.end","sca.no","sca.start","sca.end","conserve.strength","sca.ori")
#chr.map %>% dplyr::select(sca.ori) %>% distinct()
# Make into a data.frame
chr.map <- as.data.frame(chr.map)
head(chr.map)


# Load as a separate file the query scaffold names and sizes
#this was included in the folder sent by N. Chen
ZLat.scafs <- read.table("G:/My Drive/ISSJ_Genomics_150_RAD/ISSJ_Genomics_150_RAD/data/ZEFI_Annotation/final.assembly2.fasta.gz.ann.txt",header=F) %>%
  dplyr::select(c("V2","V5")) %>% 
  rename(scaffold=V2,length=V5)

head(ZLat.scafs)
chr.map[,1] <- as.character(chr.map[,1])
chr.map[,2] <- as.character(chr.map[,2])
chr.map[,5] <- as.character(chr.map[,5])
#chr.map[,6] <- as.character(chr.map[,6])
for (i in c(3:4,6:8)) {
  chr.map[,i] = as.numeric(as.character(chr.map[,i]))
}

head(chr.map)
chr.map2<-chr.map %>% mutate(sca.ori2=if_else(sca.ori=="+",1,-1)) %>% dplyr::select(-sca.ori) %>% rename(sca.ori=sca.ori2)
chr.map = chr.map2
#specific for Zlat
#chr.map[,1] <- as.character(chr.map[,1])
#chr.map[,2] <- as.character(chr.map[,2])
#chr.map[,5] <- as.character(chr.map[,5])
#chr.map[,6] <- as.character(chr.map[,6])
#for (i in c(3:4,7:10)) {
#  chr.map[,i] = as.numeric(as.character(chr.map[,i]))
#}

chr.map$len.sca.chunk <- chr.map$sca.end-chr.map$sca.start
chr.map$sca.length <- ZLat.scafs$length[match(chr.map$sca.no,ZLat.scafs$scaffold)]
chr.map$prop.of.sca <- chr.map$len.sca.chunk/chr.map$sca.length

head(chr.map)
# "chr.map" should now look something like this:
# chr.code   chr chr.start   chr.end     sca.no sca.start  sca.end conserve.strength sca.ori len.sca.chunk sca.length  prop.of.sca
# 1 CM018230.1 Chr_1 114374234 114374259 scaffold_3  16937604 16937629          1.000000       1            25   25052223 9.979154e-07
# 2 CM018230.1 Chr_1   5378513   5380784 scaffold_4  10143143 10145389          0.806253       1          2246   24224907 9.271449e-05
# 3 CM018230.1 Chr_1   5380823   5381198 scaffold_4  10145415 10145792          0.784000       1           377   24224907 1.556250e-05
# 4 CM018230.1 Chr_1   5381197   5381554 scaffold_4  10145793 10146157          0.761905       1           364   24224907 1.502586e-05
# 5 CM018230.1 Chr_1   5381575   5383833 scaffold_4  10146165 10148380          0.784544       1          2215   24224907 9.143482e-05
# 6 CM018230.1 Chr_1   5383854   5383984 scaffold_4  10148404 10148534          0.580769       1           130   24224907 5.366378e-06

# The last column, "prop.of.sca," tells you how much of the scaffold is mapped by that mapping
# It might be low, but what matters is when you add things up.

# Write to csv
#write.csv(chr.map, "chr.mapAssembled.csv")


#------------------------------
#STEP 3: CREATING SCA.CHROM.MAP

# load in necessary files:
chr.map <- read.csv("G:/My Drive/ISSJ_Genomics_150_RAD/ISSJ_Genomics_150_RAD/data/ZEFI_Annotation/chr.mapAssembled.csv",stringsAsFactors = F)
head(chr.map)

ZLat.scafs <- read.table("final.assembly2.fasta.gz.ann.txt",header=F) %>%
  dplyr::select(c("V2","V5")) %>% 
  rename(scaffold=V2,length=V5)


# necessary function:
expand.indices <- function(df) {
  if (nrow(df)>1) {
    df.expanded = apply(df,1,function(X) { return(X[1]:X[length(X)]) })
    if (class(df.expanded)=="matrix") {
      df.expanded = lapply(apply(df.expanded,2,list),unlist)
    }
    df.cat = do.call("c",df.expanded)
  } else {
    df.cat = df[1,1]:df[1,2]
  }
  return(df.cat)
}


# This loop could take a while (a minute or two or three?) but will speed up as you get to the smaller scaffolds
sca.chrom.map <- data.frame(ZLat.scaffold=ZLat.scafs$scaffold,best.ZFinch.chrom=NA,mean.loc.ZFinch.chrom=NA,ori=NA,prop.of.scaf=NA)


split.scaffolds <- list()
for (i in 1:nrow(ZLat.scafs)) {
  sca = ZLat.scafs$scaffold[i]
  chr.sca = chr.map[chr.map$sca.no==sca,]
  chr.sca = chr.sca[order(chr.sca$sca.start),]
  # get chromosome with most coverage of scaffold
  chr.cov = tapply(chr.sca$prop.of.sca,chr.sca$chr,sum)
  chr.cov = chr.cov[chr.cov>0.2]
  # If more than one chromosome mapped to over 20% of the scaffold's length, save it here for examination
  if (length(chr.cov)>1) { 
    split.scaffolds = c(split.scaffolds,list(rbind(paste(unique(chr.sca$sca)),round(chr.cov,2))))
  } 
  if (length(chr.cov)==0) {
    sca.chrom.map$best.ZFinch.chrom[i] = sca.chrom.map$mean.loc.ZFinch.chrom[i] = sca.chrom.map$ori[i] = sca.chrom.map$prop.of.scaf[i] = NA
  } else { # all good
    chr.max = names(chr.cov[chr.cov==max(chr.cov)])
    # take weighted mean of position by coverage on chromsome
    chr.sca = chr.sca[chr.sca$chr==chr.max,]
    mean.pos = mean(expand.indices(chr.sca[,c("chr.start","chr.end")]))
    stopifnot(sca.chrom.map$ZLat.scaffold[i]==sca)
    sca.chrom.map$best.ZFinch.chrom[i] = chr.max
    sca.chrom.map$mean.loc.ZFinch.chrom[i] = mean.pos
    sca.chrom.map$prop.of.scaf[i] = max(chr.cov)
    ori.mean = weighted.mean(chr.sca$sca.ori,chr.sca$prop.of.sca)
    if (ori.mean>0) { 
      sca.chrom.map$ori[i] = 1
    } else {
      sca.chrom.map$ori[i] = -1
    }
    
  }
  print(nrow(ZLat.scafs)-i)
}
sca.chrom.map <- sca.chrom.map[complete.cases(sca.chrom.map),]
head(sca.chrom.map$best.ZFinch.chrom)
unique(sca.chrom.map$best.ZFinch.chrom)
# The next lines hack together a way of ordering the names of the Zebra Finch chromosomes such that
# ones of "unknown" position ("un") come after the known ones (e.g., chr_1 vs chr_1_un)
# and the Z and LGE22 chromosomes are last
o <- substring(regmatches(sca.chrom.map$best.ZFinch.chrom,regexpr("_.*",sca.chrom.map$best.ZFinch.chrom)),2)
o <- gsub("A",".1",o)
o <- gsub("29_","",o)
o <- gsub("_ctg1","",o)
o <- gsub("pat_bTG3_random","99",o)
o <- gsub("_.*","",o)
o <- as.numeric(o)

sca.chrom.map <- sca.chrom.map[order(o,sca.chrom.map$mean.loc.ZFinch.chrom),]

sca.chrom.map$best.ZFinch.chrom <- gsub(".*pat_bTG3_random_.*","UN",sca.chrom.map$best.ZFinch.chrom)

#------------------------------
#STEP 4: ORDERING SCA.CHROM.MAP & SAVING OUTPUT FILE

# measuring length of sccessfully mapped scaffolds and comparing it to the total length of the included scaffolds
sum(ZLat.scafs$length[match(as.character(sca.chrom.map$ZLat.scaffold),ZLat.scafs$scaffold)])/sum(ZLat.scafs$length)
# For me, this value was 0.9630613 <- proportion mapped to Zebra Finch Genome
head(sca.chrom.map)
# Now, you may want to output the "sca.chrom.map" object for further use. 
# I pared it down further to just the three columns I wanted before saving it: 
scaffold.order <- sca.chrom.map[,c("best.ZFinch.chrom","ZLat.scaffold","mean.loc.ZFinch.chrom","ori")]  
colnames(scaffold.order) <- c("chr","sca","mean.loc","sca.ori")
head(scaffold.order)
#write.csv(scaffold.order,"ISSJ_scaffold_order_from_ZFinch.csv",row.names=F)

