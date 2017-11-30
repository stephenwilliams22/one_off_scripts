#This script will find all the necessary files from MARSOC that go into the barcode counting pipeline of LR

#read in finder tools for MARSOC
source("/mnt/home/stephen/crg_paper/R/utils/find_marsoc_files.R")

df <- read.delim("your lena IDs")

for (i in 1:nrow(df)) {
fileConn<-file(paste("/dir/to/put/MROs",df[i,1],".mro", sep = ""))
  writeLines(paste("@include simple_bc_count.mro \n call _CNV_EXOME( \n fragments_h5 = ",'"',find_wgs_fragments.h5(df[i,1]),'"',",", 
                   "\n fragment_phasing =", '"',find_fragment_phase(df[i,1]),'"',",",
                   "\n phased_possorted_bam =",'"',find_bam(df[i,1]),'"',",",
                   '\n queryRegions ="Heidi.bed"',",",
                   "\n bin=50",
                   "\n )", sep = ""), fileConn)
  close(fileConn)
  
}
