library(tidyverse)
library(data.table)
library(GenomicRanges)
library(Rmisc)

#read in a .bed file with gene names annotated
Heidi.bed <- fread("/mnt/home/stephen/yard/Heidi/MROs/Heidi.gene_names.bed")
#rename the columns
colnames(Heidi.bed) <- c("chrom","start","end", "gene")
#read in a list of samples
dirs <- read.table("/mnt/home/stephen/yard/Heidi/dirs")
#make the annotation file a GRanges
gr2 <- with(Heidi.bed, GRanges(chrom, IRanges(start = start, end = end), gene = gene))

#read in all the BC counting files
x=list()
for (i in 1:nrow(dirs)){
  x[i] <- list(fread(paste("/mnt/home/stephen/yard/Heidi/MROs/MRO_",dirs[i,1],"/outs/bc_count.bed", sep = ""), header = T))
}

#make all the BC counting data.tables a GRanges
gr1=list()
for (i in 1:nrow(dirs)){
  gr1[i] <- list(makeGRangesFromDataFrame(x[i], keep.extra.columns = T))
}

#find the overlaps between the annotation file and the BC counting data.tables
overlap=list()
for (i in 1:nrow(dirs)){
  overlap[i]<- list(findOverlaps(gr1[[i]], gr2))
}

# "annotate" the BC counting GRanges
for (i in 1:nrow(dirs)){
  mcols(gr1[[i]])$gene[queryHits(overlap[[i]])] <- mcols(gr2)$gene[subjectHits(overlap[[i]])]
}


#make the final annotated GRanges a data.table
final=list()
for (i in 1:nrow(dirs)){
  final[i] <- list(as.data.table(gr1[[i]]))
  final[[i]]$gene <- as.factor(final[[i]]$gene)
  final[[i]]$hap1_hap2_ratio <- final[[i]]$bc_hap1/final[[i]]$bc_hap2
  final[[i]] <- subset(final[[i]], bc_hap1>=5 & bc_hap2 >=5)
  
}

#get summary stats for each sample
final_SE=list()
for (i in 1:nrow(dirs)){
  final_SE[i] <- list(summarySE(final[[i]], measurevar = "hap1_hap2_ratio", groupvars = "gene"))
}

#make the plots
plot=list()
for (i in 1:nrow(dirs)){
 plot[i] <-list( 
      ggplot(final[[i]], aes(x=start,y=hap1_hap2_ratio, color=gene))+
        geom_point()+
        geom_hline(data=final_SE[[i]], aes(yintercept=hap1_hap2_ratio+(2*sd)))+
        geom_hline(data=final_SE[[i]], aes(yintercept=hap1_hap2_ratio-(2*sd)))+
        facet_wrap(~gene, ncol = 10, scales = "free")+
        scale_y_continuous(limits = c(0,5))+
        xlab("Position bp")+
        ylab("Hap 1/Hap 2")+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.line.x = element_line(),
              axis.line.y = element_line(),
              axis.title.y = element_text(face="bold", size = rel(2), vjust = 1),
              axis.title.x = element_text(face="bold", size = rel(2), vjust = 1),
              panel.background = element_blank(),
              legend.position="none")
        
      )
}


#save all the plots
for (i in 1:nrow(dirs)){
  ggsave(plot = plot[[i]], filename = paste("/your/dir",dirs[i,1],".pdf", sep = ""),
         height = 20, width = 20)
}
