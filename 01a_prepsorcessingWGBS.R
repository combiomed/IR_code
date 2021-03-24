## Processing WGBS files
## from Bismark output

library(dplyr)
library(GenomicRanges)

##loading raw output files
pathsMeth<-c("WGBS/Hf03_BlTN.bismark.cov",
             "WGBS/Hf04_BlTN.bismark.cov",
             "WGBS/Hf03_BlCM.bismark.cov",
             "WGBS/Hf04_BlCM.bismark.cov",
             "WGBS/Hf03_BlEM.bismark.cov",
             "WGBS/Hf04_BlEM.bismark.cov")

MethData<-list() #create empty list object##
for (i in 1:length(pathsMeth)) {
  MethData[[i]] <- read.table(pathsMeth[i], header = F)
  MethData[[i]]$Coverage=MethData[[i]]$V5 + MethData[[i]]$V6
  colnames(MethData[[i]])<-c("Chromosome", "Start", "End", "Meth.Ratio", "Meth", "Unmeth", "Coverage")
}

MethData.filtered<-list() ##filtering out low-count reads##
for (i in 1:length(MethData)) {
  MethData.filtered[[i]]<-subset(MethData[[i]], Coverage>5)
}

#Adding strand info to the meth dataset#
pathsMethStrand<-c("WGBS/Hf03_BlTN.CpG.txt",
                   "WGBS/Hf04_BlTN.CpG.txt",
                   "WGBS/Hf03_BlCM.CpG.txt",
                   "WGBS/Hf04_BlCM.CpG.txt",
                   "WGBS/Hf03_BlEM.CpG.txt",
                   "WGBS/Hf04_BlEM.CpG.txt")

MethStrand<-list()
for (i in 1:length(pathsMethStrand)) {
  MethStrand[[i]] <- read.table(pathsMethStrand[i], header = F)
}


for (i in 1:length(MethData.filtered)) {
  MethData.filtered[[i]]$Strand=MethStrand[[i]]$V3[match(paste(MethData.filtered[[i]]$Chromosome,MethData.filtered[[i]]$Start,sep=":"),paste(MethStrand[[i]]$V1,MethStrand[[i]]$V2,sep=":"))]
  MethData.filtered[[i]]$Strand<-as.character(MethData.filtered[[i]]$Strand)
  MethData.filtered[[i]]$Strand[is.na(MethData.filtered[[i]]$Strand)]="*"
}


#Ordering dataset by chromosomal position#
chrOrder<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrMT")
for (i in 1:length(MethData.filtered)) {
  MethData.filtered[[i]]$Chromosome<-sub("^", "chr", MethData.filtered[[i]]$Chromosome)
  MethData.filtered[[i]]$Chromosome<-factor(MethData.filtered[[i]]$Chromosome, levels=chrOrder)
  MethData.filtered[[i]]<-MethData.filtered[[i]][order(MethData.filtered[[i]]$Chromosome),]
  MethData.filtered[[i]]$Diff<-(lead(MethData.filtered[[i]]$Start, 1) - MethData.filtered[[i]]$Start)
}

save(MethData.filtered, file = "MethData.filtered.RData") ##tis files keeps reads from both strands
load("MethData.filtered.samp.RData")

MethData.filtered<-list() ##this files keeps reads from one strand only, to avoid bias in calculating the number of methylated CpG sites
for (i in 1:length(MethData.filtered)){
  MethData.filtered[[i]]<-subset(MethData.filtered[[i]], Diff>1)
  MethData.filtered[[i]]$Chromosome<-factor(MethData.filtered[[i]]$Chromosome, levels=chrOrder)
  MethData.filtered[[i]]<-MethData.filtered.subset[[i]][order(MethData.filtered[[i]]$Chromosome),]
}

Meth.gr<-GRangesList()
for (i in 1:length(MethData.filtered)) {
  Meth.gr[[i]]<-makeGRangesFromDataFrame(MethData.filtered[[i]], seqnames.field = "Chromosome", start.field = "Start", end.field = "End", strand.field = "Strand")
}

save(MethData.filtered, file = "MethData.filtered.RData")
load("MethData.filtered.RData")

save(Meth.gr, file = "Meth.gr.RData")
load("Meth.gr.RData")
