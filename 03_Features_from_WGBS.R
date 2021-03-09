## Generating features containing methylation information
## Pre-processed WGBS data is used for this step

library(GenomicRanges)

load("~/MONvsMAC/Meth.gr.RData")
load("~/R Analysis/MethData.filtered.RData")
load("/home/veronikap/MONvsMAC/introns.ret.samp.RData")

## creating GRanges for 5' splice site
retint.5splicepos.gr<-GRangesList()
for (i in 1:length(introns.ret.samp)) {
  retint.5splicepos.gr[[i]]<-makeGRangesFromDataFrame(introns.ret.samp[[i]], seqnames.field = "Chromosome", start.field = "Start", end.field = "Start", strand.field = "Strand")[introns.ret.samp[[i]]$Strand == "+"]
}
retint.5spliceneg.gr<-GRangesList()
for (i in 1:length(introns.ret.samp)) {
  retint.5spliceneg.gr[[i]]<-makeGRangesFromDataFrame(introns.ret.samp[[i]], seqnames.field = "Chromosome", start.field = "End", end.field = "End", strand.field = "Strand")[introns.ret.samp[[i]]$Strand == "-"]
}

retint.5splice.gr<-GRangesList()
for (i in 1:length(introns.ret.samp)) {
  retint.5splice.gr[[i]]<-sort(c(retint.5splicepos.gr[[i]], retint.5spliceneg.gr[[i]]), ignore.strand = TRUE)
}

#Expanding window around the splice site#
expandRange = function(x, upstream=2000, downstream=1000) {
  strand_is_minus = strand(x) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(x)[on_plus] = start(x)[on_plus] - upstream
  start(x)[on_minus] = start(x)[on_minus] - downstream
  end(x)[on_plus] = end(x)[on_plus] + downstream
  end(x)[on_minus] = end(x)[on_minus] + upstream
  x
}


## creating +/- 100 bp window around splice sites
retint.5ss.gr<-GRangesList()
for (i in 1:length(Meth.gr)) {
  retint.5ss.gr[[i]]<-expandRange(retint.5splice.gr[[i]], 100, 99)
}

## finding overlaps with methylation data
retint.5ss.olaps<-list()
for (i in 1:length(Meth.gr)) {
  retint.5ss.olaps[[i]]<-findOverlaps(Meth.gr[[i]], retint.5ss.gr[[i]], ignore.strand = T) #finding overlaping sites between methylation data and windows of interest, specify required Meth.gr object#
}

retint.5ss.meth<-list()
for (i in 1:length(Meth.gr)) {
  retint.5ss.meth[[i]]<-MethData.filtered[[i]][queryHits(retint.5ss.olaps[[i]]),] #filtering methylation dataset with sites that overlap with area of interest#
  retint.5ss.meth[[i]]$Group<-subjectHits(retint.5ss.olaps[[i]]) ##adding new column identifying the matches in windows##
}

retint.5ss.Cs<-list() #counting the number of CpG sites
retint.5ss.methCs<-list() #counting the number of methylated CpG sites
retint.5ss.avemeth<-list() #calculating avergae methylation value
for (i in 1:length(Meth.gr)) {
  retint.5ss.Cs[[i]]<-aggregate(Meth ~ Group, data = retint.5ss.meth[[i]], FUN = length)
  retint.5ss.methCs[[i]]<-aggregate(Meth ~ Group, data = retint.5ss.meth[[i]], FUN = length, subset = Meth.Ratio > 80)
  retint.5ss.avemeth[[i]]<-aggregate(Meth.Ratio ~ Group, data = retint.5ss.meth[[i]], FUN = mean)
}

##adding methylation features
for (i in 1:length(introns.ret.samp)) {
  introns.ret.samp[[i]]$CpGint_100bp5ss<-retint.5ss.Cs[[i]]$Meth[match(introns.ret.samp[[i]]$ID, retint.5ss.Cs[[i]]$Group)]
  introns.ret.samp[[i]]$methCpGint_100bp5ss<-retint.5ss.methCs[[i]]$Meth[match(introns.ret.samp[[i]]$ID, retint.5ss.methCs[[i]]$Group)]
  introns.ret.samp[[i]]$avemethint_100bp5ss<-retint.5ss.avemeth[[i]]$Meth[match(introns.ret.samp[[i]]$ID, retint.5ss.avemeth[[i]]$Group)]
  introns.ret.samp[[i]]$CpGint_100bp5ss<-ifelse(is.na(introns.ret.samp[[i]]$CpGint_100bp5ss), 0, introns.ret.samp[[i]]$CpGint_100bp5ss)
  introns.ret.samp[[i]]$methCpGint_100bp5ss<-ifelse(is.na(introns.ret.samp[[i]]$methCpGint_100bp5ss), 0, introns.ret.samp[[i]]$methCpGint_100bp5ss)
  introns.ret.samp[[i]]$avemethint_100bp5ss<-ifelse(is.na(introns.ret.samp[[i]]$avemethint_100bp5ss), 0, introns.ret.samp[[i]]$avemethint_100bp5ss)
}

## generating GRanges for 3' splice sites
retint.3splicepos.gr<-GRangesList()
for (i in 1:length(introns.ret.samp)) {
  retint.3splicepos.gr[[i]]<-makeGRangesFromDataFrame(introns.ret.samp[[i]], seqnames.field = "Chromosome", start.field = "End", end.field = "End", strand.field = "Strand")[introns.ret.samp[[i]]$Strand == "+"]
}
retint.3spliceneg.gr<-GRangesList()
for (i in 1:length(introns.ret.samp)) {
  retint.3spliceneg.gr[[i]]<-makeGRangesFromDataFrame(introns.ret.samp[[i]], seqnames.field = "Chromosome", start.field = "Start", end.field = "Start", strand.field = "Strand")[introns.ret.samp[[i]]$Strand == "-"]
}

retint.3splice.gr<-GRangesList()
for (i in 1:length(introns.ret.samp)) {
  retint.3splice.gr[[i]]<-sort(c(retint.3splicepos.gr[[i]], retint.3spliceneg.gr[[i]]), ignore.strand = TRUE)
}

## generating GRanges for the middle of intron
retint.middle.gr<-GRangesList()
for (i in 1:length(Meth.gr)) {
  retint.middle.gr[[i]]<-makeGRangesFromDataFrame(introns.ret.samp[[i]], seqnames.field = "Chromosome", start.field = "Middle", end.field = "Middle", strand.field = "Strand")
}

save(introns.ret.samp, file = "/home/veronikap/MONvsMAC/introns.ret.samp.RData")

## methylation data for 3' splice site and middle of intron is summarised in the same way
## the same procedure is repeated for the non-retained introns

##adding CpG density information
library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg38)


retint.5ss.CpGden<-list()
for (i in 1:length(introns.ret.samp)) {
  retint.5ss.CpGden[[i]]<-cpgDensityCalc(retint.5ss.gr[[i]], organism = Hsapiens)
}

for (i in 1:length(introns.ret.samp)) {
  introns.ret.samp[[i]]$CpGsRepi_100bp5ss<-retint.5ss.CpGden[[i]]
}

##same is repeated for 3' splice site and the middle of intron
##and non-retained introns
