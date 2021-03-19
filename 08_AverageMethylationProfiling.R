## creating a matrix representing DNA methylation profiles around splice junctions of retained introns

library(GenomicRanges)
library(EnrichedHeatmap)
library(scales)
library(dplyr)
library(reshape2)
library(ggplot2)

load("~/MONvsMAC/Meth.gr.RData")
load("~/R Analysis/MethData.filtered.RData")
load("/home/veronikap/MONvsMAC/introns.ret.samp.RData")
load("/home/veronikap/MONvsMAC/introns.nonret.samp.RData")

## creating GRanges for 5' splice sites of introns

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

##Expanding window around the splice site, in this case +/- 100 bp
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

retint.5ss.gr<-GRangesList()
for (i in 1:length(introns.ret.samp)) {
  retint.5ss.gr[[i]]<-expandRange(retint.5splice.gr[[i]], 100, 99)
}

##cutting the region into the number of required segments with k parameter
##resize is used to make window overlap
##if non-overlaping windows are required, resize step is omitted
retint.targetwin.5ss.gr<-GRangesList()
for (i in 1:length(introns.ret.samp)) {
  retint.targetwin.5ss.gr[[i]]<-makeWindows(retint.5ss.gr[[i]], k = 200)
  retint.targetwin.5ss.gr[[i]]<-resize(retint.targetwin.5ss.gr[[i]], width = 10, fix = "end", use.names = TRUE)
}

retint.5ss.olaps<-list()
for (i in 1:length(introns.ret.samp)) {
  retint.5ss.olaps[[i]]<-findOverlaps(Meth.gr[[i]], retint.targetwin.5ss.gr[[i]], ignore.strand = T) #finding overlaping sites between methylation data and windows of interest, specify required Meth.gr object#
}

retint.5ss.meth<-list()
for (i in 1:length(introns.ret.samp)) {
  retint.5ss.meth[[i]]<-MethData.filtered[[i]][queryHits(retint.5ss.olaps[[i]]),] #filtering methylation dataset with sites that overlap with area of interest#
  retint.5ss.meth[[i]]$Group<-subjectHits(retint.5ss.olaps[[i]]) ##adding new column identifying the matches in windows##
}

retint.5ss.avemeth<-list()
for (i in 1:length(Meth.gr)) {
  retint.5ss.avemeth[[i]]<-aggregate(Meth.Ratio ~ Group, data = retint.5ss.meth[[i]], FUN = mean)
}

background = NA

minusStrand<-list()
for (i in 1:length(introns.ret.samp)) {
  minusStrand[[i]] <- which(as.character(strand(retint.5ss.gr[[i]])) == '-')
}

retint.5ss.mat<-list() #creating matrix for methylation values (beta-values)
for (i in 1:length(introns.ret.samp)) {
  retint.5ss.mat[[i]]<-matrix(retint.5ssvec[[i]], nrow = length(retint.5ss.gr[[i]]), ncol = 200, byrow = T)
  retint.5ss.mat[[i]][minusStrand[[i]],] <- retint.5ss.mat[[i]][minusStrand[[i]], 200:1]
  colnames(retint.5ss.mat[[i]])<-c(1:400)
  rownames(retint.5ss.mat[[i]])<-c(introns.ret.samp[[i]]$IntronID)
}

##same is repeated for non-retained introns
## for 3' ss and the middle of introns

## plot average methylation profiles between retained and non-retained introns
## similar code is used to produce Figure 6 DNA Methylation
mean1<-colMeans(retint.5ss.mat[[1]], na.rm = T)
mean2<-colMeans(nonretint.5ss.mat[[1]], na.rm = T)
averages<-rbind(mean1, mean2)
data<-melt(averages)
ggplot(data, aes(Var2, value, color = Var1)) +
  theme_bw(base_size = 20) +
  geom_line(size = 2) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  labs(y = "DNA methylation") +
  ylim(0,100)+
  theme(axis.title.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 15), legend.background = element_blank()) +
  scale_x_continuous(breaks = c(31, 51, 70), labels = c("-200 bp", "5' ss", "+200 bp")) +
  scale_colour_manual(values = c("#CB2314", "#273046"))+
  geom_vline(xintercept = 31, color = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = 51, color = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = 70, color = "darkgrey", linetype = "dashed") +
  theme(plot.margin=unit(c(0.5,0.7,0.5,0.5),"cm"))
