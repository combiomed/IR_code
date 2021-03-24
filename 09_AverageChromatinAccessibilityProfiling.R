## creating a matrix representing DNA methylation profiles around splice junctions of retained introns

library(GenomicRanges)
library(EnrichedHeatmap)
library(scales)
library(dplyr)
library(reshape2)
library(ggplot2)

load("NOMeData.Hf03.filtered.RData")
load("NOMeData.Hf03.gr.RData")
load("introns.ret.repl.RData")
load("introns.nonret.repl.RData")

retint.5splicepos.gr<-GRangesList()
for (i in 1:length(introns.ret.repl)) {
  retint.5splicepos.gr[[i]]<-makeGRangesFromDataFrame(introns.ret.repl[[i]], seqnames.field = "Chromosome", start.field = "End", end.field = "End", strand.field = "Strand")[introns.ret.repl[[i]]$Strand == "+"]
}
retint.5spliceneg.gr<-GRangesList()
for (i in 1:length(introns.ret.repl)) {
  retint.5spliceneg.gr[[i]]<-makeGRangesFromDataFrame(introns.ret.repl[[i]], seqnames.field = "Chromosome", start.field = "Start", end.field = "Start", strand.field = "Strand")[introns.ret.repl[[i]]$Strand == "-"]
}

retint.5splice.gr<-GRangesList()
for (i in 1:length(introns.ret.repl)) {
  retint.5splice.gr[[i]]<-sort(c(retint.5splicepos.gr[[i]], retint.5spliceneg.gr[[i]]), ignore.strand = TRUE)
}

retint.5ss.gr<-GRangesList()
for (i in 1:length(introns.ret.repl)) {
  retint.5ss.gr[[i]]<-expandRange(retint.5splice.gr[[i]], 100, 99)
}

retint.targetwin.5ss.gr<-GRangesList()
for (i in 1:length(introns.ret.repl)) {
  retint.targetwin.5ss.gr[[i]]<-makeWindows(retint.5ss.gr[[i]], k = 200)
  retint.targetwin.5ss.gr[[i]]<-resize(retint.targetwin.5ss.gr[[i]], width = 10, fix = "end", use.names = TRUE)
}

retint.5ss.olaps<-list()
for (i in 1:length(NOMe.Hf03.gr)) {
  retint.5ss.olaps[[i]]<-findOverlaps(NOMe.Hf03.gr[[i]], retint.targetwin.5ss.gr[[i]], ignore.strand = T) #finding overlaping sites between methylation data and windows of interest, specify required Meth.gr object#
}

retint.5ss.occup<-list()
for (i in 1:length(NOMe.Hf03.gr)) {
  retint.5ss.occup[[i]]<-NOMeData.Hf03.filtered[[i]][queryHits(retint.5ss.olaps[[i]]),] #filtering methylation dataset with sites that overlap with area of interest#
  retint.5ss.occup[[i]]$Group<-subjectHits(retint.5ss.olaps[[i]]) ##adding new column identifying the matches in windows##
}

retint.5ss.cov<-list()
retint.5ss.GCHmeth<-list()
retint.5ss.GCHoccup<-list()
retint.5ss.Ts<-list()
for (i in 1:length(NOMe.Hf03.gr)) {
  retint.5ss.cov[[i]]<-aggregate(Coverage ~ Group, data = retint.5ss.occup[[i]], FUN = mean)
  retint.5ss.GCHmeth[[i]]<-aggregate(GCH_Meth ~ Group, data = retint.5ss.occup[[i]], FUN = mean)
  retint.5ss.GCHoccup[[i]]<-aggregate(Occupancy ~ Group, data = retint.5ss.occup[[i]], FUN = mean)
  retint.5ss.Ts[[i]]<-aggregate(Unmeth ~ Group, data = retint.5ss.occup[[i]], FUN = mean)
}

background = NA
minusStrand<-list()
for (i in 1:length(retint.5ss.gr)) {
  minusStrand[[i]] <- which(as.character(strand(retint.5ss.gr[[i]])) == '-')
}

retint.5ssvec<-list() #creating vector for coverage
for (i in 1:3) {
  retint.5ssvec[[i]] = rep(background, length(retint.targetwin.5ss.gr[[i]])) ##creating a vector with default value NA##
  retint.5ssvec[[i]][ as.numeric(retint.5ss.cov[[i]]$Group) ] = retint.5ss.cov[[i]]$Coverage ##change a default value to the the sum of methylated cytosines value##
}

retint.5ss.mat.cov<-list() #creating matrix for coverage
for (i in 1:3) {
  for (j in c(3,6,9)){
  retint.5ss.mat.cov[[i]]<-matrix(retint.5ssvec[[i]], nrow = length(retint.5ss.gr[[i]]), ncol = 200, byrow = T)
  retint.5ss.mat.cov[[i]][minusStrand[[i]],] <- retint.5ss.mat.cov[[i]][minusStrand[[i]], 200:1]
  colnames(retint.5ss.mat.cov[[i]])<-c(1:200)
  rownames(retint.5ss.mat.cov[[i]])<-c(introns.ret.repl[[j]]$IntronID)
}}


Hf03.BlTN.retint.500bp.5ss.occup.mat<-retint.5ss.mat.GCHoccup[[1]]
Hf03.BlCM.retint.500bp.5ss.occup.mat<-retint.5ss.mat.GCHoccup[[2]]
Hf03.BlEM.retint.500bp.5ss.occup.mat<-retint.5ss.mat.GCHoccup[[3]]
save(Hf03.BlTN.retint.500bp.5ss.occup.mat, file = "Hf03.BlTN.retint.500bp.5ss.occup.mat.RData") ###matrix is for +/- 100 bp around a splice site
save(Hf03.BlCM.retint.500bp.5ss.occup.mat, file = "Hf03.BlCM.retint.500bp.5ss.occup.mat.RData")
save(Hf03.BlEM.retint.500bp.5ss.occup.mat, file = "Hf03.BlEM.retint.500bp.5ss.occup.mat.RData")

Hf03.BlTN.retint.500bpsmooth.5ss.NOMeGCH.mat<-retint.5ss.mat.NOMeGCH[[1]]
Hf03.BlCM.retint.500bpsmooth.5ss.NOMeGCH.mat<-retint.5ss.mat.NOMeGCH[[2]]
Hf03.BlEM.retint.500bpsmooth.5ss.NOMeGCH.mat<-retint.5ss.mat.NOMeGCH[[3]]
save(Hf03.BlTN.retint.500bpsmooth.5ss.NOMeGCH.mat, file = "Hf03.BlTN.retint.500bpsmooth.5ss.NOMeGCH.mat.RData") ###matrix is for +/- 500 bp surround a splice
save(Hf03.BlCM.retint.500bpsmooth.5ss.NOMeGCH.mat, file = "Hf03.BlCM.retint.500bpsmooth.5ss.NOMeGCH.mat.RData")
save(Hf03.BlEM.retint.500bpsmooth.5ss.NOMeGCH.mat, file = "Hf03.BlEM.retint.500bpsmooth.5ss.NOMeGCH.mat.RData")

Hf03.BlTN.retint.500bpsmooth.5ss.NOMeTs.mat<-retint.5ss.mat.Ts[[1]]
Hf03.BlCM.retint.500bpsmooth.5ss.NOMeTs.mat<-retint.5ss.mat.Ts[[2]]
Hf03.BlEM.retint.500bpsmooth.5ss.NOMeTs.mat<-retint.5ss.mat.Ts[[3]]
save(Hf03.BlTN.retint.500bpsmooth.5ss.NOMeTs.mat, file = "Hf03.BlTN.retint.500bpsmooth.5ss.NOMeTs.mat.RData") ###matrix is for +/- 500 bp surround a splice
save(Hf03.BlCM.retint.500bpsmooth.5ss.NOMeTs.mat, file = "Hf03.BlCM.retint.500bpsmooth.5ss.NOMeTs.mat.RData")
save(Hf03.BlEM.retint.500bpsmooth.5ss.NOMeTs.mat, file = "Hf03.BlEM.retint.500bpsmooth.5ss.NOMeTs.mat.RData")


Hf03.BlTN.retint.500bpsmooth.5ss.NOMecov.mat<-retint.5ss.mat.cov[[1]]
Hf03.BlCM.retint.500bpsmooth.5ss.NOMecov.mat<-retint.5ss.mat.cov[[2]]
Hf03.BlEM.retint.500bpsmooth.5ss.NOMecov.mat<-retint.5ss.mat.cov[[3]]
save(Hf03.BlTN.retint.500bpsmooth.5ss.NOMecov.mat, file = "Hf03.BlTN.retint.500bpsmooth.5ss.NOMecov.mat.RData") ###matrix is for +/- 500 bp surround a splice
save(Hf03.BlCM.retint.500bpsmooth.5ss.NOMecov.mat, file = "Hf03.BlCM.retint.500bpsmooth.5ss.NOMecov.mat.RData")
save(Hf03.BlEM.retint.500bpsmooth.5ss.NOMecov.mat, file = "Hf03.BlEM.retint.500bpsmooth.5ss.NOMecov.mat.RData")

##same is repeated for non-retained introns
## for 3' ss and the middle of introns

##plot average GCH methylation profiles
mean1<-colMeans(retint.5ss.mat.NOMeGCH[[1]], na.rm = T)
mean2<-colMeans(nonretint.5ss.mat.NOMeGCH[[1]], na.rm = T)
averages<-rbind(mean1, mean2)
data<-melt(averages)
ggplot(data, aes(Var2, value, color = Var1, linetype = NFR)) +
  theme_bw(base_size = 20) +
  geom_line(size = 2) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  labs(y = "GCH methylation") +
  ylim(0,50)+
  theme(axis.title.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 15), legend.background = element_blank()) +
  scale_x_continuous(breaks = c(31, 51, 70), labels = c("-200 bp", "5' ss", "+200 bp")) +
  scale_colour_manual(values = c("#CB2314", "#273046"))+
  geom_vline(xintercept = 31, color = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = 51, color = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = 70, color = "darkgrey", linetype = "dashed") +
  theme(plot.margin=unit(c(0.5,0.7,0.5,0.5),"cm"))
