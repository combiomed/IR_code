## ChIP-seq files processed in MACS2 are used in this step
## GRanges for 5' splice site are generated as described in 03_Features_from_WGBS

load("introns.ret.samp.RData")

Hm01.BlMo.ChIPseq.files<-c('ChIP/Hm01_BlMo_H3K27ac_peaks.xls',
                           'ChIP/Hm01_BlMo_H3K27me3_peaks.xls',
                           'ChIP/Hm01_BlMo_H3K36me3_peaks.xls',
                           'ChIP/Hm01_BlMo_H3K4me1_peaks.xls',
                           'ChIP/Hm01_BlMo_H3K4me3_peaks.xls',
                           'ChIP/Hm01_BlMo_H3K9me3_peaks.xls')

Hm01.BlMo.ChIP.narrowpeaks<-list()
for (i in 1:length(Hm01.BlMo.ChIPseq.files)) {
  Hm01.BlMo.ChIP.narrowpeaks[[i]]<-read.table(Hm01.BlMo.ChIPseq.files[i], header = T)
  Hm01.BlMo.ChIP.narrowpeaks[[i]]$chr<-sub("^", "chr", Hm01.BlMo.ChIP.narrowpeaks[[i]]$chr)
  Hm01.BlMo.ChIP.narrowpeaks[[i]]$strand<-c("*")
  Hm01.BlMo.ChIP.narrowpeaks[[i]]$chr<-factor(Hm01.BlMo.ChIP.narrowpeaks[[i]]$chr, levels=chrOrder)
  Hm01.BlMo.ChIP.narrowpeaks[[i]]<-Hm01.BlMo.ChIP.narrowpeaks[[i]][order(Hm01.BlMo.ChIP.narrowpeaks[[i]]$chr),]
}

Hm01.BlMo.pileup.cutoff<-list()
for (i in 1:length(Hm01.BlMo.ChIPseq.files)) {
  Hm01.BlMo.pileup.cutoff[i]<-mean(Hm01.BlMo.ChIP.narrowpeaks[[i]]$pileup)+sd(Hm01.BlMo.ChIP.narrowpeaks[[i]]$pileup)
}

Hm01.BlMo.ChIP.narrowpeaks.gr<-GRangesList()
Hm01.BlMo.ChIP.summitpeaks.gr<-GRangesList()
for (i in 1:length(Hm01.BlMo.ChIPseq.files)) {
  Hm01.BlMo.ChIP.narrowpeaks.gr[[i]]<-makeGRangesFromDataFrame(Hm01.BlMo.ChIP.narrowpeaks[[i]], seqnames.field = "chr", start.field = "start", end.field = "end", strand.field = "strand")
  Hm01.BlMo.ChIP.summitpeaks.gr[[i]]<-makeGRangesFromDataFrame(Hm01.BlMo.ChIP.narrowpeaks[[i]], seqnames.field = "chr", start.field = "abs_summit", end.field = "abs_summit", strand.field = "strand")
  Hm01.BlMo.ChIP.summitpeaks.gr[[i]]<-expandRange(Hm01.BlMo.ChIP.summitpeaks.gr[[i]], 50, 50)
}

Hm01.BlMo.ChIPPeaks.olaps<-list()
Hm01.BlMo.ChIPPeaks.summit.olaps<-list()
for (i in 1:length(Hm01.BlMo.ChIPseq.files)) {
  Hm01.BlMo.ChIPPeaks.olaps[[i]]<-findOverlaps(Hm01.BlMo.ChIP.narrowpeaks.gr[[i]], retint.5ss.gr[[1]], ignore.strand = T)
  Hm01.BlMo.ChIPPeaks.summit.olaps[[i]]<-findOverlaps(Hm01.BlMo.ChIP.summitpeaks.gr[[i]], retint.5ss.gr[[1]], ignore.strand = T)
}

Hm01.BlMo.ChIPPeaks.retint<-list()
Hm01.BlMo.ChIPPeaks.summit.retint<-list()
for (i in 1:length(Hm01.BlMo.ChIPseq.files)) {
  Hm01.BlMo.ChIPPeaks.retint[[i]]<-Hm01.BlMo.ChIP.narrowpeaks[[i]][queryHits(Hm01.BlMo.ChIPPeaks.olaps[[i]]),] #filtering methylation dataset with sites that overlap with area of interest#
  Hm01.BlMo.ChIPPeaks.retint[[i]]$Group<-subjectHits(Hm01.BlMo.ChIPPeaks.olaps[[i]]) ##adding new column identifying the matches in windows##
  Hm01.BlMo.ChIPPeaks.summit.retint[[i]]<-Hm01.BlMo.ChIP.narrowpeaks[[i]][queryHits(Hm01.BlMo.ChIPPeaks.summit.olaps[[i]]),] #filtering methylation dataset with sites that overlap with area of interest#
  Hm01.BlMo.ChIPPeaks.summit.retint[[i]]$Group<-subjectHits(Hm01.BlMo.ChIPPeaks.summit.olaps[[i]]) ##adding new column identifying the matches in windows##
}

Hm01.BlMo.ChIPPeaks.score<-list()
for (i in 1:length(Hm01.BlMo.ChIPseq.files)) {
  Hm01.BlMo.ChIPPeaks.score[[i]]<-aggregate(pileup ~ Group, data = Hm01.BlMo.ChIPPeaks.summit.retint[[i]], FUN = mean)
}

introns.ret.samp[[1]]$H3K27ac_olap_100bp5ss<-ifelse(match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.retint[[1]]$Group), 1, 0)
introns.ret.samp[[1]]$H3K27ac_pileup_100bp5ss<-Hm01.BlMo.ChIPPeaks.score[[1]]$pileup[match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.score[[1]]$Group)]
introns.ret.samp[[1]]$H3K27me3_olap_100bp5ss<-ifelse(match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.retint[[2]]$Group), 1, 0)
introns.ret.samp[[1]]$H3K27me3_pileup_100bp5ss<-Hm01.BlMo.ChIPPeaks.score[[2]]$pileup[match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.score[[2]]$Group)]
introns.ret.samp[[1]]$H3K36me3_olap_100bp5ss<-ifelse(match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.retint[[3]]$Group), 1, 0)
introns.ret.samp[[1]]$H3K36me3_pileup_100bp5ss<-Hm01.BlMo.ChIPPeaks.score[[3]]$pileup[match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.score[[3]]$Group)]
introns.ret.samp[[1]]$H3K4me1_olap_100bp5ss<-ifelse(match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.retint[[4]]$Group), 1, 0)
introns.ret.samp[[1]]$H3K4me1_pileup_100bp5ss<-Hm01.BlMo.ChIPPeaks.score[[4]]$pileup[match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.score[[4]]$Group)]
introns.ret.samp[[1]]$H3K4me3_olap_100bp5ss<-ifelse(match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.retint[[5]]$Group), 1, 0)
introns.ret.samp[[1]]$H3K4me3_pileup_100bp5ss<-Hm01.BlMo.ChIPPeaks.score[[5]]$pileup[match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.score[[5]]$Group)]
introns.ret.samp[[1]]$H3K9me3_olap_100bp5ss<-ifelse(match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.retint[[6]]$Group), 1, 0)
introns.ret.samp[[1]]$H3K9me3_pileup_100bp5ss<-Hm01.BlMo.ChIPPeaks.score[[6]]$pileup[match(introns.ret.samp[[1]]$ID, Hm01.BlMo.ChIPPeaks.score[[6]]$Group)]

introns.ret.samp[[1]]$H3K27ac_olap_100bp5ss<-ifelse(introns.ret.samp[[1]]$H3K27ac_pileup_100bp5ss >= Hm01.BlMo.pileup.cutoff[[1]], 2, introns.ret.samp[[1]]$H3K27ac_olap_100bp5ss)
introns.ret.samp[[1]]$H3K27me3_olap_100bp5ss<-ifelse(introns.ret.samp[[1]]$H3K27me3_pileup_100bp5ss >= Hm01.BlMo.pileup.cutoff[[2]], 2, introns.ret.samp[[1]]$H3K27me3_olap_100bp5ss)
introns.ret.samp[[1]]$H3K36me3_olap_100bp5ss<-ifelse(introns.ret.samp[[1]]$H3K36me3_pileup_100bp5ss >= Hm01.BlMo.pileup.cutoff[[3]], 2, introns.ret.samp[[1]]$H3K36me3_olap_100bp5ss)
introns.ret.samp[[1]]$H3K4me1_olap_100bp5ss<-ifelse(introns.ret.samp[[1]]$H3K4me1_pileup_100bp5ss >= Hm01.BlMo.pileup.cutoff[[4]], 2, introns.ret.samp[[1]]$H3K4me1_olap_100bp5ss)
introns.ret.samp[[1]]$H3K4me3_olap_100bp5ss<-ifelse(introns.ret.samp[[1]]$H3K4me3_pileup_100bp5ss >= Hm01.BlMo.pileup.cutoff[[5]], 2, introns.ret.samp[[1]]$H3K4me3_olap_100bp5ss)
introns.ret.samp[[1]]$H3K9me3_olap_100bp5ss<-ifelse(introns.ret.samp[[1]]$H3K9me3_pileup_100bp5ss >= Hm01.BlMo.pileup.cutoff[[6]], 2, introns.ret.samp[[1]]$H3K9me3_olap_100bp5ss)

introns.ret.samp[[1]]$H3K27ac_olap_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K27ac_olap_100bp5ss), 0, introns.ret.samp[[1]]$H3K27ac_olap_100bp5ss)
introns.ret.samp[[1]]$H3K27ac_pileup_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K27ac_pileup_100bp5ss), 0, introns.ret.samp[[1]]$H3K27ac_pileup_100bp5ss)
introns.ret.samp[[1]]$H3K27me3_olap_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K27me3_olap_100bp5ss), 0, introns.ret.samp[[1]]$H3K27me3_olap_100bp5ss)
introns.ret.samp[[1]]$H3K27me3_pileup_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K27me3_pileup_100bp5ss), 0, introns.ret.samp[[1]]$H3K27me3_pileup_100bp5ss)
introns.ret.samp[[1]]$H3K36me3_olap_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K36me3_olap_100bp5ss), 0, introns.ret.samp[[1]]$H3K36me3_olap_100bp5ss)
introns.ret.samp[[1]]$H3K36me3_pileup_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K36me3_pileup_100bp5ss), 0, introns.ret.samp[[1]]$H3K36me3_pileup_100bp5ss)
introns.ret.samp[[1]]$H3K4me1_olap_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K4me1_olap_100bp5ss), 0, introns.ret.samp[[1]]$H3K4me1_olap_100bp5ss)
introns.ret.samp[[1]]$H3K4me1_pileup_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K4me1_pileup_100bp5ss), 0, introns.ret.samp[[1]]$H3K4me1_pileup_100bp5ss)
introns.ret.samp[[1]]$H3K4me3_olap_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K4me3_olap_100bp5ss), 0, introns.ret.samp[[1]]$H3K4me3_olap_100bp5ss)
introns.ret.samp[[1]]$H3K4me3_pileup_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K4me3_pileup_100bp5ss), 0, introns.ret.samp[[1]]$H3K4me3_pileup_100bp5ss)
introns.ret.samp[[1]]$H3K9me3_olap_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K9me3_olap_100bp5ss), 0, introns.ret.samp[[1]]$H3K9me3_olap_100bp5ss)
introns.ret.samp[[1]]$H3K9me3_pileup_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$H3K9me3_pileup_100bp5ss), 0, introns.ret.samp[[1]]$H3K9me3_pileup_100bp5ss)

save(introns.ret.samp, file = "introns.ret.samp.RData")

##same is repeated for 3' splice site and the middle of intron
##and non-retained introns
## and the remaining biological replicates
