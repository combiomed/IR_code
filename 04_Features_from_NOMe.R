## adding features containing information from NOMe data
## pre-processed NOMe-seq dtaa is used for this step
## GRanges for the splice sites and the middle of introns are generated in the same way as in 03_Features_from_WGBS

load("/spectrum/GSCT/veronikap/DEEP/NOMe/monocytes/NOMeData.BlMo.Hm01.filtered.RData") #due to the file size each biological replicate was processed separately
load("/spectrum/GSCT/veronikap/DEEP/NOMe/monocytes/NOMeData.BlMo.Hm01.gr.RData")
load("/home/veronikap/MONvsMAC/introns.ret.samp.RData")

NOMeData.BlMo.Hm01.filtered$Diff<-(lead(NOMeData.BlMo.Hm01.filtered$Start, 1) - NOMeData.BlMo.Hm01.filtered$Start)
NOMeData.BlMo.Hm01.subset<-subset(NOMeData.BlMo.Hm01.filtered, Diff>1)
NOMeData.BlMo.Hm01.subset.gr<-makeGRangesFromDataFrame(NOMeData.BlMo.Hm01.subset, seqnames.field = "Chromosome", start.field = "Start", end.field = "End", strand.field = "Strand")

retint.5ss.gr<-GRangesList()
for (i in 1:5) {
  retint.5ss.gr[[i]]<-expandRange(retint.5splice.gr[[i]], 100, 99)
}

retint.5ss.olaps<-findOverlaps(NOMeData.BlMo.Hm01.subset.gr, retint.5ss.gr[[1]], ignore.strand = T) #finding overlaping sites between methylation data and windows of interest, specify required Meth.gr object#
retint.5ss.occup<-NOMeData.BlMo.Hm01.subset[queryHits(retint.5ss.olaps),] #filtering methylation dataset with sites that overlap with area of interest#
retint.5ss.occup$Group<-subjectHits(retint.5ss.olaps) ##adding new column identifying the matches in windows##
retint.5ss.GCHs<-aggregate(Coverage ~ Group, data = retint.5ss.occup, FUN = length) # counting the number of GCH sites
retint.5ss.GCHmeth<-aggregate(GCH_Meth ~ Group, data = retint.5ss.occup, FUN = mean) #calculating average GCH methylation
retint.5ss.GCHoccup<-aggregate(Occupancy ~ Group, data = retint.5ss.occup, FUN = mean) $calculating average nucleosome occupancy

introns.ret.samp[[1]]$GCHs_100bp5ss<-retint.5ss.GCHs$Coverage[match(introns.ret.samp[[1]]$ID, retint.5ss.GCHs$Group)]
introns.ret.samp[[1]]$GCHmeth_100bp5ss<-retint.5ss.GCHmeth$GCH_Meth[match(introns.ret.samp[[1]]$ID, retint.5ss.GCHmeth$Group)]
introns.ret.samp[[1]]$GCHoccup_100bp5ss<-retint.5ss.GCHoccup$Occupancy[match(introns.ret.samp[[1]]$ID, retint.5ss.GCHoccup$Group)]
introns.ret.samp[[1]]$GCHs_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$GCHs_100bp5ss), 0, introns.ret.samp[[1]]$GCHs_100bp5ss) # NA values (absence of the signal) are replaced with 0
introns.ret.samp[[1]]$GCHmeth_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$GCHmeth_100bp5ss), 0, introns.ret.samp[[1]]$GCHmeth_100bp5ss)
introns.ret.samp[[1]]$GCHoccup_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$GCHoccup_100bp5ss), 0, introns.ret.samp[[1]]$GCHoccup_100bp5ss)

####Repeating for BlMo Hm05
load("/spectrum/GSCT/veronikap/DEEP/NOMe/monocytes/NOMeData.BlMo.Hm05.filtered.RData")
NOMeData.BlMo.Hm05.filtered$Diff<-(lead(NOMeData.BlMo.Hm05.filtered$Start, 1) - NOMeData.BlMo.Hm05.filtered$Start)
NOMeData.BlMo.Hm05.subset<-subset(NOMeData.BlMo.Hm05.filtered, Diff>1)
NOMeData.BlMo.Hm05.subset.gr<-makeGRangesFromDataFrame(NOMeData.BlMo.Hm05.subset, seqnames.field = "Chromosome", start.field = "Start", end.field = "End", strand.field = "Strand")


####Repeating for BlMa Hm05
load("/spectrum/GSCT/veronikap/DEEP/NOMe/macrophages/NOMeData.BlMa.Hm05.filtered.RData")
NOMeData.BlMa.Hm05.filtered$Diff<-(lead(NOMeData.BlMa.Hm05.filtered$Start, 1) - NOMeData.BlMa.Hm05.filtered$Start)
NOMeData.BlMa.Hm05.subset<-subset(NOMeData.BlMa.Hm05.filtered, Diff>1)
NOMeData.BlMa.Hm05.subset.gr<-makeGRangesFromDataFrame(NOMeData.BlMa.Hm05.subset, seqnames.field = "Chromosome", start.field = "Start", end.field = "End", strand.field = "Strand")


##same is repeated for 3' splice site and the middle of intron
##and non-retained introns


## identifying overlaps between the regions of interest and nucleosome localisation, as predicted using gNOMePeaks
retint.samp.gr<-GRangesList() ####creating genomic coordinates for introns
for (i in 1:length(introns.ret.samp)) {
  retint.samp.gr[[i]]<-makeGRangesFromDataFrame(introns.ret.samp[[i]], seqnames.field = "Chromosome", start.field = "Start", end.field = "End", strand.field = "Strand")
}

NOMe.peaks.files<-c("/spectrum/GSCT/veronikap/DEEP/NOMe/monocytes/gnomePeaks/Hm01.BlMo.NOMe.peaks",
                    "/spectrum/GSCT/veronikap/DEEP/NOMe/monocytes/gnomePeaks/Hm05.BlMo.NOMe.peaks",
                    "/spectrum/GSCT/veronikap/DEEP/NOMe/macrophages/gnomePeaks/Hm05.BlMa.NOMe.peaks") ###loading NOMe peaks data

NOMe.peaks<-list()
for (i in 1:length(NOMe.peaks.files)) {
  NOMe.peaks[[i]]<-read.table(NOMe.peaks.files[i], sep = '\t', header = TRUE)
  NOMe.peaks[[i]]$chr<-sub("^", "chr", NOMe.peaks[[i]]$chr)
  NOMe.peaks[[i]]$width<-NOMe.peaks[[i]]$stop-NOMe.peaks[[i]]$start+1
}

NOMe.peaks.subset<-list()
for (i in 1:length(NOMe.peaks.files)) {
  NOMe.peaks.subset[[i]]<-subset(NOMe.peaks[[i]], p.value<0.05 & width>140) ## nucleosome region definition
}

NOMe.peaks.subset.gr<-GRangesList()
for (i in 1:length(NOMe.peaks.files)) {
  NOMe.peaks.subset.gr[[i]]<-makeGRangesFromDataFrame(NOMe.peaks.subset[[i]], seqnames.field = "chr", start.field = "start", end.field = "stop", strand.field = "*")
  NOMe.peaks.subset.gr[[i]]<-unique(NOMe.peaks.subset.gr[[i]])
}

Hm01.BlMo.gnomePeaks.olaps<-findOverlaps(NOMe.peaks.subset.gr[[1]], retint.5ss.gr[[1]], ignore.strand = T)
Hm01.BlMo.gnomePeaks.hits<- pintersect(NOMe.peaks.subset.gr[[1]][queryHits(Hm01.BlMo.gnomePeaks.olaps)], retint.5ss.gr[[1]][subjectHits(Hm01.BlMo.gnomePeaks.olaps)])
Hm01.BlMo.gnomePeaks.percentOverlap <- width(Hm01.BlMo.gnomePeaks.hits) / width(retint.5ss.gr[[1]][subjectHits(Hm01.BlMo.gnomePeaks.olaps)])

Hm01.BlMo.gnomePeaks.tab<-cbind(as.data.frame(Hm01.BlMo.gnomePeaks.olaps), Hm01.BlMo.gnomePeaks.percentOverlap)
introns.ret.samp[[1]]$gnomePeak_pcOlap_100bp5ss<-round(Hm01.BlMo.gnomePeaks.tab$Hm01.BlMo.gnomePeaks.percentOverlap, 4)[match(introns.ret.samp[[1]]$ID, Hm01.BlMo.gnomePeaks.tab$subjectHits)]
introns.ret.samp[[1]]$gnomePeak_olap_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$gnomePeak_pcOlap_100bp5ss), 0, 1)


## identifying overlaps between the regions of interest and nucleosome free regions, as predicted using gNOMePeaks
NOMe.peaks.NFR.files<-c("/spectrum/GSCT/veronikap/DEEP/NOMe/monocytes/gnomePeaks/Hm01.BlMo.NOMe.NFRpeaks",
                        "/spectrum/GSCT/veronikap/DEEP/NOMe/monocytes/gnomePeaks/Hm05.BlMo.NOMe.CpGmethpeaks",
                        "/spectrum/GSCT/veronikap/DEEP/NOMe/macrophages/gnomePeaks/Hm05.BlMa.NOMe.NFRpeaks") ###loading NOMe peaks data

NOMe.peaks.NFR<-list()
for (i in 1:length(NOMe.peaks.NFR.files)) {
  NOMe.peaks.NFR[[i]]<-read.table(NOMe.peaks.NFR.files[i], sep = '\t', header = TRUE)
  NOMe.peaks.NFR[[i]]$chr<-sub("^", "chr", NOMe.peaks.NFR[[i]]$chr)
  NOMe.peaks.NFR[[i]]$width<-NOMe.peaks.NFR[[i]]$stop-NOMe.peaks.NFR[[i]]$start+1
}

NOMe.peaks.NFR.subset<-list()
for (i in 1:length(NOMe.peaks.NFR)) {
  NOMe.peaks.NFR.subset[[i]]<-subset(NOMe.peaks.NFR[[i]], p.value<0.05 & width>100) ## nucleosome free region definition
}

NOMe.peaks.NFR.gr<-GRangesList()
for (i in 1:length(NOMe.peaks.NFR)) {
  NOMe.peaks.NFR.gr[[i]]<-makeGRangesFromDataFrame(NOMe.peaks.NFR.subset[[i]], seqnames.field = "chr", start.field = "start", end.field = "stop", strand.field = "*")
  NOMe.peaks.NFR.gr[[i]]<-unique(NOMe.peaks.NFR.gr[[i]])
}

Hm01.BlMo.gnomePeaks.olaps<-findOverlaps(NOMe.peaks.NFR.gr[[1]], retint.5ss.gr[[1]], ignore.strand = T)
Hm01.BlMo.gnomePeaks.hits<- pintersect(NOMe.peaks.NFR.gr[[1]][queryHits(Hm01.BlMo.gnomePeaks.olaps)], retint.5ss.gr[[1]][subjectHits(Hm01.BlMo.gnomePeaks.olaps)])
Hm01.BlMo.gnomePeaks.percentOverlap <- width(Hm01.BlMo.gnomePeaks.hits) / width(retint.5ss.gr[[1]][subjectHits(Hm01.BlMo.gnomePeaks.olaps)])
Hm01.BlMo.gnomePeaks.tab<-cbind(as.data.frame(Hm01.BlMo.gnomePeaks.olaps), Hm01.BlMo.gnomePeaks.percentOverlap)
introns.ret.samp[[1]]$gnomePeak_pcOlap_NFR_100bp5ss<-round(Hm01.BlMo.gnomePeaks.tab$Hm01.BlMo.gnomePeaks.percentOverlap, 4)[match(introns.ret.samp[[1]]$ID, Hm01.BlMo.gnomePeaks.tab$subjectHits)]
introns.ret.samp[[1]]$gnomePeak_olap_NFR_100bp5ss<-ifelse(is.na(introns.ret.samp[[1]]$gnomePeak_pcOlap_NFR_100bp5ss), 0, 1)

save(introns.ret.samp, file = "/home/veronikap/MONvsMAC/introns.ret.samp.RData")

##same is repeated for 3' splice site and the middle of intron
##and non-retained introns
## and the remaining biological replicates
