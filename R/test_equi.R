library(Rcpp)
library(GenomicRanges)
library(PICS)
load("dat.rda")
source("AllClasses.R")
source("seg.R")
sourceCpp("seg.cpp")
res<-SEGCPP(gr)

seg<-segmentPICS(gr)


#CPPICS
nyF1<-sum(unlist(lapply(res$chr21$yF, length)))
ncF1<-sum(unlist(lapply(res$chr21$cF, length)))
first1<-res[[1]]$yF[[1]]


#PICS
nyF2<-sum(unlist(lapply(seg[[1]], function(x){length(x@yF)})))
ncF2<-sum(unlist(lapply(seg[[1]], function(x){length(x@cF)})))
first2<-seg[[1]][[1]]@yF

#
nyF1==nyF2
ncF1==ncF2

