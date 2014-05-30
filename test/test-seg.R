library(Rcpp)
library(GenomicRanges)
library(PICS)
sourceCpp("~/workspace/git/CPPICS/R/seg.cpp")
source("R/seg.R")
gr <- readRDS(file="test/gr.rds")
grC <- readRDS(file="test/grC.rds")

seg <- segmentPICS(gr)
res <- SEGCPP(gr)

unlist(res$chr21$yF)
unlist(lapply(seg@List, function(x){x@yF}))

yFs <- function(x){
  if(class(x)  == "segReadsList"){
    return(unlist(lapply(seg@List, function(X){X@yF)))
  } else {
    return(unlist(x[["chr21"]][["yF"]]))
  }
}

plot(yFs(res), yFs(seg))


source("R/AllClasses.R")
sourceCpp("R/segS4.cpp")
load("call_CPP_seg.rda")
segList<-CPP_seg(dList, cList, minReadsInRegion, minLregion, minReads, step, width, maxStep, minDist, verbose)


min       lq   median       uq      max neval
454.9944 459.4684 463.9398 472.3572 556.4697   100

min       lq   median       uq      max neval
437.1006 439.4241 442.1002 449.4767 566.0396   100

min       lq   median       uq      max neval
453.1137 458.4022 464.6844 477.6191 575.0005   100



