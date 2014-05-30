CPPICS
======

PICS using Rcpp

    library(Rcpp)
    load("dat.rda")
    gr
    source("seg.R")
    sourceCpp("seg.cpp")
    res<-RCPP(gr)

