SEGCPP<-function(data, dataC=NULL, minReads=2, minReadsInRegion=3, minLregion=100, maxLregion=0, minDist=0, step=20, width=250, verbose=2){
  #TODO:
  #minDist = 25 in PING
  #maxLregion
  #package? Use a boolean instead of a string
  .checkInput(data, dataC)
  
  if(maxLregion > 0){
    maxStep <- (maxLregion-2*width)/step
  } else{
    maxStep <- 0
  }
  dList<-.format4cpp(data)
  if(is.null(dataC)){
    cList<-vector('list', length(seqlevels(data)))
    names(cList)<-seqlevels(data)
    cList<-lapply(cList, function(x){x<-list("F"=numeric(), "R"=numeric())})
  } else{
    cList<-.format4cpp(dataC)
  }
  #segList list of chr
    #print(dList)
    #print(cList)
  segList<-CPP_seg(dList, cList, minReadsInRegion, minLregion, minReads, step, width, maxStep, minDist, verbose)
  #data, dataC, start, end, jitter, paraSW, maxStep, minLregion, pPackage, PACKAGE
  #             ----map----         >step-width-minReads
  if(is.null(unlist(segList))){
    stop("No candidate region found, you should decrease 'minReads'")
  }
  paraSW<-list(step=step, width=width, minReads=minReads)
  #segList<-cpp_cons(segList, paraSW, length(data), length(dataC))
  return(segList)
}

#Prepare a GRanges object to be passed to CPP_seg
.format4cpp<-function(GRobject){
  chrs<-levels(seqnames(GRobject))
  ret<-vector('list', length(chrs))
  names(ret) <- chrs
  for(chr in chrs){
    GRchr<-GRobject[seqnames(GRobject)==chr]
    ret[[chr]][["F"]]<-start(GRchr[strand(GRchr)=="+"])
    ret[[chr]][["R"]]<-end(GRchr[strand(GRchr)=="-"])
  }
  return(ret)
}

.checkInput<-function(data, dataC){
  if(!is(data, "GRanges")){
    stop("The input data should be an object of class 'GRanges'. You provided: ", class(data))
  }
  if(!is.null(dataC)){
    if(!is(dataC, "GRanges") ){
      stop("The input dataC should be an object of class 'GRanges' or NULL. You provided: ", class(dataC))
    }
    if(length(levels(seqnames(dataC))) != length(levels(seqnames(data)))){
      stop("The data and control do not have the same number of chromosomes.\n",
        "data: ", length(levels(seqnames(data))), "  dataC: ", length(levels(seqnames(dataC))))
    }
  }
}

#Switch the level of the list to have chr - regIdx - yF/cF/cR
.switchList<-function(obj){
  total_reg_cnt<-sum(unlist(lapply(obj, function(x){ length(x$yF) })))
  ret <- vector('list', total_reg_cnt)
  retIdx <- 1
  for(chrIdx in 1:length(obj)){
    reg <- obj[[chrIdx]]
    nReg <- length(reg$yF)
    for(i in 1:nReg){
      #ret[[retIdx]]<-list(yF=reg$yF[[i]], yR=reg$yR[[i]], cF=reg$cF[[i]], cR=reg$cR[[i]], map=reg$map[[i]], chr=names(obj)[[chrIdx]])
      ret[[retIdx]]<-list(yF=reg$yF[[i]], yR=reg$yR[[i]], cF=reg$cF[[i]], cR=reg$cR[[i]], chr=names(obj)[[chrIdx]])
      retIdx <- retIdx+1
    }
  }
  return(ret)
}
