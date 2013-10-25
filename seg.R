RCPP<-function(data, dataC=NULL, minReads=2, minDist=0, step=20, minLength=100, width=250, verbose=2){
  #minDist = 25 in PING
  .checkInput(data, dataC)
  
  maxStep=0;


  #if(is.null(dataC)){
    #cF
  #} else{
    #cF<-start(data[strand(dataC)=="+"])
    #cR<-end(data[strand(dataC)=="-"])
  #}
  #cList<-list(list("cF"=cF, "cR"=cR))
  dList<-.format4cpp(data)
  cList<-.format4cpp(dataC)
  resList<-CPP_seg(dList, cList, minReads, step, width, maxStep, minLength, minDist, verbose)
  #resList<-CPP_seg(yF, yR, minReads, step, width, maxStep, minLength, minDist)
  return(resList)
}

#Prepare a GRanges object to be passed to CPP_seg
.format4cpp<-function(GRobject){
  chrs<-levels(seqnames(GRobject))
  ret<-vector('list', length(chrs))
  names(ret) <- chrs
  for(chr in chrs){
    GRchr<-GRobject[seqnames(GRobject)==chr]
    ret[[chr]]<-vector('list',2)
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
