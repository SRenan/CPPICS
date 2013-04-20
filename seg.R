RCPP<-function(data, minReads=2, minDist=0, step=20, minLength=100, width=250){
  #minDist = 25 in PING
  #data is GRanges with 1chr
  
  maxStep=0;

  ##TODO: Later, make a list of yF and a list of yR. For chrs. (add a char vector for chrnames)
  yF<-start(data[strand(data)=="+"])
  yR<-end(data[strand(data)=="-"])
  #resList<-.Call("CPP_seg", yF, yR)
  resList<-CPP_seg(yF, yR, minReads, step, width, maxStep, minLength, minDist)
  return(resList)
}
