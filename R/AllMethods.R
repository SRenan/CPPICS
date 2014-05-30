## show and summary methods
setMethod("show", "segReads",
          function(object)
      {
          cat("Object of class ",as.character(class(object)),"\n")
          cat("This object has the following slots: \n")
          cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
          #cat("yR, yF, cF, cR, map\n")
      })

setMethod("show", "segReadsList",
          function(object)
      {
          cat("Object of class",as.character(class(object))," with the following slots: \n")
          cat(paste(names(getSlots(class(object))),collapse=", "),"\n")
          #cat("List, paraSW, N, Nc\n")
          cat("List is a list of 'segReads' ojects, each of which has the following slots:\n")
          cat("yR, yF, cR, cF, map, chr\n")
      })


setGeneric("map", function(x, ...) standardGeneric("map"))
setMethod("map", "segReads",
          function(x)
          {
            if(is.null(x@map) | (nrow(x@map)==0))
            {
              return(0);
            }
            else
            {
              n<-nrow(x@map)
              m<-min(x@yF[1],x@yR[1],x@map[1,1]);M<-max(tail(x@yF,1),tail(x@yR,1),x@map[n,2]);
              return(sum(diff(t(x@map)))/max(M-m,1));
            }
})

setMethod("map", "segReadsList",
          function(x)
          {
            ans<-.Call("getMap", x@List, PACKAGE="PICS");
            return(ans)
          }
)

setMethod("length", "segReadsList",
          function(x)
          {
            return(length(x@List))
})

setMethod("summary", "segReadsList",
          function(object)
      {
          cat("** Experiment information ** \n")
          cat("Chromosomes interogated: ")
          cat(unique(unlist(lapply(object@List,function(obj){obj@chr}))),"\n")
          cat("Number of reads")
          cat(" in IP: ",object@N," and in control: ",object@Nc,"\n")
          cat("** Segmentation parameters ** \n")
          cat("The following settings were used:\n")          
          cat("  Sliding window half width: ", object@paraSW$width,"\n")
          cat("  Step size: ", object@paraSW$step,"\n")          
          cat("  Minimum number of reads: ", object@paraSW$minReads,"\n")
          cat("** Segmentation summary ** \n")                    
          cat("Number of segmented regions:",length(object@List),"\n")
          cat("Summary on the number of Forward/Reverse reads per region:\n")
          cat("  Forward:\n") 
          cat("  ")
          tempF<-lapply(object@List,function(obj){length(obj@yF)})
          print(summary(as.integer(unlist(tempF))))
          cat("  Reverse:\n") 
          cat("  ")
          tempR<-lapply(object@List,function(obj){length(obj@yR)})
          print(summary(as.integer(unlist(tempR))))
          cat("Summary on the number of control Forward/Reverse reads per region:\n")
          cat("  Forward:\n") 
          cat("  ")
          tempF<-lapply(object@List,function(obj){length(obj@cF)})
          print(summary(as.integer(unlist(tempF))))
          cat("  Reverse:\n") 
          cat("  ")
          tempR<-lapply(object@List,function(obj){length(obj@cR)})
          print(summary(as.integer(unlist(tempR))))                    
          tempMap<-map(object)
          cat("** Mappability summary **\n")
          cat("Non mappable intervals cover an average ", mean(unlist(tempMap)),"% of all regions \n")                  
      })


setMethod("summary", "segReads",
      function(object)
      {
        m<-min(object@yF[1],object@yR[1])
        M<-max(tail(object@yF,1),tail(object@yR,1))
        cat("** Region summary ** \n")
        cat("Summary on Forward reads:\n")
        print(summary(object@yF,digits=100))
        cat("Summary on Reverse reads:\n")
        print(summary(object@yR,digits=100))
        cat("Summary on control Forward reads:\n")
        print(summary(object@cF,digits=100))
        cat("Summary on control Reverse reads:\n")
        print(summary(object@cR,digits=100))
        cat("Non mappable intervals cover ", sum(diff(t(object@map)))/(M-m),"% of the region \n")
      })


setMethod("[","segReadsList",
		function(x,i, j,..., drop=FALSE)
		{
			if(missing(i))
			{
				return(x)
			}
			if(!missing(j))
			{
			  stop("incorrect number of dimensions")
			}
      else
      {
        segReadsList(x@List[i],x@paraSW,x@N,x@Nc)
      }
		})

setMethod("[[","segReadsList",
    function(x, i, j, ..., exact = TRUE)
    {
      if(length(i) != 1)
      {
        stop("subscript out of bounds (index must have length 1)")
      }
      if(missing(i))
      {
        return(x)
      }
      if(!missing(j))
      {
        stop("incorrect number of dimensions")
      }
      x@List[[i]]
})
