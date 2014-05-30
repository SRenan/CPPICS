SEXP segRPING(SEXP chr, SEXP dataF, SEXP dataR, SEXP contF, SEXP contR, SEXP StartRegion, SEXP EndRegion, SEXP StartMap, SEXP EndMap, SEXP jitter, int nRegions)
{
	int nF = length(dataF), nR = length(dataR), ncF = length(contF), ncR = length(contR), nMap = length(StartMap);
	int i=0, j=0, pF=0, pR=0, pcF=0, pcR=0, pM=0;
	int indStart=0, indEnd=0, indStartM=0, indEndM=0, indStartF=0, indEndF=0,indStartR=0, indEndR=0, *rmap, nProtect;
	int minLocF, maxLocF, minLocR, maxLocR, temp; //temparoryly save boundary each regions
	SEXP ans, map, seg, yF, yR, cF, cR, classDef;
	SEXP name;

	int ext=50; // extend each segments by ext bps on each side, only forward/reverse reads exist in the extended regions on the left/right
	
	PROTECT(name=NEW_CHARACTER(1));
	SET_STRING_ELT(name,0,mkChar(CHAR(chr)));
	GetRNGstate();
	/** Define the list for the output **/
	PROTECT(ans = NEW_LIST(nRegions));
	for(i=0;i<nRegions;i++)
	{
		/** Initialize the protects **/
		nProtect=0;
		
		if (pM>0) pM--; //if PM>0 we want to move the counter one step back since last unmappable regions might overlap with two segments
		/*Process mappability profile*/
		if((pM<nMap) & (nMap>0))
		{
			while((pM<nMap) & (INTEGER(EndMap)[pM]<INTEGER(StartRegion)[i]))
			{
				pM++;
			}
			/** Makes sure the index does not go out of bound */
			indStartM=imin2(pM,nMap);
			minLocF=imax2(INTEGER(StartMap)[indStartM],INTEGER(StartRegion)[i]-ext);
			//minLocR=imax2(INTEGER(StartMap)[indStartM],INTEGER(StartRegion)[i]+ext);
			minLocR=INTEGER(StartRegion)[i]+ext;
			
			
			/** Keep looking **/
			while((pM<nMap) & (INTEGER(StartMap)[pM]<INTEGER(EndRegion)[i]))
			{
				pM++;
			}
			indEndM=imin2(pM,nMap);
			//maxLocF=imin2(INTEGER(EndMap)[indEndM],INTEGER(EndRegion)[i]-ext);
			maxLocR=imin2(INTEGER(EndMap)[indEndM],INTEGER(EndRegion)[i]+ext);
			maxLocF=INTEGER(EndRegion)[i]-ext;
			
			
			PROTECT(map = allocMatrix(INTSXP,indEndM-indStartM,2));
			nProtect++;
			rmap=INTEGER(map);
			
			for(j=indStartM;j<indEndM;j++)
			{
				rmap[j-indStartM+(indEndM-indStartM)*0]=imax2(INTEGER(StartMap)[j],INTEGER(StartRegion)[i]);
				rmap[j-indStartM+(indEndM-indStartM)*1]=imin2(INTEGER(EndMap)[j],INTEGER(EndRegion)[i]);
			}
			
			
		}
		else
		{
			minLocF=INTEGER(StartRegion)[i]-ext;
			minLocR=INTEGER(StartRegion)[i]+ext;
			maxLocF=INTEGER(EndRegion)[i]-ext;
			maxLocR=INTEGER(EndRegion)[i]+ext;
			PROTECT(map = allocMatrix(INTSXP,0,2));
			nProtect++;
		}
		
		//Rprintf("No Truncation: \t\t minLoc[%i]=%i,\t maxLoc[%i]=%i \n",i,INTEGER(StartRegion)[i],i,INTEGER(EndRegion)[i]);
		//Rprintf("Map Truncation: \t minLoc[%i]=%i,\t maxLoc[%i]=%i \n",i,minLoc,i,maxLoc);
		
		/* find the index of 1st yF bounded by regions start,
		 and define "minLoc=min(yF[1], StartMap[1])" */
		while((pF<nF) && (INTEGER(dataF)[pF]<((INTEGER(StartRegion)[i])-ext)))
		{
			pF++;
		}
		indStartF=pF;
		temp=minLocF;
		if (temp>(INTEGER(StartRegion)[i]-ext)) 
		{
			minLocF=imin2(temp, INTEGER(dataF)[indStartF]);
		}else {
			minLocF=imax2(temp, INTEGER(dataF)[indStartF]);
		}

		
	
		/* find the index of 1st yR bounded by minLoc */
		while((pR<nR) && (INTEGER(dataR)[pR]<minLocR))
		{
			pR++;
		}
		indStartR=pR;
		
		
		/* find the index of last yR bounded by regions ends,
		 and define "maxLoc=min(yF[max], EndMap[max])" */
		while((pR<nR) && (INTEGER(dataR)[pR]<=(INTEGER(EndRegion)[i]+ext)))
		{
			pR++;
		}
		indEndR=imin2(pR-1,nR);
		indEndR=imax2(indEndR,indStartR);
		temp=maxLocR;
		if (temp<INTEGER(EndRegion)[i])
		{
			maxLocR=imax2(temp, INTEGER(dataR)[indEndR]);
		}else {
			maxLocR=imin2(temp, INTEGER(dataR)[indEndR]);
		}

		//Rprintf("Reads Truncation: \t minLoc[%i]=%i,\t minLoc[%i]=%i \n",i,minLoc,i,maxLoc);
		
		/* find the index of last yF bounded by maxLoc */
		while((pF<nF) && (INTEGER(dataF)[pF]<=maxLocF))
		{
			pF++;
		}
		indEndF=imin2(pF-1,nF);
		indEndF=imax2(indEndF,indStartF);
		
		
		/*
		 Rprintf("Start: yF[%i]=%i, \n", indStartF,  INTEGER(dataF)[indStartF]);
		 Rprintf("Start: yR[%i]=%i, \n", indStartR,  INTEGER(dataR)[indStartR]);
		 Rprintf("End: yF[%i]=%i,   \n", indEndF,    INTEGER(dataF)[indEndF]);
		 Rprintf("End: yR[%i]=%i,   \n", indEndR,    INTEGER(dataR)[indEndR]);		
		 */
		
		/** Split the data using the start/end index **/
		PROTECT(yF = allocVector(REALSXP,indEndF-indStartF+1));
		nProtect++;
		for(j=indStartF;j<=indEndF;j++)
		{
			REAL(yF)[j-indStartF]=INTEGER(dataF)[j]+rnorm(0,2)*INTEGER(jitter)[0];
		}
		
		PROTECT(yR = allocVector(REALSXP,indEndR-indStartR+1));
		nProtect++;    
		for(j=indStartR;j<=indEndR;j++)
		{
			REAL(yR)[j-indStartR]=INTEGER(dataR)[j]+rnorm(0,2)*INTEGER(jitter)[0];
		}
		
		
		/*Process control data*/
		if((ncF>0) & (ncR>0))
		{
			while((pcF<ncF) & (INTEGER(contF)[pcF]<minLocF))
			{
				pcF++;
			}
			indStart=pcF;
			
			while((pcF<ncF) & (INTEGER(contF)[pcF]<=maxLocF))
			{
				pcF++;
			}
			indEnd=imin2(pcF,ncF);
			
			PROTECT(cF = allocVector(REALSXP,indEnd-indStart));
			nProtect++;      
			for(j=indStart;j<indEnd;j++)
			{
				REAL(cF)[j-indStart]=INTEGER(contF)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
			}
			
			while((pcR<ncR) & (INTEGER(contR)[pcR]<minLocR))
			{
				pcR++;
			}
			indStart=pcR;
			
			while((pcR<ncR) & (INTEGER(contR)[pcR]<=maxLocR))
			{
				pcR++;
			}
			indEnd=imin2(pcR,ncR);
			
			PROTECT(cR = allocVector(REALSXP,indEnd-indStart));
			nProtect++;      
			for(j=indStart;j<indEnd;j++)
			{
				REAL(cR)[j-indStart]=INTEGER(contR)[j]+rnorm(0,.1)*INTEGER(jitter)[0];
			}
		}
		else
		{
			cF = R_NilValue;
			cR = R_NilValue;
		}
		
		classDef=MAKE_CLASS("segReads");
		PROTECT(seg=NEW_OBJECT(classDef));
		nProtect++;    
		SET_SLOT(seg,mkChar("yF"),yF);
		SET_SLOT(seg,mkChar("yR"),yR);
		SET_SLOT(seg,mkChar("cF"),cF);
		SET_SLOT(seg,mkChar("cR"),cR);
		SET_SLOT(seg,mkChar("map"),map);
		SET_SLOT(seg,mkChar("chr"),name);
		SET_VECTOR_ELT(ans,i,seg);
		UNPROTECT(nProtect);
	}
	UNPROTECT(2);
	PutRNGstate();
	return(ans);
}
