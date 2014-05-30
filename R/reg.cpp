void callRegions(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions)
{
  int i=0,p=0;
  int w=0,ww=0,max=0,maxF=0,maxR=0;
  *nRegions=0;
  Rprintf("Width in callRegions: %d\n", *width);

  while(p<*nProbes)
  {
    if((scoreF[p]>=*cutoff) & (scoreR[p]>=*cutoff))
    {
     maxF=scoreF[p];
     maxR=scoreR[p];
      (*nRegions)++;
      StartRegion[*nRegions-1]=center[p]-*width/2;
      w=p+1;
      max=p;
      ww=p;

      while((w<*nProbes) && ((center[w]-center[ww])<=*width))
      {
        if((scoreF[w]>=*cutoff) & (scoreR[w]>=*cutoff))
        {
         if(scoreF[w]>maxF)
         {
           maxF=scoreF[w];
         }

         if(scoreR[w]>maxR)
         {
           maxR=scoreR[w];
         }
          max=w;
          ww=max;
        }
        w++;
      }
     scoreRegionF[*nRegions-1]=maxF;
     scoreRegionR[*nRegions-1]=maxR;
     EndRegion[*nRegions-1]=center[max]+*width/2;
     p=w;
    }
    else
    {
      p++;
    }
  }
}

