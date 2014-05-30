void callRegionsL(int *center, int *nProbes, int *width, int *scoreF, int *scoreR, int *scoreRegionF, int *scoreRegionR, int *cutoff, int *StartRegion, int *EndRegion, int *nRegions, int maxStep, int kStep, int minL)
{
        //Rprintf("\n Regions shorter than %i bps, will not be called. \n", minL);
        int p=0, w=0, ww=0, start=0, minScore=0,  min=0,  cutCenter=0, pp;
        //p             :the index of current center (outer loop)
        //ww    :the index of temperory last center in the regions
        //w             :the index of current center (inner loop)
        //start :the index of first center in current segment
        //minScore: save the minimium score of current region
        //min   :save the index of center correponding to the min score, used for cut window in ceter
        //cutCenter: a flag indicating if last region was cut inthe center of window
        //kStep :INTEGER(width)/INTEGER(step), when cut in the center, I want the index of start window is cut point+kStep

        *nRegions=0;

        //Rprintf("nProbes=%i \n",*nProbes);
        while(p < *nProbes)
        {
                //Rprintf("Start loop: p=%i \n",p);
                if(((scoreF[p]>=*cutoff) && (scoreR[p]>=*cutoff))|| (cutCenter==1)) //current window has enough reads to be recorded
                {
                        //maxF=scoreF[p];
                        //maxR=scoreR[p];     
                        (*nRegions)++;

                        /*set start point of the region according to if last region was cut in the center*/
                        if (cutCenter==0)
                        {
                                StartRegion[*nRegions-1]=center[p]-*width/2;
                                //update new minscore
                                start=p;
                                min=start;
                                minScore=imin2(scoreF[start], scoreR[start]);
                                //Rprintf("\n Start a new segment \n");
                        }else
                        {
                                StartRegion[*nRegions-1]=EndRegion[*nRegions-2]+1;
                                //update new mini score
                                start=min+kStep;
                                min=start;
                                minScore=imin2(scoreF[start], scoreR[start]);
                                for (pp=start; pp<=p; pp++) {
                                        if (scoreF[pp]<minScore)
                                        {
                                                minScore= scoreF[pp];
                                                min=pp;
                                        }
                                        if (scoreR[pp]<minScore)
                                        {
                                                minScore= scoreR[pp];
                                                min=pp;
                                        }
                                        pp++;

                                }
                                //Rprintf("\n Continue from last cut window \n");
                        }

                        //Rprintf("StartRegion[%i]=%i \n", *nRegions-1,StartRegion[*nRegions-1]);
                        //Rprintf("min=%i, \t start=%i, \t p=%i, \t kStep=%i \n", min, start, p, kStep);


                        w=p+1;
                        ww=p;


                        while(((w-start)<=maxStep) && ((center[w]-center[ww])<=*width) && (w < *nProbes))
                        {
                                if((scoreF[w]>=*cutoff) & (scoreR[w]>=*cutoff))
                                {
                                        ww=w;
                                        if (scoreF[w]<minScore)
                                        {
                                                minScore= scoreF[w];
                                                min=w;
                                        }
                                        if (scoreR[w]<minScore)
                                        {
                                                minScore= scoreR[w];
                                                min=w;
                                        }

                                }
                                w++;
                        }

                        if (w == *nProbes)
                        {
                                EndRegion[*nRegions-1]=center[ww]+*width/2;
                        }else if ((w-start)>maxStep)
                        {
                                EndRegion[*nRegions-1]=center[min];
                                cutCenter=1;
                        }else {
                                EndRegion[*nRegions-1]=center[ww]+*width/2;
                                cutCenter=0;
                        }
                        //Rprintf("EndRegion[%i]=%i \n", *nRegions-1,EndRegion[*nRegions-1]);
                        if ((EndRegion[*nRegions-1]-StartRegion[*nRegions-1]) < minL)
                        {
                                //Rprintf("\n %i - %i < %i, go back one step \n", EndRegion[*nRegions-1], StartRegion[*nRegions-1], minL);
                                (*nRegions)--;
                        }


                        p=w;
                } else // do nothing, look at next window
                {
                        p++;
                }/*else if (cutCenter==1) //current window have too few reads, but last regions was cut
                  { 
                  (*nRegions)++;
                  StartRegion[*nRegions-1]=EndRegion[*nRegions-2]+1;
                  EndRegion[*nRegions-1]=center[p-1]+*width/2;
                  cutCenter=0;
                  Rprintf("StartRegion[%i]=%i \n", *nRegions-1,StartRegion[*nRegions-1]);
                  Rprintf("EndRegion[%i]=%i \n", *nRegions-1,EndRegion[*nRegions-1]);
                  p=p-1+kStep;
                  if ((EndRegion[*nRegions-1]-StartRegion[*nRegions-1]) < minL) 
                  {
                  Rprintf("\n %i - %i < %i, go back one step \n", EndRegion[*nRegions-1], StartRegion[*nRegions-1], minL);
                  (*nRegions)--;
                  }
                  Rprintf("\n End last cut window \n");
                  }*/
                //Rprintf("end loop: p=%i \n",p);
        }
}
