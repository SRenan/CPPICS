/**
* TODO:
*   PING > maxLregion=1200 > maxStep>0 > callRegionsL
*   add mappability profile reads. Used only in segR
*   replace minDist by package
*   include PING variations
*   Maybe return a list of chr$regionIdx$yF/yR/cF..
*   Remove the R part in seg_char. Maybe return just the region vector.
*/
#include <Rcpp.h>
#include <iostream>

using namespace std;
struct window{ //Coverage of the sequence by sliding windows
  int center;
  int scoreF;
  int scoreR;
};

struct region{ //The regions
  int start;
  int end;
  int scoreRegionF;
  int scoreRegionR;
  vector<int> regReadsF;
  vector<int> regReadsR;
  vector<int> regContF;
  vector<int> regContR;
};

SEXP seg_chr(vector<int>& yF, vector<int>& yR, vector<int>& contF, vector<int>& contR, int minReadsRegion, int minLregion, int minReads, int step, int width, int maxStep, int minDist, int verbose);
void getWindowsScores(vector<window>& windows, vector<int>& yF, vector<int>& yR, int width, int minDist);
int callRegions(vector<window>& windows, int width, int minReads, vector<region>& regions);
int callLongRegions(vector<window>& windows, int width, int minReads, vector<region>& regions, int maxStep, int kStep, int minLregion);
void segR(vector<int>& yF, vector<int>& yR, vector<int>& contF, vector<int>& contR, vector<window>& windows, vector<region>& regions);


//[[Rcpp::export]]
SEXP CPP_seg(Rcpp::List data, Rcpp::List dataC, int minReadsRegion, int minLregion, int minReads, int step, int width, int maxStep, int minDist, int verbose){
//minReadsRegion : The min number of reads for a region to be kept as candidate
//minLregion     : The min size  of a region
//minReads : The min number of reads for a window to be included in a region
//step     : The step for the sliding window
//width    : Half the length of the sliding window
//maxStep  : Used in callRegionsL
//minDist  : in PING, when reads are too close to the center of the windows, they are not counted in the window score
//verbose  : Print debug values
  if(verbose>1){
    cout<<"step     = "<<step<<endl;
    cout<<"maxStep  = "<<maxStep<<endl;
    cout<<"width    = "<<width<<endl;
    cout<<"minReads = "<<minReads<<endl;
    cout<<"minLregion = "<<minLregion<<endl;
  }
  string chr;
  vector<string> chrs =  data.names();

  Rcpp::List ret(data.size());

  for(int i=0; i < data.size(); i++){
    chr = chrs[i];
    if(verbose>1){
      cout<<"Processing chromosome "<<chr<<endl;
    }
    //access sub list
    Rcpp::List dataList = data[chr];
    Rcpp::List contList = dataC[chr];
    //R to c++ vectors
    std::vector<int> yF = Rcpp::as<std::vector<int> >(dataList["F"]);
    std::vector<int> yR = Rcpp::as<std::vector<int> >(dataList["R"]);
    std::vector<int> contF = Rcpp::as<std::vector<int> >(contList["F"]);
    std::vector<int> contR = Rcpp::as<std::vector<int> >(contList["R"]);
    
    ret[i] = seg_chr(yF, yR, contF, contR, minReadsRegion, minLregion, minReads, step, width, maxStep, minDist, verbose);
  }
  ret.names() = chrs;
  return(Rcpp::wrap(ret));
}

SEXP seg_chr(vector<int>& yF, vector<int>& yR, vector<int>& contF, vector<int>& contR, int minReadsRegion, int minLregion, int minReads, int step, int width, int maxStep, int minDist, int verbose){

  sort(yF.begin(), yF.end());
  sort(yR.begin(), yR.end());
  sort(contF.begin(), contF.end());
  sort(contR.begin(), contR.end());
  
  int m=std::min(yF.at(0), yR.at(0));
  int M=std::max(yF.at(yF.size()-1), yR.at(yR.size()-1));

  int nWindows = (M-m)/step; //Number of windows
  if(verbose>1){
    cout<<"nWindows : "<<nWindows<<endl;
  }

  vector<window> windows(nWindows);
  int i;
  for(i=0; i<nWindows; i++){
    windows.at(i).center=m+i*step;
  }

  //Vector of the maximal possible size for the bigger regions (merged adjacent windows)
//  std::vector<region> regions((M-m)/(2*width));
  std::vector<region> regions(nWindows);

  getWindowsScores(windows, yF, yR, width, minDist); //Sliding window
  int nRegions;
  if(maxStep==0){
    nRegions = callRegions(windows, width, minReads, regions); //Get the actual regions
  } else{
    int kStep = width/step;
    nRegions = callLongRegions(windows, width, minReads, regions, maxStep, kStep, minLregion);
  }
    
  regions.resize(nRegions);

  if(verbose>1){
    cout<<"nRegions : "<<nRegions<<endl;
  }

  if(nRegions>0){
    segR(yF, yR, contF, contR, windows, regions);
  }
  else{
  }

  vector<vector<int> > yFret(nRegions); //vector of vectors
  vector<vector<int> > yRret(nRegions);
  vector<vector<int> > cFret(nRegions);
  vector<vector<int> > cRret(nRegions);

  int idxRet = 0;
  int regLen = 0;
  for(int idxReg = 0; idxReg < regions.size(); ++idxReg){
    regLen = max( regions.at(idxReg).regReadsF.at(regions.at(idxReg).regReadsF.size()-1),
                  regions.at(idxReg).regReadsR.at(regions.at(idxReg).regReadsR.size()-1))
           - min(regions.at(idxReg).regReadsF.at(0),regions.at(idxReg).regReadsR.at(0));
    if(regions.at(idxReg).regReadsF.size() < minReadsRegion || regions.at(idxReg).regReadsR.size() < minReadsRegion || regLen < minLregion){
    } else {
      yFret.at(idxRet) = regions.at(idxReg).regReadsF;
      yRret.at(idxRet) = regions.at(idxReg).regReadsR;
      cFret.at(idxRet) = regions.at(idxReg).regContF;
      cRret.at(idxRet) = regions.at(idxReg).regContR;
      idxRet++;
    }
  }
  if(verbose>1){
    cout<<nRegions-idxRet<<" regions filtered out."<<endl; 
  }
  yFret.resize(idxRet);
  yRret.resize(idxRet);
  cFret.resize(idxRet);
  cRret.resize(idxRet);

  Rcpp::List ret = Rcpp::List::create( //return the forward and reverse reads
      Rcpp::Named("yF", yFret),
      Rcpp::Named("yR", yRret),
      Rcpp::Named("cF", cFret),
      Rcpp::Named("cR", cRret)
      );

  return(ret);
}

/* Fill the scoreF/R vectors
*   Count the number of reads for each pos of the sliding window
*   F reads on the left within width and R on the right within width
*   Total length of the window = 2*width
*/
void getWindowsScores(vector<window>& windows, vector<int>& yF, vector<int>& yR, int width, int minDist){
  int i; //Read index iterator
  int startF=0, startR=0; //The read index for the next window
  int nbF=0, nbR=0; //Number of reads for one window
  int nF=yF.size(), nR=yR.size();

  std::vector<window>::iterator it;
  for(it = windows.begin(); it != windows.end(); ++it){
    nbF=0;
    i=startF;
    while((i<nF) && ((it->center-yF.at(i))>width)){
      i++;
    }
    startF=i;
    while((i<nF) && ((it->center-yF.at(i))<=width) && ((it->center-yF.at(i))>=minDist)){ //Fwd reads in the winow
      nbF++;
      i++;
    }
    nbR=0;
    i=startR;
    while((i<nR) && ((yR.at(i)-it->center)<minDist)){  //Rev reads that are out of the window on the left (before the window center)
      i++;
    }
    startR=i;
    while((i<nR) && ((yR.at(i)-it->center)<=width) && ((yR.at(i)-it->center)>=minDist)){ //Rwd reads in the winow
      nbR++;
      i++;
    }
    it->scoreF=nbF;
    it->scoreR=nbR;
  }
}

/* Fill ScoreRegionF/R vectors
*   Returns the number of regions
*   Merge the close windows and get the maximum of scoreF/R in these window
*/
int callRegions(vector<window>& windows, int width,int minReads, vector<region>& regions){
  int winIdx=0, nextWin=0, curWin=0; //iterators
  int maxF=0, maxR=0;
  int nRegions=0;
  int nWindows=windows.size();

  while(winIdx < nWindows){
    R_CheckUserInterrupt(); //In long loops 
    if((windows.at(winIdx).scoreF >= minReads) & (windows.at(winIdx).scoreR >= minReads)){
      maxF=windows.at(winIdx).scoreF;
      maxR=windows.at(winIdx).scoreR;
      regions.at(nRegions).start=windows.at(winIdx).center-width;

      nextWin=winIdx+1;
      curWin=winIdx;
      while((nextWin < nWindows) && ((windows.at(nextWin).center-windows.at(curWin).center)<=width*2)){
        if((windows.at(nextWin).scoreF >= minReads) && (windows.at(nextWin).scoreR >= minReads)){
          if(windows.at(nextWin).scoreF > maxF){
            maxF=windows.at(nextWin).scoreF;
          }
          if(windows.at(nextWin).scoreR > maxR){
            maxF=windows.at(nextWin).scoreR;
          }
          curWin=nextWin;
        }
        nextWin++;
      }
      regions.at(nRegions).scoreRegionF=maxF;
      regions.at(nRegions).scoreRegionR=maxR;
      regions.at(nRegions).end=windows.at(curWin).center+width;
      nRegions++;
      winIdx=nextWin;
      
    } else{
      winIdx++;
    }
  }
  return(nRegions);
}

int callLongRegions(vector<window>& windows, int width, int minReads, vector<region>& regions, int maxStep, int kStep, int minLregion){
  int pp=0, winIdx=0, nextWin=0, curWin=0, start=0, minScore=0, minIdx=0, cutCenter=0;
  //start    : the index of first center in current segment
  //minScore : save the minimium score of current region
  //min      : save the index of center correponding to the min score, used for cut window in ceter
  //cutCenter: a flag indicating if last region was cut inthe center of window
  //kStep    : width/step, when cut in the center, I want the index of start window is cut point+kStep
  int nRegions=0;
  int nWindows=windows.size();
  
  int cc = 0;
  cout<<"max regions"<<regions.size()<<endl;
  
  while(winIdx < nWindows){
    R_CheckUserInterrupt();
    if(((windows.at(winIdx).scoreF >= minReads) && (windows.at(winIdx).scoreR >= minReads)) || (cutCenter==1)){
      nRegions++;
      if(cutCenter == 0){
//        regions.at(nRegions).start=windows.at(winIdx).center-width;
        regions.at(nRegions-1).start=windows.at(winIdx).center-width;
        start = winIdx;
        minIdx = start;
        minScore = min(windows.at(winIdx).scoreF, windows.at(winIdx).scoreR);
      } else{
//        regions.at(nRegions).start = regions.at(nRegions-1).end+1; //isn't that shitty?
//        cout<<"here:"<<nRegions-1<<endl;
//        cout<<"start "<<" = "<<regions.at(nRegions-1).start<<endl;
//        cout<<"end   "<<nRegions-2<<" = "<<regions.at(nRegions-2).end+1<<endl;
        regions.at(nRegions-1).start = regions.at(nRegions-2).end+1;
        start = minIdx+kStep;
        minIdx = start;
        minScore = min(windows.at(start).scoreF, windows.at(start).scoreR);
        for(pp=start; pp<=winIdx; pp++){
          if(windows.at(pp).scoreF < minScore){
            minScore = windows.at(pp).scoreF;
            minIdx = pp;
          }
          if(windows.at(pp).scoreR < minScore){
            minScore = windows.at(pp).scoreR;
            minIdx=pp;
          }
          pp++; //Why?
        }
      }
      nextWin=winIdx+1; //w=p+1
      curWin=winIdx;    //ww=p
      while(((nextWin-start)<=maxStep) && ((windows.at(nextWin).center-windows.at(curWin).center) <= width*2) && (nextWin < nWindows)){
        if((windows.at(nextWin).scoreF>=minReads) && (windows.at(nextWin).scoreR>=minReads)){
          curWin = nextWin;
          if(windows.at(nextWin).scoreF < minScore){
            minScore = windows.at(nextWin).scoreF;
            minIdx = nextWin;
          }
          if(windows.at(nextWin).scoreR < minScore){
            minScore = windows.at(nextWin).scoreR;
            minIdx = nextWin;
          }
        }
        nextWin++;
      }
      if(nextWin == nWindows){
//        regions.at(nRegions).end = windows.at(curWin).center+width;
        regions.at(nRegions-1).end = windows.at(curWin).center+width;
      } else if((nextWin-start) > maxStep){
//      } else if((curWin-start) >= maxStep){
//        regions.at(nRegions).end = windows.at(minIdx).center;
        regions.at(nRegions-1).end = windows.at(minIdx).center;
        cutCenter = 1;
        cc++;
      } else{
//        regions.at(nRegions).end = windows.at(curWin).center+width;
        regions.at(nRegions-1).end = windows.at(curWin).center+width;
        cutCenter = 0;
      }
      if((regions.at(nRegions-1).end - regions.at(nRegions-1).start) < minLregion){
        nRegions--;
      }
//      nRegions++;
      winIdx=nextWin;
    } else{
      winIdx++;
    }
//    cout<<winIdx<<endl;
  }
  cout<<"p="<<winIdx<<"; pp="<<pp<<endl;
  cout<<"cc = "<<cc<<endl;
  return(nRegions);
}


/* Find which reads belongs to which region
*  TODO: The reads will change when I add map & control
*/
void segR(vector<int>& yF, vector<int>& yR, vector<int>& contF, vector<int>& contR, vector<window>& windows, vector<region>& regions){
  int minLoc, maxLoc; //tmp save regions boundaries
  int nF=yF.size(), nR=yR.size(), ncF=contF.size(), ncR=contR.size(); //TODO: prolly useless <-
  int nRegions=regions.size();
  
  int yStartF=0, yEndF=0;
  int yStartR=0, yEndR=0;
  int cStartF=0, cEndF=0;
  int cStartR=0, cEndR=0;

  int pF=0, pR=0, pcF=0, pcR=0, indStart=0, indEnd=0, indStartF=0, indStartR=0, indEndF=0, indEndR=0, temp=0;
  int j=0;
  
  vector<region>::iterator it; //vector<class>::iterator
  for(it = regions.begin(); it != regions.end(); ++it){
    minLoc = it->start;
    maxLoc = it->end;

    while((pF<nF) && (yF[pF]<it->start)){
      pF++;
    }
    indStartF=pF;
    temp=minLoc;
    if(temp>it->start){
      minLoc = min(temp, yF[indStartF]);
    } else{
      minLoc = max(temp, yF[indStartF]);
    }
    while((pR<nR) && (yR[pR]<minLoc)){ //idx of first yR bound by minLoc
      pR++;
    }
    indStartR=pR;
    while((pR<nR) && (yR[pR]<=it->end)){ //idx of last yR bound by region end
      pR++;
    }
    indEndR=min(pR-1, nR);
    indEndR=max(indEndR, indStartR);
    temp=maxLoc;
    if(temp<=it->end){
      maxLoc = min(temp, yR[indEndR]);
    } else{
      maxLoc = max(temp, yR[indEndR]);
    }
    while((pF<nF) && (yF[pF]<=maxLoc)){
      pF++;
    }
    indEndF=min(pF-1, nF);
    indEndF=max(indEndF, indStartF);

    //Controls
    if((ncF>0) && (ncR>0)){
      //cF
      while((pcF<ncF) && (contF[pcF]<minLoc)){
        pcF++;
      }
      indStart=pcF;
      while((pcF<ncF) && (contF[pcF]<=maxLoc)){
        pcF++;
      }
      indEnd=min(pcF, ncF);
      vector<int> regionCF(indEnd-indStart);
      std::copy(contF.begin()+indStart, contF.begin()+indEnd, regionCF.begin());
      it->regContF = regionCF;
      //cR
      while((pcR<ncR) && (contR[pcR]<minLoc)){
        pcR++;
      }
      indStart=pcR;
      while((pcR<ncR) && (contR[pcR]<=maxLoc)){
        pcR++;
      }
      indEnd=min(pcR, ncR);
      vector<int> regionCR(indEnd-indStart);
      std::copy(contR.begin()+indStart, contR.begin()+indEnd, regionCR.begin());
      it->regContR = regionCR;
    } else{
      it->regContF = vector<int> (0);
      it->regContR = vector<int> (0);
    }

    //Create the vectors
    vector<int> regionYF(indEndF-indStartF+1);
    vector<int> regionYR(indEndR-indStartR+1);
    std::copy(yF.begin()+indStartF, yF.begin()+indEndF+1, regionYF.begin());
    std::copy(yR.begin()+indStartR, yR.begin()+indEndR+1, regionYR.begin());
    it->regReadsF = regionYF;
    it->regReadsR = regionYR;
    //cout<<"minLoc = "<<minLoc<<", maxLoc = "<<maxLoc<<endl;
  }
}    
