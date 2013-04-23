/**
* TODO:
*   callRegionsL
*   add mappability profile and control reads. Used only in segR
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
  vector<int> readsF;
  vector<int> readsR;
};

SEXP seg_chr(vector<int>& yF, vector<int>& yR, vector<int>& contF, vector<int>& contR, int minReads, int step, int width, int maxStep, int minLength, int minDist, int verbose);
void getWindowsScores(vector<window>& windows, vector<int>& yF, vector<int>& yR, int width, int minDist);
int callRegions(vector<window>& windows, int width, int minReads, vector<region>& regions);
void segR(vector<int> yF, vector<int>yR, vector<window>& windows, vector<region>& regions);

//[[Rcpp::export]]
SEXP CPP_seg(Rcpp::List data, Rcpp::List dataC, int minReads, int step, int width, int maxStep, int minLength, int minDist, int verbose){
//minReads : The min number of reads for a window to be included in a region
//step     : The step for the sliding window
//width    : Half the length of the sliding window
//maxStep  : Used in callRegionsL
//minLength: Used in callRegionsL
//minDist  : in PING, when reads are too close to the center of the windows, they are not counted in the window score
//verbose  : Print debug values
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
    ret[i] = seg_chr(yF, yR, yF, yR, minReads, step, width, maxStep, minLength, minDist, verbose);
  }
  ret.names() = chrs;
  return(Rcpp::wrap(ret));
}


SEXP seg_chr(vector<int>& yF, vector<int>& yR, vector<int>& contF, vector<int>& contR, int minReads, int step, int width, int maxStep, int minLength, int minDist, int verbose){

  std::sort(yF.begin(), yF.end());
  std::sort(yR.begin(), yR.end());
  if(length(contF)>0 && length(contR)>0){
    sort(contF.begin(), contF.end());
    sort(contR.begin(), contR.end());
  }
  
  int m=std::min(yF.at(0), yR.at(0));
  int M=std::max(yF.at(yF.size()-1), yR.at(yR.size()-1));

  int nWindows = (M-m)/step; //Number of windows

  vector<window> windows(nWindows);
  int i;
  for(i=0; i<nWindows; i++){
    windows.at(i).center=m+i*step;
  }

  //Vector of the maximal possible size for the bigger regions (merged adjacent windows)
  std::vector<region> regions((M-m)/(2*width));

  getWindowsScores(windows, yF, yR, width, minDist); //Sliding window
  int nRegions;
  if(maxStep==0){
    nRegions = callRegions(windows, width, minReads, regions); //Get the actual regions
  } else{
    //int nRegions = callLongRegions(windows, width, minReads, regions);
  }
    
  regions.resize(nRegions);

  if(verbose>1){
    cout<<"nWindows : "<<nWindows<<endl;
    cout<<"nReg: "<<nRegions<<endl;
  }

  if(nRegions>0){
    segR(yF, yR, windows, regions);
  }
  else{
  }

  vector<vector<int> > Freads(nRegions); //vector of vectors
  vector<vector<int> > Rreads(nRegions);

  for(int idxRet = 0; idxRet < regions.size(); ++idxRet){
      Freads.at(idxRet) = regions.at(idxRet).readsF;
      Rreads.at(idxRet) = regions.at(idxRet).readsR;
  }

  Rcpp::List ret = Rcpp::List::create( //return the forward and reverse reads
      Rcpp::Named("Freads", Freads),
      Rcpp::Named("Rreads", Rreads)
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
  int i=0, winIdx=0, nextWin=0, curWin=0; //iterators
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
      
    }
    else{
      winIdx++;
    }
  }
  return(nRegions);
}



/* Find which reads belongs to which region
*  TODO: The reads will change when I add map & control
*/
void segR(vector<int> yF, vector<int>yR, vector<window>& windows, vector<region>& regions){
  int minLoc, maxLoc; //tmp save regions boundaries
  int nF=yF.size(), nR=yR.size();
  int nRegions=regions.size();
  
  int startF=0, endF=0;
  int startR=0, endR=0;

  vector<region>::iterator it; //vector<class>::iterator
  for(it = regions.begin(); it != regions.end(); ++it){
    minLoc = it->start; 
    maxLoc = it->end;
    
    //subset yF & yR based on the start & end of the regions
    while(yF[startF]<minLoc){
      startF++;
    }
    endF=startF;
    while(yF[endF]<maxLoc){
      endF++;
    }
    endF--;
    while(yR[startR]<minLoc){
      startR++;
    }
    endR=startR;
    while(yR[endR]<maxLoc){
      endR++;
    }
    endR--;

    vector<int> regionF(endF-startF);
    vector<int> regionR(endR-startR);
    std::copy(yF.begin()+startF, yF.begin()+endF, regionF.begin());
    std::copy(yR.begin()+startR, yR.begin()+endR, regionR.begin());

    it->readsF = regionF;
    it->readsR = regionR;

    //Regions do not overlap
    startF = endF;
    startR = endR;

    //process control data
    if(contF.size()>0 && contR.size()>0){

    }
    
  }
}
