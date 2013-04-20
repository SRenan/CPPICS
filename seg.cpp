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

void getWindowsScores(vector<window>& windows, vector<int>& yyF, vector<int>& yyR, int width, int minDist);
int callRegions(vector<window>& windows, int width, int minReads, vector<region>& regions);
void segR(vector<int> yyF, vector<int>yyR, vector<window>& windows, vector<region>& regions);



//[[Rcpp::export]]
SEXP CPP_seg(SEXP yF, SEXP yR, int minReads, int step, int width, int maxStep, int minLength, int minDist){ //yF=dataF
//minReads : The min number of reads for a window to be included in a region
//step     : The length of the sliding window
//width    : The max length of the regions
//maxStep  : Used in callRegionsL
//minLength: Used in callRegionsL
//minDist  : in PING, when reads are too close to the center of the windows, they are not counted in the window score

  int i;
  int verbose = 2;

  //To c++ vectors
  std::vector<int> yyF = Rcpp::as<std::vector<int> >(yF);
  std::vector<int> yyR = Rcpp::as<std::vector<int> >(yR);
  
  int nF = yyF.size();
  int nR = yyR.size();

  //sort yyF & yyR
  std::sort(yyF.begin(), yyF.end());
  std::sort(yyR.begin(), yyR.end());
  int m=std::min(yyF.at(0), yyR.at(0));
  int M=std::max(yyF.at(nF-1), yyR.at(nR-1));

  int nWindows = (M-m)/step; //Number of windows

  vector<window> windows(nWindows);
  for(i=0; i<nWindows; i++){
    windows.at(i).center=m+i*step;
  }

  //Vector for the bigger regions (merged adjacent windows)
  std::vector<region> regions((M-m)/(2*width));

  getWindowsScores(windows, yyF, yyR, width, minDist); //Sliding window
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
    segR(yyF, yyR, windows, regions);
  }
  else{
  }

  vector<vector<int> > Freads(nRegions);
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
void getWindowsScores(vector<window>& windows, vector<int>& yyF, vector<int>& yyR, int width, int minDist){
  int i; //Read index iterator
  int startF=0, startR=0; //The read index for the next window
  int nbF=0, nbR=0; //Number of reads for one window
  int nF=yyF.size(), nR=yyR.size();

  std::vector<window>::iterator it;
  for(it = windows.begin(); it != windows.end(); ++it){
    nbF=0;
    i=startF;
    while((i<nF) && ((it->center-yyF.at(i))>width)){
      i++;
    }
    startF=i;
    while((i<nF) && ((it->center-yyF.at(i))<=width) && ((it->center-yyF.at(i))>=minDist)){ //Fwd reads in the winow
      nbF++;
      i++;
    }
    nbR=0;
    i=startR;
    while((i<nR) && ((yyR.at(i)-it->center)<minDist)){  //Rev reads that are out of the window on the left (before the window center)
      i++;
    }
    startR=i;
    while((i<nR) && ((yyR.at(i)-it->center)<=width) && ((yyR.at(i)-it->center)>=minDist)){ //Rwd reads in the winow
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
void segR(vector<int> yyF, vector<int>yyR, vector<window>& windows, vector<region>& regions){
  int minLoc, maxLoc; //tmp save regions boundaries
  int nF=yyF.size(), nR=yyR.size();
  int nRegions=regions.size();
  
  int startF=0, endF=0;
  int startR=0, endR=0;

  vector<region>::iterator it; //vector<class>::iterator
  for(it = regions.begin(); it != regions.end(); ++it){
    minLoc = it->start; 
    maxLoc = it->end;
    
    //subset yyF & yyR based on the start & end of the regions
    while(yyF[startF]<minLoc){
      startF++;
    }
    endF=startF;
    while(yyF[endF]<maxLoc){
      endF++;
    }
    endF--;
    while(yyR[startR]<minLoc){
      startR++;
    }
    endR=startR;
    while(yyR[endR]<maxLoc){
      endR++;
    }
    endR--;

    vector<int> regionF(endF-startF);
    vector<int> regionR(endR-startR);
    std::copy(yyF.begin()+startF, yyF.begin()+endF, regionF.begin());
    std::copy(yyR.begin()+startR, yyR.begin()+endR, regionR.begin());

    it->readsF = regionF;
    it->readsR = regionR;

    //Regions do not overlap
    startF = endF;
    startR = endR;
  }
}
