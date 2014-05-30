/**
* TODO:
*   verbose in the args
*   maybe add a few header comments for maintainer
*/
#include <Rcpp.h>
#include <iostream>

using namespace std;

class pics_param{
  private:
  //paraEM
  int minK;
  int maxK;
  int tol;
  int B;
  std::string mSelect;
  bool mergePeaks;
  bool mapCorrect;
  //paraPrior
  int xi;     //avg DNA fragment size
  int rho;    //DNA fragment size distrib variance
  int alpha;  //First hyperparameter of the inverse Gamma distribution for sigma^2 in
  int beta;   //Sexond hyperparameter of the inverse Gamma distribution for sigma^2 in
  int lambda; //Precision of the prior for mu (used for histone data)
  int dMu;    //The distance between two neighbouring nucleosomes (should be 0 for TF)
  public:
  pics_param(Rcpp::List paraEM, Rcpp::List paraPrior); //Constructor prototype
  void display();
};

SEXP fitModel_region(vector<int>& yF, vector<int>& yR, vector<int>& cF, vector<int>& cR, const pics_param& para);

pics_param::pics_param(Rcpp::List paraEM, Rcpp::List paraPrior){
  //paraEM
  minK = paraEM["minK"];
  maxK = paraEM["maxK"];
  tol  = paraEM["tol"];
  B    = paraEM["B"];
  std::string str = paraEM["mSelect"]; //R character to std::string. can't pass it directly
  mSelect = str; 
  mergePeaks = paraEM["mergePeaks"];
  mapCorrect = paraEM["mapCorrect"];
  //paraPrior
  xi     = paraPrior["xi"];
  rho    = paraPrior["rho"];
  alpha  = paraPrior["alpha"];
  beta   = paraPrior["beta"];
  lambda = paraPrior["lambda"];
  dMu    = paraPrior["dMu"];
}

void pics_param::display(){
  //paraEM
  cout<<"EM parameters:"<<endl;
  cout<<" minK = "<<minK<<endl;
  cout<<" maxK = "<<minK<<endl;
  cout<<" tol  = "<<tol<<endl;
  cout<<" B    = "<<B<<endl;
  cout<<" mSelect    = "<<mSelect<<endl;
  cout<<" mergePeaks = "<<mergePeaks<<endl;
  cout<<" mapCorrect = "<<mapCorrect<<endl;
  //paraPrior
  cout<<"prior parameters:"<<endl;
  cout<<" xi    = "<<xi<<endl;
  cout<<" rho   = "<<rho<<endl;
  cout<<" alpha = "<<alpha<<endl;
  cout<<" beta  = "<<beta<<endl;
  cout<<" lambda= "<<lambda<<endl;
  cout<<" dMu   = "<<dMu<<endl;
}
  

//[[Rcpp::export]]
SEXP cpp_pics(Rcpp::List bigAssList, Rcpp::List paraEM, Rcpp::List paraPrior, int minReads){
  int verbose = 2;
  pics_param para (paraEM, paraPrior);
  if(verbose>1){
    para.display();
    cout<<"nRegion = "<<bigAssList.size()<<endl;
  }

  Rcpp::List ret(bigAssList.size());
  for(int regIdx = 0; regIdx<bigAssList.size(); regIdx++){
    Rcpp::List cur_region = bigAssList[regIdx];
    vector<int> yF = Rcpp::as<vector<int> >(cur_region["yF"]);
    vector<int> yR = Rcpp::as<vector<int> >(cur_region["yR"]);
    vector<int> cF = Rcpp::as<vector<int> >(cur_region["cF"]);
    vector<int> cR = Rcpp::as<vector<int> >(cur_region["cR"]);
    ret[regIdx] = fitModel_region(yF, yR, cF, cR, para);
  }
  return(ret);
}

/*
* Fit the model for a single region
*/
SEXP fitModel_region(vector<int>& yF, vector<int>& yR, vector<int>& cF, vector<int>& cR, const pics_param& para){
  //Should be in pics object
  int rangeMax = max(yF.at(yF.size()-1), yR.at(yR.size()-1));
  int rangeMin = min(yF.at(0), yR.at(0));
  //
  int range = rangeMax - rangeMin;

  int maxKK = min(para->maxK, range/(para->xi+4.*sqrt(para->beta/para->alpha))+1);
  

  return(Rcpp::wrap(range));
}


