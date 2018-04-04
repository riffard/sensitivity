#ifndef ModelBuilder_hh
#define ModelBuilder_hh 1

#include <map>

#include "TFile.h"

#include "TH2D.h"

#include "RooWorkspace.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgList.h"
#include "RooAbsArg.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

using namespace std;

class ModelBuilder{

public: 
  enum ComponantType{Signal, Background};

  
public:


  ModelBuilder(RooWorkspace *w, string name, double exposure, string histo_name, string variables, int verbose = 1);
  ~ModelBuilder();


  void AddComponant(ComponantType componantType, string contrib_name, string file_name);

  

  void BuildModel(bool IsNPActivated);
  
  void BuildData();
  
private:
  
  static string getStringFromArray(vector<string> v){
    string str = "";
    for(size_t i=0; i<v.size(); ++i){
      if(i > 0) str += ",";
      str += v[i];
    }
    return str;
  }

  
  
private:

  RooWorkspace *fw;
  string fName;
  double fExposure;
  int fVerbose;
  char *fHistoName;
  char *fVariables;
  
  vector<string> fModelComponant;
  vector<string> fNuisanceParameters;
  vector<string> fNuisanceParametersConstrain;
  vector<string> fNuisanceParametersGlobalObs;
  vector<string> fInterestParameters;


  
  string fAnalysisVariables;
  
  map<string, TH2D*> fComponants;
  map<string, RooDataHist*> fComponants_data;
  map<string, RooHistPdf*> fComponants_pdf;
  map<string, RooExtendPdf*> fExtendedPdf;

  RooArgList* fPdfList;
  RooAddPdf* fFullPdf;

  double fTotalRate; 
  //vector<double> fRates;
  map<string, double> fRates;
  

  RooDataSet *fFakeDataSet;


  
};




#endif
