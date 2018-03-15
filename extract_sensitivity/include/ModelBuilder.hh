#ifndef ModelBuilder_hh
#define ModelBuilder_hh 1

#include <map>

#include "TFile.h"

#include "TH2D.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgList.h"
#include "RooAbsArg.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"

using namespace std;

class ModelBuilder{

public:


  ModelBuilder(string name, string analysis_variable, map<string, RooRealVar* >* variables);
  ~ModelBuilder();

  void AddComponant(string contrib_name, string file_name, string histo_name, RooArgList arglist);
  
  void BuildPdf();
  
  RooAddPdf* GetPdf(){return fFullPdf;};

  RooArgList* GetPdfList(){return fPdfList;};

  RooDataSet* GenerateFakeData(double exposure, RooArgList arglist);

  RooExtendPdf* GetExtendPdf(string name){return fExtendedPdf["extented_pdf_" + name];};

  RooRealVar* GetAmplitude(string name){return (*fVariablesCollection)["A_" + name];};
  
  
private:

  string fName;
  string fAnalysisVariables;
  map<string, RooRealVar* >* fVariablesCollection;
  
  map<string, TH2D*> fComponants;
  map<string, RooDataHist*> fComponants_data;
  map<string, RooHistPdf*> fComponants_pdf;
  map<string, RooExtendPdf*> fExtendedPdf;
  
  RooArgList* fPdfList;
  RooAddPdf* fFullPdf;

  double fTotalRate; 
  vector<double> fRates;


  double fExposure;
  RooDataSet *fFakeDataSet;
  
};




#endif
