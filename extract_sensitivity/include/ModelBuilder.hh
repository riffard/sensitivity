#ifndef ModelBuilder_hh
#define ModelBuilder_hh 1

#include <map>

#include "TFile.h"

#include "TH2D.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"

using namespace std;

class ModelBuilder{

public:


  ModelBuilder(map<string, RooRealVar* >* variables);
  ~ModelBuilder();

  void AddComponant(string file_name, string histo_name, RooArgList arglist);

  void BuildPdf();
  
private:

  map<string, RooRealVar* >* fVariables;

  map<string, TH2D*> fComponants;
  map<string, RooDataHist*> fComponants_data;
  map<string, RooHistPdf*> fComponants_pdf;
  

};




#endif
