//----------------------------------------------------
// Standard C/C++
//----------------------------------------------------
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>

//----------------------------------------------------
// ROOT
//----------------------------------------------------
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "THStack.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TMinuit.h"
#include "TColor.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TApplication.h"

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooHistPdf.h"


#include "../GlobalParameters.hh"

#include "ModelBuilder.hh"




using namespace std;




int main(int argc,char** argv){
  
  //----------------------------------------------------
  //----------------------------------------------------
  TApplication theApp("App", &argc, argv);


  //----------------------------------------------------
  // Global parameters
  //----------------------------------------------------  
  GlobalParameters* parameters = GlobalParameters::GetInstance();

  //----------------------------------------------------
  // Set the observables variables
  //----------------------------------------------------
  RooRealVar E_recoil("E_recoil","E_recoil",parameters->E_min, 100) ;
  RooRealVar s1("s1", "s1", 0, parameters->s1_max) ;
  RooRealVar s2("s2", "s2", 0, parameters->s2_max) ;
  RooRealVar logs2s1("logs2s1", "logs2s1", parameters->logs2s1_min, parameters->logs2s1_max) ;

  map<string, RooRealVar* > variables;
  variables["E_recoil"] = &E_recoil;
  variables["s1"]       = &s1;
  variables["s2"]       = &s2;
  variables["logs2s1"]  = &logs2s1;
  

  ModelBuilder* signal = new ModelBuilder(&variables);
  signal->AddComponant("../generate_pdf/pdf/wimp_100.root","s1_s2_wimp_100", RooArgList(s1, s2));
  
  

  ModelBuilder* background = new ModelBuilder(&variables);
  


  


  
  

  TFile* f = TFile::Open("../generate_pdf/pdf/wimp_100.root");
  


  
 
  TH2D* histo = (TH2D*) f->Get("s1_s2_wimp_100");

  
  RooDataHist data("bdata","bdata",RooArgList(s1,s2), histo);

  RooHistPdf dataPdf("name","title",RooArgList(s1, s2), data);



  RooPlot* s1frame = s1.frame();
  data.plotOn(s1frame);
  s1frame->Draw();


  RooPlot* s2frame = s2.frame();
  data.plotOn(s2frame);
  s2frame->Draw();



  cout<<"Done "<<endl;
  //theApp.Run();

  
} 
