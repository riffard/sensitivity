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
#include "RooDataSet.h"
#include "RooFitResult.h"

//----------------------------------------------------
//----------------------------------------------------
#include "GlobalParameters.hh"
#include "ModelBuilder.hh"
//----------------------------------------------------
//----------------------------------------------------

using namespace std;

//----------------------------------------------------
//----------------------------------------------------
int main(int argc,char** argv){
  
  //----------------------------------------------------
  //----------------------------------------------------
  TApplication theApp("App", &argc, argv);

  double TotalExpousre = 1; // ton.year
  int N_pseudo_exp = 100;
  int vervose_level = -1;
  vector<double> sigma_test = {0.01, 0.1, 1};

  
  vector<TH1D*> hq;
  for(size_t i=0; i< sigma_test.size(); ++i){
    hq.push_back(new TH1D("","",1000,-1,10));
  }

  
  //----------------------------------------------------
  // Global parameters
  //----------------------------------------------------  
  GlobalParameters* parameters = GlobalParameters::GetInstance();

  //----------------------------------------------------
  // Set the observables variables
  //----------------------------------------------------
  RooRealVar E_recoil("E_recoil","E_recoil",parameters->E_min, 100) ;
  RooRealVar s1("s1", "s1 [npe]", 0, parameters->s1_max) ;
  RooRealVar s2("s2", "s2", 0, parameters->s2_max) ;
  RooRealVar logs2s1("logs2s1", "logs2s1", parameters->logs2s1_min, parameters->logs2s1_max) ;

  map<string, RooRealVar* > variables;
  variables["E_recoil"] = &E_recoil;
  variables["s1"]       = &s1;
  variables["s2"]       = &s2;
  variables["logs2s1"]  = &logs2s1;

  string analysis_variables = "s1,logs1s2";

  //ModelBuilder* signal = new ModelBuilder("Signal", analysis_variables, &variables);
  //signal->AddComponant("wimp","../generate_pdf/pdf/wimp_100_GeV_100_Vcm.root", "s1_logs1s2", RooArgList(s1, logs2s1));  
  //signal->BuildPdf();
  //RooAddPdf* sinal_pdf = signal->GetPdf();

  ModelBuilder* background = new ModelBuilder("Background", analysis_variables,  &variables);
  background->AddComponant("NR_neutrino","../generate_pdf/pdf/neutrino_NR_100_Vcm.root", "s1_logs1s2", RooArgList(s1, logs2s1));
  background->AddComponant("ER_Flat", "../generate_pdf/pdf/flat_ER_100Vcm.root", "s1_logs1s2", RooArgList(s1, logs2s1));
  background->AddComponant("wimp","../generate_pdf/pdf/wimp_100_GeV_100_Vcm.root", "s1_logs1s2", RooArgList(s1, logs2s1));
  background->BuildPdf();

  RooAddPdf* background_pdf = background->GetPdf();


  
  
  for(int n =0; n <N_pseudo_exp; ++n){
       
    //----------------------------------------------------
    // Generate data without WIMP
    //----------------------------------------------------
    background->GetAmplitude("wimp")->setVal(0);
    RooDataSet *dataBkgOnly = background->GenerateFakeData(TotalExpousre, RooArgList(s1, logs2s1));

    //----------------------------------------------------
    // Fit un constrain Likelyhood
    //----------------------------------------------------    
    background->GetAmplitude("wimp")->setConstant(kFALSE);
    RooFitResult* fit_unconstrain = background_pdf->fitTo(*dataBkgOnly,
							  RooFit::Extended(kTRUE),
							  RooFit::PrintLevel(vervose_level),
							  RooFit::Save());

    double lnL_unconstrain = fit_unconstrain->minNll();
    
    //----------------------------------------------------
    // Loop over all the test cross sections
    //----------------------------------------------------
    for(size_t i = 0 ; i< sigma_test.size(); i++){
      
      background->GetAmplitude("wimp")->setVal(sigma_test[i]);   
      background->GetAmplitude("wimp")->setConstant(kTRUE);
      RooFitResult* fit_constrain = background_pdf->fitTo(*dataBkgOnly,
							  RooFit::Extended(kTRUE),
							  RooFit::PrintLevel(vervose_level),
							  RooFit::Save());
      
      double lnL_constrain = fit_constrain->minNll();
      

      double q = 2*(lnL_constrain - lnL_unconstrain);
      
      cout << n+1 << " / " << N_pseudo_exp << " --> " << q << endl;
      //background->GetAmplitude("wimp")->Print();

    
      hq[i]->Fill(q);

    }
    

    /*
    if(n==100000){
      
      TH2D* hh_pdf = (TH2D*)data->createHistogram("s1,logs2s1",500,500) ;
      
      TCanvas* c2 = new TCanvas;
      hh_pdf->Draw("colz") ;
  
      c2->Update();
      
      
      TCanvas* c3 = new TCanvas;
      RooPlot* mesframe = s1.frame() ;
      data->plotOn(mesframe) ;
      background_pdf->plotOn(mesframe) ;
      background_pdf->plotOn(mesframe,
			     RooFit::Components( *(background->GetExtendPdf("ER_Flat"))),
			     RooFit::LineStyle(kDashed));
      
      background_pdf->plotOn(mesframe,
			     RooFit::Components( *(background->GetExtendPdf("NR_neutrino"))),
			     RooFit::LineStyle(kDashed));
      
      background_pdf->plotOn(mesframe,
			     RooFit::Components( *(background->GetExtendPdf("wimp"))),
			     RooFit::LineStyle(kDashed),
			     RooFit::LineColor(kRed));
      mesframe->Draw() ;
      c3->Update();
      }*/
  }
  

  
  TCanvas* c = new TCanvas;
  c->Divide(sigma_test.size(),1);
  for(size_t i =0; i< sigma_test.size(); ++i){
    c->cd(i+1);
    hq[i]->Draw();
  }
  
  //hq->Draw();
  c->Update();

  

  
  
  cout<<"Done "<<endl;
  theApp.Run();

  
} 
