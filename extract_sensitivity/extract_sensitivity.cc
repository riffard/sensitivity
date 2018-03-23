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

#include "Math/PdfFuncMathCore.h"

//----------------------------------------------------
//----------------------------------------------------
#include "GlobalParameters.hh"
#include "ModelBuilder.hh"
#include "Style.hh"

//----------------------------------------------------
//----------------------------------------------------

using namespace std;

//----------------------------------------------------
//----------------------------------------------------
int main(int argc,char** argv){
  LoadStyle();
  //----------------------------------------------------
  //----------------------------------------------------
  TApplication theApp("App", &argc, argv);

  double TotalExpousre = 1; // ton.year
  int N_pseudo_exp = 400;
  int vervose_level = -1;
  //vector<double> sigma_test = {0, 5, 10, 20, 25, 30};
  vector<double> sigma_test = {0, 5, 10, 20, 22, 25, 30, 40};
  //vector<double> sigma_test = {0, 5};

  vector<vector<double>> q_collection;
  
  vector<TH1D*> hq;
  for(size_t i=0; i< sigma_test.size(); ++i){
    hq.push_back(new TH1D("","",100,-1,10));

    q_collection.push_back(vector<double> ());

  }
  TH1D* hmu_estimator = new TH1D("","",10000,-10, 40);
  
  
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

  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);
  
  cout<<endl<<"Start pseudo experiments"<<endl<<endl;;
  
  for(int n =0; n <N_pseudo_exp; ++n){

    cout << n+1 << " / " << N_pseudo_exp << endl;
       
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

    hmu_estimator->Fill(background->GetAmplitude("wimp")->getVal());

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
      

      hq[i]->Fill(q);
      q_collection[i].push_back(q);
    }
    //----------------------------------------------------
    //----------------------------------------------------

  }

  TGraph* gr_chi2 = new TGraph();
  gr_chi2->SetLineColor(12);
  gr_chi2->SetLineWidth(3);
  for(int b=1; b<=hq[0]->GetNbinsX(); b++){
    double x = hq[0]->GetBinCenter(b);
    gr_chi2->SetPoint(gr_chi2->GetN(), x , ROOT::Math::chisquared_pdf(x,1));
  }
  
  
  
  vector<double> median_collection;
  for(size_t i = 0 ; i< sigma_test.size(); i++){
    median_collection.push_back(TMath::Median(q_collection[i].size(), &q_collection[i][0]));
    hq[i]->Scale(1/hq[i]->Integral("width"));

    for(int b=1; b<=hq[i]->GetNbinsX(); b++) hq[i]->SetBinError(b, 0);
    
  }
      
  hq[0]->SetLineWidth(3);
  hq[0]->SetLineColor(10);
  hq[0]->SetFillColorAlpha(10, 0.5);

  TLine* line = new TLine();
  line->SetLineWidth(3);

  TGraph* gr_pvalue = new TGraph();
  
  for(size_t i =1; i< sigma_test.size(); ++i){


    int bin_median = hq[0]->GetXaxis()->FindBin(median_collection[i]);
    gr_pvalue->SetPoint(gr_pvalue->GetN(), sigma_test[i], hq[0]->Integral(bin_median, hq[0]->GetNbinsX(), "width"));

    
    hq[i]->SetLineWidth(3);
    hq[i]->SetLineColor(11);
    hq[i]->SetFillColorAlpha(11, 0.5);

    
    THStack* hs = new THStack;
    hs->Add(hq[i]);
    hs->Add(hq[0]);
    
    TCanvas* c = new TCanvas;
    hs->Draw("nostack");
    gr_chi2->Draw("samel");
    hs->SetTitle("");
    c->SetLogy();

    line->DrawLine(median_collection[i], 0, median_collection[i], hq[i]->GetMaximum());

    c->Update();
  }
  
  //hq->Draw();

  TCanvas* c = new TCanvas;
  hmu_estimator->Draw();
  c->Update();
  
  
  gr_pvalue->SetMarkerStyle(21);
  c = new TCanvas;
  gr_pvalue->Draw("apl");
  c->Update();
  
  
  cout<<"Done "<<endl;
  theApp.Run();

  
} 
