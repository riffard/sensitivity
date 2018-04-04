#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>

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

using namespace std;

TH1D* hTot = NULL;

void AddComponant(string filename, string label, THStack* hs, TLegend* leg, int color = 10, double scale=1){
  
  TFile *f = TFile::Open(filename.c_str());
  TH1D* hE = (TH1D*)f->Get("E_keVee");

  hE->SetTitle(label.c_str());

  hE->SetLineColor(color);
  hE->SetLineWidth(2);


  
  for(int i=1; i<=hE->GetNbinsX(); ++i)hE->SetBinError(i,0);

  if(!hTot) hTot = (TH1D*)hE->Clone();
  else hTot->Add(hE);

  hs->Add(hE);

  leg->AddEntry(hE, label.c_str(), "l");
  
}

void plot_background_model(){
  
  THStack* hs = new THStack;
  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  
  map<string, string> componants;
  componants["WIMP #sigma_{#chi-N} = 10^{-45}cm^{2}"] = "../generate_pdf/pdf/wimp_100_GeV_100_Vcm.root";
  componants["Neutrino NR"] = "../generate_pdf/pdf/neutrino_NR_100_Vcm.root";
  componants["Rn220"] = "../generate_pdf/pdf/Rn220_ER_100Vcm.root";
  componants["Rn222"] = "../generate_pdf/pdf/Rn222_ER_100Vcm.root";


  map<string, int> colors;
  colors["WIMP"] = 15;
  colors["Neutrino NR"] = 10;
  colors["Rn220"] = 11;
  colors["Rn222"] = 12;

  
  for(auto it: componants){
    
    AddComponant(it.second, it.first, hs, leg, colors[it.first]);
  }
  

  hTot->SetLineColor(1);

  TCanvas* c = new TCanvas;
  hs->Draw("nostack");

  hs->GetXaxis()->SetTitle("Recoil energy [keVee]");
  hs->GetYaxis()->SetTitle("Rate [evt/t/y/keV]");
  
  hTot->Draw("same");
  leg->Draw("same");

  c->SetLogy();
  
} 
