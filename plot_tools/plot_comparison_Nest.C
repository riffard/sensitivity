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

void Get_bands(TH2D*h2, TGraph* gr_mean, TGraph* gr_sup, TGraph* gr_inf, TMultiGraph* mg, int color){

  h2->SetStats(0);
  
  for(int i=1; i<= h2->GetNbinsX(); ++i){

    double s1 = h2->GetXaxis()->GetBinCenter(i);
    if(s1 >120)break;

    TH1D* h_proj = h2->ProjectionY("h_proj",i,i);

    double mean = h_proj->GetMean();
    double rms = h_proj->GetRMS();
    
    gr_mean->SetPoint(i-1, s1, mean);
    gr_sup->SetPoint(i-1, s1, mean+rms);
    gr_inf->SetPoint(i-1, s1, mean-rms);
    
  }

  gr_mean->SetLineColor(10);
  gr_sup->SetLineColor(10);
  gr_inf->SetLineColor(10);

  gr_mean->SetLineWidth(3);
  gr_sup->SetLineWidth(3);
  gr_inf->SetLineWidth(3);

  gr_sup->SetLineStyle(7);
  gr_inf->SetLineStyle(7);

  
  TCanvas* c = new TCanvas;

  h2->DrawClone("colz");
  gr_mean->DrawClone("samepl");
  gr_sup->DrawClone("samepl");
  gr_inf->DrawClone("samepl");


  gr_mean->SetLineColor(color);
  gr_sup->SetLineColor(color);
  gr_inf->SetLineColor(color);
  
  mg->Add(gr_mean);
  mg->Add(gr_sup);
  mg->Add(gr_inf);
  

  
}



void plot_comparison_Nest(){


  TFile *fnest = TFile::Open("../generate_pdf/pdf/dd_gun_100_Vcm.root");
  TH2D* hs1s2_nest = (TH2D*)fnest->Get("s1_s2");
  hs1s2_nest->SetTitle("NEST");

  hs1s2_nest->RebinX(20);
  hs1s2_nest->RebinY(10);
  
  TFile *fdata = TFile::Open("NR_Band.root");
  TTree *data_tree = (TTree*)fdata->Get("data");


  TH2D*hs1s2_data = (TH2D*)hs1s2_nest->Clone("hs1s2_data");
  hs1s2_data->SetTitle("Data");
  data_tree->Draw("s2:s1>>hs1s2_data","field > 90 && field < 120", "goff");
  data_tree->SetMarkerColor(kRed);

  

  TGraph* gr_mean_data = new TGraph();
  TGraph* gr_sup_data = new TGraph();
  TGraph* gr_inf_data = new TGraph();

  TGraph* gr_mean_nest = new TGraph();
  TGraph* gr_sup_nest = new TGraph();
  TGraph* gr_inf_nest = new TGraph();


  TMultiGraph* mg = new TMultiGraph();
  
  Get_bands(hs1s2_data, gr_mean_data, gr_sup_data, gr_inf_data, mg, 10);
  
  Get_bands(hs1s2_nest, gr_mean_nest, gr_sup_nest, gr_inf_nest, mg, 1);


  TLegend* leg_color = new TLegend(0.15, 0.6,0.4,0.85);

  leg_color->AddEntry(gr_mean_data,"Data", "l");
  leg_color->AddEntry(gr_mean_nest,"NEST", "l");

  leg_color->AddEntry(gr_mean_nest,"mean", "l");
  leg_color->AddEntry(gr_sup_nest,"#pm1#sigma", "l");
  
  TCanvas* c = new TCanvas;

  hs1s2_nest->Draw("colz");
  hs1s2_data->Draw("same");
  
  hs1s2_nest->SetStats(0);
  hs1s2_nest->SetXTitle("s1");
  hs1s2_nest->SetYTitle("s2");
  hs1s2_data->SetMarkerColor(kRed);
  hs1s2_data->SetMarkerStyle(2);


  
  c = new TCanvas;
  mg->Draw("a");

  mg->GetXaxis()->SetTitle("s1");
  mg->GetYaxis()->SetTitle("s2");
  
  leg_color->Draw("same");
  
} 
