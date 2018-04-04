#ifndef HomeMadePLRData_hh
#define HomeMadePLRData_hh

#include<iostream>

#include "TH1D.h"
#include "TMath.h"
#include "TLine.h"
#include "TAxis.h"
#include "TPad.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCanvas.h"

using namespace std;

class HomeMadePLRData{

public:

  HomeMadePLRData(){
    hqH0 = NULL;
    hqH1 = NULL;
  }
  
  HomeMadePLRData(string name){
    Init(name); 
  }
  
  ~HomeMadePLRData(){
  }

  
  
  
  void AddData(double qH0, double qH1){
    
    hqH0->Fill(qH0);
    hqH1->Fill(qH1);

    qH0Collection.push_back(qH0);
    qH1Collection.push_back(qH1);
    
  }


  double GetPValue(){

    double median_qH1 = TMath::Median(qH1Collection.size(), &qH1Collection[0]);

    hqH0->Scale(1/hqH0->Integral("width"), "nosw2");
    hqH1->Scale(1/hqH1->Integral("width"), "nosw2");

    int bin_start = hqH0->GetXaxis()->FindBin(median_qH1);
    int bin_end = hqH0->GetNbinsX();;
    
    double pvalue = hqH0->Integral(bin_start, bin_end, "width");

    return pvalue;
  }
  
  void Draw(){

    hqH0->Scale(1/hqH0->Integral("width"), "nosw2");
    hqH1->Scale(1/hqH1->Integral("width"), "nosw2");
    
    double median_qH1 = TMath::Median(qH1Collection.size(), &qH1Collection[0]);

    hqH0->SetLineColor(10);
    hqH1->SetLineColor(11);
    hqH0->SetFillColorAlpha(10, 0.5);
    hqH1->SetFillColorAlpha(11, 0.5);
    
    THStack* hs = new THStack();
    TLegend* leg = new TLegend(0.55, 0.75, 0.85, 0.85);
    leg->SetTextSize(0.05);
    leg->SetNColumns(2);
    
    hs->Add(hqH0);
    leg->AddEntry(hqH0, "f(q|H_{0})","l");
    hs->Add(hqH1);
    leg->AddEntry(hqH1, "f(q|H_{1})","l");    

    TLine* line = new TLine();
    TCanvas* c = new TCanvas();
    c->SetLogy();
    hs->Draw("nostack");
    hs->GetXaxis()->SetTitle("t_{#mu} = -2#ln#lambda(#mu)");
    c->Update();
    
    
    double ymin = pow(10,gPad->GetUymin());
    double ymax = pow(10,gPad->GetUymax());
    line->DrawLine(median_qH1, ymin, median_qH1, ymax);

    leg->Draw("same");
  
    c->Update();
    
  }
  
private:
  
  void Init(string name){
    hqH0 = new TH1D(Form("hqH0_%s",name.c_str()),name.c_str(), 100,0,10);
    hqH1 = new TH1D(Form("hqH1_%s",name.c_str()),name.c_str(), 100,0,10);
  }
  
  TH1D *hqH0;
  TH1D *hqH1;
  
  vector<double> qH0Collection;
  vector<double> qH1Collection;

};




#endif
