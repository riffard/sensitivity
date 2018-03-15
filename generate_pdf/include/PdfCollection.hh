#ifndef PdfCollection_hh
#define PdfCollection_hh

#include <iostream>

#include "TH1D.h"
#include "TH2D.h"

#include "GlobalParameters.hh"


using namespace std;

struct PdfCollection{

  PdfCollection(string pdf_name){
    
    GlobalParameters* param = GlobalParameters::GetInstance();

    
    /*    if(param->is_E_log_scale)
      hEnergy = new TH1D(Form("E_%s", pdf_name.c_str()), "", param->E_nbins, param->E_log_scale_bins);
    else
      hEnergy = new TH1D(Form("E_%s", pdf_name.c_str()), "", param->E_nbins, 0, param->E_max);
    
    hs1 = new TH1D(Form("s1_%s", pdf_name.c_str()), "", param->s1_nbins, 0, param->s1_max);
    hs2 = new TH1D(Form("s2_%s", pdf_name.c_str()), "", param->s2_nbins, 0, param->s2_max);
    htdrift = new TH1D(Form("tdrift_%s", pdf_name.c_str()), "", 1000, 0, 500);

    h2_s1_s2 = new TH2D( Form("s1_s2_%s", pdf_name.c_str()), "", param->s1_nbins, 0, param->s1_max, param->s2_nbins, 0, param->s2_max);
    h2_s1_logs2s1 = new TH2D( Form("s1_logs1s2_%s", pdf_name.c_str()), "", param->s1_nbins, 0, param->s1_max, param->logs2s1_nbins, param->logs2s1_min, param->logs2s1_max);
    
    */
    if(param->is_E_log_scale)
      hEnergy = new TH1D("E", "", param->E_nbins, param->E_log_scale_bins);
    else
      hEnergy = new TH1D("E", "", param->E_nbins, 0, param->E_max);
    
    hs1 = new TH1D("s1", "", param->s1_nbins, 0, param->s1_max);
    hs2 = new TH1D("s2", "", param->s2_nbins, 0, param->s2_max);
    htdrift = new TH1D("tdrift" , "", 1000, 0, 500);

    h2_s1_s2 = new TH2D( "s1_s2" , "", param->s1_nbins, 0, param->s1_max, param->s2_nbins, 0, param->s2_max);
    h2_s1_logs2s1 = new TH2D( "s1_logs1s2", "", param->s1_nbins, 0, param->s1_max, param->logs2s1_nbins, param->logs2s1_min, param->logs2s1_max);
    

  };


  void Apply_low_Energy_th(double){}

  
  
  TH1D* hEnergy, *hs1, *hs2, *htdrift;
  TH2D* h2_s1_s2, *h2_s1_logs2s1;

  

  
  
  
};


#endif
