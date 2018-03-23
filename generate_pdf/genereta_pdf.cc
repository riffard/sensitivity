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

//----------------------------------------------------
// NEST
//----------------------------------------------------
//#include "NEST.hh"
//#include "detector.hh"
//#include "detector_LUX.hh"
#include "Nester.hh"


//----------------------------------------------------
//----------------------------------------------------
//#include "PdfCollection.hh"
#include "Xe_wimp.hh"


#include "NeutrinoRate.hh"
#include "NeutrinoCrossSection_coherent_NR.hh"
#include "NeutrinoFlux.hh"

#include "DetectorEfficiency_ER_LUX_Run3.hh"
#include "DetectorEfficiency_NR_LUX_Run3.hh"

//----------------------------------------------------
//----------------------------------------------------
using namespace std;
//using namespace NEST;


//----------------------------------------------------
//----------------------------------------------------
//----------------------------------------------------
//----------------------------------------------------



//----------------------------------------------------
//----------------------------------------------------
int main(int argc,char** argv){
  
  //----------------------------------------------------
  //----------------------------------------------------
  TApplication theApp("App", &argc, argv);

  //----------------------------------------------------
  // General parameters
  //----------------------------------------------------
  string output_path = "pdf";
  int N_nevent_generation = 1000000; //10000000;
  double sigma0Si = 1e-45;    // cm^2
  double low_E_th = 1;        // keV
  
  
  //vector<double> fields = {100, 270, 350}; // V/cm
  vector<double> fields = {100}; // V/cm

  
  //----------------------------------------------------
  // Signal: WIMP parameters
  //----------------------------------------------------
  bool build_wimp = true;
  vector<double> wimp_masses = {10, 50, 100, 500, 1000};
  //vector<double> wimp_masses = {500};
  


  //----------------------------------------------------
  // Calibration: DD gun NR parameters
  //----------------------------------------------------
  bool build_DD_gun_NR = false;
  int N_nevent_generation_DD_gun = 1000000;
  //----------------------------------------------------
  // Background: Neutrino NR parameters
  //----------------------------------------------------
  bool build_neutrino_NR = true;

  //----------------------------------------------------
  // Background: Flat ER parameters
  //----------------------------------------------------
  bool build_flat_ER = true;
  double flat_ER_rate = 1; // unit: Event/ton/year/keV



  Nester* nester = new Nester("LuxRun4TB1");
  
  //----------------------------------------------------
  // Detector efficiency
  //----------------------------------------------------
  DetectorEfficiency* LUX_ER_efficiency = new DetectorEfficiency_ER_LUX_Run3();
  LUX_ER_efficiency->Draw();
  DetectorEfficiency* LUX_NR_efficiency = new DetectorEfficiency_NR_LUX_Run3();
  LUX_NR_efficiency->Draw();

  
  //----------------------------------------------------
  // Target configuration
  //----------------------------------------------------
  FormFactorDataBase* FormEvenNucleus = new FormFactorDataBase("Helm","Helm");
  FormFactorDataBase* Form129Xe = new FormFactorDataBase("Helm","129Xe_p");
  FormFactorDataBase* Form131Xe = new FormFactorDataBase("Helm","131Xe_p");
  
  double fraction, Mass, A, Z, Sp, Sn,J;

  Target* target = new Target("Xenon");  
  target->AddNucleus(fraction = 1e-3,   new Nucleus("Xe124", Mass = 123.9058958,  A = 124, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
  target->AddNucleus(fraction = 9e-4,   new Nucleus("Xe126", Mass = 125.9042689,  A = 126, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
  target->AddNucleus(fraction = 0.0191, new Nucleus("Xe128", Mass = 127.9035304,  A = 128, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
    
  target->AddNucleus(fraction = 0.262,  new Nucleus("Xe129", Mass = 128.9047795,  A = 129, Z = 54, J = 0.5,   Sp = 0.010,  Sn =  0.329, Form129Xe));
  target->AddNucleus(fraction = 0.041,  new Nucleus("Xe130", Mass = 129.9035079,  A = 130, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
    
  target->AddNucleus(fraction = 0.218,  new Nucleus("Xe131", Mass = 130.9050819,  A = 77, Z = 54, J = 3./2., Sp = -0.009, Sn =  -0.272,Form131Xe));

  target->AddNucleus(fraction = 0.269,  new Nucleus("Xe132", Mass = 131.9041545,  A = 132, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
  target->AddNucleus(fraction = 0.104,  new Nucleus("Xe134", Mass = 133.9053945,  A = 134, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
  target->AddNucleus(fraction = 0.089,  new Nucleus("Xe136", Mass = 135.9072195,  A = 136, Z = 54, J = 0,     Sp = 0,      Sn =  0, FormEvenNucleus));
  
  target->CheckNormalization();
  
  
  //----------------------------------------------------
  // Neutrino fluxes
  //----------------------------------------------------
  NeutrinoFlux* neutrino_fluxes = new NeutrinoFlux();
  
  //----------------------------------------------------
  // Neutrino cross section
  //----------------------------------------------------
  NeutrinoCrossSection_coherent_NR* cross_section = new NeutrinoCrossSection_coherent_NR;
    
  
  //----------------------------------------------------
  // Create the output director
  //----------------------------------------------------
  system(Form("mkdir -p %s", output_path.c_str()));
  
  //----------------------------------------------------
  // Build WIMP models
  //----------------------------------------------------
  if(build_wimp){

    cout<<"Generation of WIMP NR signal"<<endl;

    Xe_wimp* xe_wimp = new Xe_wimp(target);

    for(size_t i=0; i<wimp_masses.size(); ++i){
      for(size_t ifield = 0; ifield < fields.size(); ++ifield){

	cout<<"    - m_WIMP = "<< wimp_masses[i] << " GeV, " << fields[ifield] <<" V/cm and csx = "<< sigma0Si <<" cm^2  --> " << flush;

	ostringstream oss;
	oss << wimp_masses[i];
	ostringstream oss_field;
	oss_field << fields[ifield];
	string pdf_name = "wimp_" + oss.str() + "_GeV_" + oss_field.str() + "_Vcm";

	TFile* f = new TFile((output_path + "/" +pdf_name + ".root").c_str(), "RECREATE");
	
	PdfCollection pdf_wimp(pdf_name);
	
	xe_wimp->GetRate(pdf_wimp, wimp_masses[i], sigma0Si, LUX_NR_efficiency);
	
	nester->GetNestedPdf(pdf_wimp, N_nevent_generation, "NR", fields[ifield]);
	
	f->Write();

	cout<< pdf_wimp.hEnergy->Integral("width") << " evt/ton/year  ---- Done"<<endl;
	
      
      }
    }
  }
  //----------------------------------------------------
  //----------------------------------------------------
  if(build_DD_gun_NR){


    
    for(size_t ifield = 0; ifield < fields.size(); ++ifield){
      cout<<" - Generation of DD gun NR calibration @ " << fields[ifield]<<" V/cm  "<<flush;
      
      ostringstream oss_field;
      oss_field << fields[ifield];
      string pdf_name = "dd_gun_" + oss_field.str() + "_Vcm";
      
      TFile* f = new TFile((output_path + "/" +pdf_name + ".root").c_str(), "RECREATE");
      
      PdfCollection pdf_dd(pdf_name);
            
      nester->GetNestedPdf(pdf_dd, N_nevent_generation, "DD", fields[ifield]);
      
      f->Write();
      
      cout<< pdf_dd.hEnergy->Integral("width") << " evt/ton/year  ---- Done"<<endl;
	
      
      }

  }
  

  //----------------------------------------------------
  //----------------------------------------------------
  if(build_neutrino_NR){

    cout<<"Generation of neutrino NR background"<<endl;
    
    for(size_t ifield = 0; ifield < fields.size(); ++ifield){
      
      cout<<"    - @"<<fields[ifield]<<" V/cm " <<flush;

      ostringstream oss_field;
      oss_field << fields[ifield];
      string pdf_name = "neutrino_NR_" + oss_field.str() + "_Vcm";
      TFile* f = new TFile((output_path + "/" +pdf_name + ".root").c_str(), "RECREATE");
      
      
      NeutrinoRate* neutrino_rate = new NeutrinoRate(target, neutrino_fluxes, "All", cross_section, LUX_NR_efficiency);
      //NeutrinoRate* neutrino_rate = new NeutrinoRate(target, neutrino_fluxes, "All", cross_section);
      PdfCollection pdf_neutrino_NR("neutrino_NR");
      neutrino_rate->GetRate(pdf_neutrino_NR.hEnergy);
      nester->GetNestedPdf(pdf_neutrino_NR, N_nevent_generation, "NR", fields[ifield]);

      f->Write();
      
      cout<< pdf_neutrino_NR.hEnergy->Integral("width") << " evt/ton/year"<<" ---- Done"<<endl;

    }
  }
  //----------------------------------------------------
  //----------------------------------------------------

  
  //----------------------------------------------------
  //----------------------------------------------------
  if(build_flat_ER){

    cout<<"Generation of Flat ER background"<<endl;
    
    for(size_t ifield = 0; ifield < fields.size(); ++ifield){

      cout<<"    - @" << fields[ifield] <<"V/cm rate = "<< flat_ER_rate<<" evt/ton/y/keV --> " <<flush;
    
      
      ostringstream oss_field;
      oss_field << fields[ifield];
      string pdf_name = "flat_ER_" + oss_field.str() + "Vcm";
      TFile* f = new TFile((output_path + "/" +pdf_name + ".root").c_str(), "RECREATE");
      
      PdfCollection pdf_flat_ER("flat_ER");
      
      for(int i =1; i<= pdf_flat_ER.hEnergy->GetNbinsX(); ++i) {
	double E = pdf_flat_ER.hEnergy->GetBinCenter(i);
	double rate = LUX_ER_efficiency->GetEfficiency(E) * flat_ER_rate;
	pdf_flat_ER.hEnergy->SetBinContent(i, rate);
	
      }
      
      nester->GetNestedPdf(pdf_flat_ER, N_nevent_generation, "NR", fields[ifield]);

    f->Write();
    
    cout<< pdf_flat_ER.hEnergy->Integral("width") << " evt/ton/year ---- Done"<<endl;

    }
  }
  //----------------------------------------------------
  //----------------------------------------------------



  cout<<"Done !"<<endl;
  //theApp.Run();
  
  return 0;
}
//----------------------------------------------------
//----------------------------------------------------






