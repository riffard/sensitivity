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
#include "Nester.hh"

//----------------------------------------------------
//----------------------------------------------------
#include "GlobalParameters.hh"
#include "ReadFromFile.hh"
#include "Discrimination.hh"

#include "WimpRate.hh"
#include "NeutrinoRate.hh"
#include "NeutrinoCrossSection_coherent_NR.hh"
#include "NeutrinoFlux.hh"

#include "DetectorEfficiency_ER_LUX_Run3.hh"
#include "DetectorEfficiency_NR_LUX_Run3.hh"


//----------------------------------------------------
//----------------------------------------------------
using namespace std;

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
  int N_nevent_generation = 10000000;
  double sigma0Si = GlobalParameters::GetInstance()->ref_wimp_SI_cx;    // cm^2
  double low_E_th = 1;        // keV
  double ER_rejection = 0.995;
  double NR_rejection = (1-0.5);
  
  
  
  //vector<double> fields = {100, 270, 350}; // V/cm
  vector<double> fields = {100}; // V/cm
  
  //----------------------------------------------------
  // Signal: WIMP parameters
  //----------------------------------------------------
  bool build_wimp = true;
  vector<double> wimp_masses = {10, 50, 100, 500, 1000};
  //vector<double> wimp_masses = {500};


  //-------------------------------------------------------------------------
  // Halo model & parameters
  //-------------------------------------------------------------------------
  double v0 = 220; /*km.s-1*/
  double vesc = 554; /*km.s-1*/
  double sig0 = v0/sqrt(2.);
  
  HaloModel* halo = new HaloModel(v0,sig0 ,vesc, "SMH");
  

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
  double flat_ER_rate = 0.03; // unit: Event/ton/year/keV

  //----------------------------------------------------
  // Background: Rn22* parameters
  //----------------------------------------------------
  bool build_Rn222 = true;
  bool build_Rn220 = true;
  
  string Rn222_file_name = "data_base/Rn222.txt";
  string Rn220_file_name = "data_base/Rn220.txt";

  //----------------------------------------------------
  //----------------------------------------------------
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
  // WIMP rate
  //----------------------------------------------------
  WimpRate* wimpRate = new WimpRate("wimp_calc", target, halo, LUX_NR_efficiency);
  
  
  //----------------------------------------------------
  // Neutrino Claculator parameters
  //----------------------------------------------------
  // Neutrino fluxes
  NeutrinoFlux* neutrino_fluxes = new NeutrinoFlux();
  
  // Neutrino cross section
  NeutrinoCrossSection_coherent_NR* cross_section = new NeutrinoCrossSection_coherent_NR;    

  // Neutrino rate calculator
  NeutrinoRate* neutrino_rate = new NeutrinoRate(target, neutrino_fluxes, "All", cross_section, LUX_NR_efficiency);

  
  //----------------------------------------------------
  // Create the output director
  //----------------------------------------------------
  system(Form("mkdir -p %s", output_path.c_str()));
  
  //----------------------------------------------------
  // Build WIMP models
  //----------------------------------------------------
  if(build_wimp){

    cout<<"Generation of WIMP NR signal"<<endl;

    //Xe_wimp* xe_wimp = new Xe_wimp(target);

    for(size_t i=0; i<wimp_masses.size(); ++i){
      for(size_t ifield = 0; ifield < fields.size(); ++ifield){

	cout<<"    - m_WIMP = "<< wimp_masses[i] << " GeV, " << fields[ifield] <<" V/cm and csx = "<< sigma0Si <<" cm^2  --> " << flush;

	ostringstream oss;
	oss << wimp_masses[i];
	ostringstream oss_field;
	oss_field << fields[ifield];
	string pdf_name = "wimp_" + oss.str() + "_GeV_" + oss_field.str() + "_Vcm";

	TFile* f = new TFile((output_path + "/" +pdf_name + ".root").c_str(), "RECREATE");
	
	PdfCollection pdf(pdf_name);
	
	wimpRate->GetRate(pdf.hEnergy, wimp_masses[i], sigma0Si);

	AddFlatDiscrimination(NR_rejection, pdf);
	
	nester->GetNestedPdf(pdf, N_nevent_generation, "NR", fields[ifield]);
	
	f->Write();

	cout<< pdf.hEnergy->Integral("width") << " evt/ton/year  ---- Done"<<endl;
	
      
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
      
      PdfCollection pdf(pdf_name);
            
      nester->GetNestedPdf(pdf, N_nevent_generation, "DD", fields[ifield]);
      
      f->Write();
      
      cout<< pdf.hEnergy->Integral("width") << " evt/ton/year  ---- Done"<<endl;
	
      
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
      
      
      PdfCollection pdf("neutrino_NR");
      neutrino_rate->GetRate(pdf.hEnergy);

      AddFlatDiscrimination(NR_rejection, pdf);
      
      nester->GetNestedPdf(pdf, N_nevent_generation, "NR", fields[ifield]);

      f->Write();
      
      cout<< pdf.hEnergy->Integral("width") << " evt/ton/year"<<" ---- Done"<<endl;

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
      
      PdfCollection pdf("flat_ER");
      
      for(int i =1; i<= pdf.hEnergy->GetNbinsX(); ++i) {
	double E = pdf.hEnergy->GetBinCenter(i);
	double rate = LUX_ER_efficiency->GetEfficiency(E) * flat_ER_rate;
	pdf.hEnergy->SetBinContent(i, rate);
	
      }
      AddFlatDiscrimination(ER_rejection, pdf);

      nester->GetNestedPdf(pdf, N_nevent_generation, "ER", fields[ifield]);

    f->Write();
    
    cout<< pdf.hEnergy->Integral("width") << " evt/ton/year ---- Done"<<endl;

    }
  }
  //----------------------------------------------------
  //----------------------------------------------------

  //----------------------------------------------------
  //----------------------------------------------------
  if(build_Rn220){
    cout<<"Generation of Rn220"<<endl;
    
    for(size_t ifield = 0; ifield < fields.size(); ++ifield){
      
      ostringstream oss_field;
      oss_field << fields[ifield];
      string pdf_name = "Rn220_ER_" + oss_field.str() + "Vcm";
      TFile* f = new TFile((output_path + "/" +pdf_name + ".root").c_str(), "RECREATE");
      
      PdfCollection pdf("Rn220_ER");

      ReadFromFile(Rn220_file_name, pdf);
      AddFlatDiscrimination(ER_rejection, pdf);

      nester->GetNestedPdf(pdf, N_nevent_generation, "ER", fields[ifield]);

      f->Write();
      
      cout<< pdf.hEnergy->Integral("width") << " evt/ton/year ---- Done"<<endl;

    }
  }
  //----------------------------------------------------
  //----------------------------------------------------

  //----------------------------------------------------
  //----------------------------------------------------
  if(build_Rn222){
    cout<<"Generation of Rn222"<<endl;
    
    for(size_t ifield = 0; ifield < fields.size(); ++ifield){
      
      ostringstream oss_field;
      oss_field << fields[ifield];
      string pdf_name = "Rn222_ER_" + oss_field.str() + "Vcm";
      TFile* f = new TFile((output_path + "/" +pdf_name + ".root").c_str(), "RECREATE");
      
      PdfCollection pdf("Rn222_ER");

      ReadFromFile(Rn222_file_name, pdf);
      AddFlatDiscrimination(ER_rejection, pdf);
      
      nester->GetNestedPdf(pdf, N_nevent_generation, "ER", fields[ifield]);

      f->Write();
      
      cout<< pdf.hEnergy->Integral("width") << " evt/ton/year ---- Done"<<endl;

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






