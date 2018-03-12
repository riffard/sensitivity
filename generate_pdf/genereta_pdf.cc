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
#include "NEST.h"


//----------------------------------------------------
//----------------------------------------------------
#include "PdfCollection.hh"
#include "Xe_wimp.hh"


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
void get_nested_pdf(PdfCollection& pdfs, NEST &myNEST, int Nevent, string pt){
  

  int particleType;
  if (pt =="DD" || pt =="NR"){
    particleType = 0;
  }
  else if (pt =="CH3T" || pt == "ER"){
    particleType =1;
  }
  else{
    std::cout << "WARNING YOU HAVE SET AND ILLEGAL PARTICLE TYPE!" << std::endl;
    std::cout << "Please use DD or NR" << std::endl;
    exit(-1);
    
  }

  myNEST.SetParticleType(particleType);

  for(int i=0; i<Nevent; i++){
    
    double energy = pdfs.hEnergy->GetRandom();
    
    myNEST.SetEnergy(energy);
    
    myNEST.DetectorResponse();
    
    double drift = myNEST.GetDriftLocation();
    double S1c = myNEST.GetS1c();
    double S1raw = myNEST.GetS1();
    double S2c = myNEST.GetS2c();
    double S2raw = myNEST.GetS2();

    pdfs.hs1->Fill(S1c);
    pdfs.hs2->Fill(S2c);
    pdfs.htdrift->Fill(drift);
    pdfs.h2_s1_s2->Fill(S1c, S2c);
    pdfs.h2_s1_logs2s1->Fill(S1c, log10(S2c/S1c));

       
  }

  //Spectra normalization
  
  double Total_rate = pdfs.hEnergy->Integral("width");
  
  pdfs.hs1->Scale(Total_rate / pdfs.hs1->Integral("width"));
  pdfs.hs2->Scale(Total_rate / pdfs.hs2->Integral("width"));
  pdfs.htdrift->Scale(Total_rate / pdfs.htdrift->Integral("width"));
  pdfs.h2_s1_s2->Scale(Total_rate / pdfs.h2_s1_s2->Integral("width"));
  pdfs.h2_s1_logs2s1->Scale(Total_rate / pdfs.h2_s1_logs2s1->Integral("width"));
  
  
}
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
  int N_nevent_generation = 10000000;
  double sigma0Si = 1e-45;    // cm^2
  double low_E_th = 1;        // keV
  
  //----------------------------------------------------
  // Signal: WIMP parameters
  //----------------------------------------------------
  bool build_wimp = true;
  vector<double> wimp_masses = {50, 100, 1000};
  
  //----------------------------------------------------
  // Background: Neutrino NR parameters
  //----------------------------------------------------
  bool build_neutrino_NR = true;

  //----------------------------------------------------
  // Background: Flat ER parameters
  //----------------------------------------------------
  bool build_flat_ER = true;
  double flat_ER_rate = 1; // unit: Event/ton/year/keV

  

  //----------------------------------------------------
  // NEST configuration
  //----------------------------------------------------
  Detector myDetector;
  myDetector.LZSettings();
  myDetector.cathodeVoltage=100; //Example to modify the cathodeVoltage for Skin Response

  double df = 200; // drift field
  double dt = -1; // drift time init
  NEST myNEST(0, 1, df, 2.88, dt); //particle type (0==NR 1 ==ER), energy deposited (keV), e-field (V/cm), xe density (g/cm^3)
  
  myNEST.SetDetectorParameters(myDetector);
  myNEST.SetDTmin(40.); // minimum drift in us
  myNEST.SetDTmax(300.); // maximum drift in us
  
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
      
      ostringstream oss;
      oss << wimp_masses[i];
      string pdf_name = "wimp_" + oss.str();
      
      TFile* f = new TFile((output_path + "/" +pdf_name + ".root").c_str(), "RECREATE");
      
      PdfCollection pdf_wimp(pdf_name);
      
      xe_wimp->GetRate(pdf_wimp, wimp_masses[i], sigma0Si, LUX_NR_efficiency);
      
      get_nested_pdf(pdf_wimp, myNEST, N_nevent_generation, "NR");

	
      f->Write();


      cout<<"    - m_WIMP = "<< wimp_masses[i] << " GeV and csx = "<< sigma0Si <<" cm^2  --> " << pdf_wimp.hEnergy->Integral("width") << " evt/ton/year ---- Done"<<endl;
      
  

    }
  }
  //----------------------------------------------------
  //----------------------------------------------------


  //----------------------------------------------------
  //----------------------------------------------------
  if(build_neutrino_NR){

    cout<<"Generation of neutrino NR background"<<endl;

    string pdf_name = "neutrino_NR";
    TFile* f = new TFile((output_path + "/" +pdf_name + ".root").c_str(), "RECREATE");
        
    
    NeutrinoRate* neutrino_rate = new NeutrinoRate(target, neutrino_fluxes, "All", cross_section, LUX_NR_efficiency);
    
    PdfCollection pdf_neutrino_NR("neutrino_NR");

    neutrino_rate->GetRate(pdf_neutrino_NR.hEnergy);
    
    get_nested_pdf(pdf_neutrino_NR, myNEST, N_nevent_generation, "NR");
    
    f->Write();

    cout<<"    - " << pdf_neutrino_NR.hEnergy->Integral("width") << " evt/ton/year ---- Done"<<endl;
    
  }
  //----------------------------------------------------
  //----------------------------------------------------

  
  //----------------------------------------------------
  //----------------------------------------------------
  if(build_flat_ER){

    cout<<"Generation of Flat ER background"<<endl;
    
    string pdf_name = "flat_ER";
    TFile* f = new TFile((output_path + "/" +pdf_name + ".root").c_str(), "RECREATE");
       
    PdfCollection pdf_flat_ER("flat_ER");

    for(int i =1; i<= pdf_flat_ER.hEnergy->GetNbinsX(); ++i) {
      double E = pdf_flat_ER.hEnergy->GetBinCenter(i);
      double rate = LUX_ER_efficiency->GetEfficiency(E) * flat_ER_rate;
      pdf_flat_ER.hEnergy->SetBinContent(i, rate);

    }
      
    get_nested_pdf(pdf_flat_ER, myNEST, N_nevent_generation, "ER");
    
    f->Write();
   
    cout<<"    - rate = "<< flat_ER_rate<<" evt/ton/y/keV --> " << pdf_flat_ER.hEnergy->Integral("width") << " evt/ton/year ---- Done"<<endl;
    
  }
  //----------------------------------------------------
  //----------------------------------------------------



  cout<<"Done !"<<endl;
  //theApp.Run();
  
  return 0;
}
//----------------------------------------------------
//----------------------------------------------------
