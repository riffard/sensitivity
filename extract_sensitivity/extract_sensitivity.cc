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
#include "RooWorkspace.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/NumEventsTestStat.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"


#include "Math/PdfFuncMathCore.h"

//----------------------------------------------------
//----------------------------------------------------
#include "GlobalParameters.hh"
#include "ModelBuilder.hh"
#include "Style.hh"
#include "HomeMadePLR.hh"

#include "ConfigReader.hh"
//----------------------------------------------------
//----------------------------------------------------

using namespace std;
using namespace RooFit;
using namespace RooStats;

void PrintWelcome(){

  cout<<"------------------------------------------------------------"<<endl;
  cout<<" _     _   _                     _       _                  "<<endl;
  cout<<"|_)|  |_) | \\ __|_ _  __|_ _ __ |_)_ ___|__ ____  _ __  _ _ "<<endl;
  cout<<"|  |__| \\ |_/(/_|_(/_(_ |_(_)|  | (/_|  |(_)| |||(_|| |(_(/_"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  cout<<"Dev: Q.Riffard (qriffard@lbl.gov)"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
}

void PrintHelp(){

  cout<<"------------------------------------------------------------"<<endl;
  cout<<"Wrong usage of the software !"<<endl;
  cout<<"Usage:"<<endl;
  cout<<"   ./extract_sensitivity.bin -c config_filename -m wimp_mass"<<endl;
  cout<<"Arguments:"<<endl;
  cout<<" - c: set the config file"<<endl;
  cout<<" - m: set the WIMP mass"<<endl;
  cout<<" - f: set the electric field"<<endl;
  cout<<" - b: enable the batch mode (optional)"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  exit(1);
}


//----------------------------------------------------
//----------------------------------------------------
int main(int argc,char** argv){
  
  //----------------------------------------------------
  //----------------------------------------------------
  PrintWelcome();
  LoadStyle();

  string config_file = "";
  string Wimp_mass = "";
  string electric_field = "";
  bool batch_mode = false;
  
  int i=1;
  while(i < argc){

    if((string)argv[i] == "-c"){
      config_file = argv[i+1];
      i++;
    }else if((string)argv[i] == "-m"){
      Wimp_mass =  argv[i+1];
      i++;
    }else if((string)argv[i] == "-f"){
      electric_field =  argv[i+1];
      i++;
    }else if((string)argv[i] == "-b"){
      batch_mode = true;
    }
      
    i++;
  }

  
  
  if(config_file == "" || Wimp_mass == "" || electric_field == "") PrintHelp();
  

  if(batch_mode) gROOT->SetBatch();
  
  //----------------------------------------------------
  //----------------------------------------------------
  
  ConfigReader cfg(config_file);
  //----------------------------------------------------
  //----------------------------------------------------
  TApplication theApp("App", &argc, argv);

  //----------------------------------------------------
  // Global parameters
  //----------------------------------------------------  
  GlobalParameters* parameters = GlobalParameters::GetInstance();
  
  //----------------------------------------------------
  // Model parameters
  //----------------------------------------------------
  const double TotalExpousre = cfg.Get_TotalExpousre(); // ton.year
  const bool UseNuisanceParam = cfg.Get_UseNuisanceParam();
  
  //----------------------------------------------------
  // PLR parameters
  //----------------------------------------------------
  const int N_pseudo_exp = cfg.Get_N_pseudo_exp();
  const int vervose_level = cfg.Get_vervose_level();

  const double p_value_target = cfg.Get_p_value_target(); 
  const double p_value_tolerence = cfg.Get_p_value_tolerence();
  
  //----------------------------------------------------
  // Scaning mode
  //----------------------------------------------------
  //const ScaningMode scaningMode = Dichotomy;
  const ScaningMode scaningMode = cfg.Get_scaningMode();
  
  const double ref_cross_section = parameters->ref_wimp_SI_cx;
  
  const double grid_cross_section_start = ref_cross_section;
  const double grid_cross_section_end = cfg.Get_grid_cross_section_end();
  const int grid_Nsteps = cfg.Get_grid_Nsteps();

  const double dico_cross_section_start = ref_cross_section;
  const int dico_Nsteps_max = cfg.Get_dico_Nsteps_max();
  double dico_cross_section_step = 0.5*ref_cross_section;


  ostringstream oss;
  oss << N_pseudo_exp ;
  string output_filename = "PLR_results_" + Wimp_mass + "GeV_" + electric_field + "Vcm_"+ oss.str() +"Exp.root";

  TFile* f = new TFile(output_filename.c_str(), "RECREATE");

  cout<<"------------------------------------------------------------"<<endl;
  cout<<"Calculation parametes:"<<endl;
  cout<<"------------------------------------------------------------"<<endl;
  
  
  //----------------------------------------------------
  // Preparation of the scaning grid
  //----------------------------------------------------  
  vector<double> scaning_grid;
 
  if(scaningMode == Grid){

    double log_start = log10(grid_cross_section_start);
    double log_end = log10(grid_cross_section_end);
    double log_step = (log_end - log_start) / double(grid_Nsteps);

    for(double log_cx = log_start; log_cx >= log_end; log_cx += log_step)
      scaning_grid.push_back( pow(10, log_cx) );
    
  }else if(scaningMode == Dichotomy){
    scaning_grid.push_back(dico_cross_section_start);
  }
  
  const int MaxStep = scaningMode == Grid ? scaning_grid.size() : dico_Nsteps_max;


 
  
  //----------------------------------------------------
  // Output
  //----------------------------------------------------
  vector<HomeMadePLRData> plrData;  
    
  TGraph* gr_pvalue = new TGraph();
  TGraph* gr_pvalue_inv = new TGraph();
  
  
  //----------------------------------------------------
  //----------------------------------------------------
  RooWorkspace *w = new RooWorkspace("w");

  //----------------------------------------------------
  //----------------------------------------------------
  RooRealVar exposure("exposure","exposure [t.y]",TotalExpousre);
  w->import(exposure);
  
  RooRealVar s1("s1", "s1 [npe]", 0, parameters->s1_max) ;
  w->import(s1);
    
  RooRealVar logs2s1("logs2s1", "logs2s1", parameters->logs2s1_min, parameters->logs2s1_max) ;
  w->import(logs2s1);

  //----------------------------------------------------
  //---------------------------------------------------
  ModelBuilder* model = new ModelBuilder(w, "Model", TotalExpousre, "s1_logs1s2", "s1,logs2s1");
  
  model->AddComponant(ModelBuilder::Background, "NR_neutrino", Form("../generate_pdf/pdf/neutrino_NR_%s_Vcm.root", electric_field.c_str()));  
  model->AddComponant(ModelBuilder::Background, "Rn222", Form("../generate_pdf/pdf/Rn222_ER_%sVcm.root", electric_field.c_str()));
  model->AddComponant(ModelBuilder::Background, "Rn220", Form("../generate_pdf/pdf/Rn220_ER_%sVcm.root", electric_field.c_str()));
  
  model->AddComponant(ModelBuilder::Signal, "wimp", Form("../generate_pdf/pdf/wimp_%s_GeV_%s_Vcm.root",Wimp_mass.c_str(), electric_field.c_str()));

  //----------------------------------------------------
  //----------------------------------------------------
  model->BuildModel(UseNuisanceParam);


  //----------------------------------------------------
  //----------------------------------------------------
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::DataHandling);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments);

  
  //----------------------------------------------------
  //----------------------------------------------------
  HomeMadePLR* plr = new HomeMadePLR(w);


  //----------------------------------------------------
  //----------------------------------------------------

  int n_scaning = scaning_grid.size();
  
  for(size_t i=0; i < n_scaning && i < MaxStep; ++i){

    cout<<"i = "<<i<<endl;
    
    f->cd();
    double test_cross_section = scaning_grid[i];
    
    string name = Form("%e", test_cross_section);
    plrData.push_back(HomeMadePLRData (name));   

    double multiplicator = test_cross_section / ref_cross_section; 
    
    plr->Process(N_pseudo_exp, multiplicator, plrData[i]);
    
    double p_value = plrData[i].GetPValue();
    gr_pvalue->SetPoint(gr_pvalue->GetN(), test_cross_section, p_value);
    gr_pvalue_inv->SetPoint(gr_pvalue_inv->GetN(), p_value, test_cross_section);

  
    if(scaningMode == Grid) continue;
    
    double dist_to_target = abs(p_value_target - p_value);
    
    if(dist_to_target < p_value_tolerence ){
      cout<<"Optimium found!! "<< p_value << " target: "<< p_value_target <<" ("<<dist_to_target<<")"<<endl;
      break;
    }
    
    if(p_value > p_value_target && dico_cross_section_step > 0){
      dico_cross_section_step = -dico_cross_section_step/2.;
    }else if(p_value < p_value_target && dico_cross_section_step < 0){
      dico_cross_section_step = -dico_cross_section_step/2.;
    }

    scaning_grid.push_back( scaning_grid.back() - dico_cross_section_step);

    
    n_scaning = scaning_grid.size();
    
  }

  gr_pvalue->Write("pvalue");
  gr_pvalue_inv->Write("pvalue_inv");
  
  f->Write();
  f->Close();
  
  cout<<"Done "<<endl;

  if(!batch_mode) theApp.Run();
  

  






} 
