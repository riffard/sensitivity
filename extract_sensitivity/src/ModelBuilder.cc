#include "ModelBuilder.hh"

#include "RooStats/ModelConfig.h"

#include "TRandom.h"


using namespace RooStats;
using namespace RooFit;

ModelBuilder::ModelBuilder(RooWorkspace *w, string name, double exposure, string histo_name, string variables, int verbose){

  fw = w;
  fName = name;
  fExposure = exposure;
  fVerbose = verbose;
  fHistoName = (char*)histo_name.c_str(); 
  fVariables = (char*)variables.c_str();
   
}



ModelBuilder::~ModelBuilder(){


  
}


//----------------------------------------------------
//----------------------------------------------------
void ModelBuilder::AddComponant(ComponantType componantType, string _contrib_name, string _file_name){

  const char* contrib_name = _contrib_name.c_str();
  const char* file_name = _file_name.c_str();
  
  double error_factor = 0.1;

  
  cout<<"<Info::Builder> ("<< fName<< ") "<<"Add new componant from "<<file_name<< endl;

  //----------------------------------------------------
  //----------------------------------------------------
  TFile *f = TFile::Open(file_name);
  if(!f){
    cout<<"Error: "<<file_name<<" cannot be open !"<<endl;
    exit(-1);
  }  

  //----------------------------------------------------
  // Get the TH1D for the normalization
  // ** This step must be discarded at some point and the TH2D/obs must be normalized correctly ** 
  //----------------------------------------------------
  TH1D* hE = (TH1D*)f->Get("E");
  if(!hE){ cout<<"Error: The histo E cannot be found !"<<endl;  exit(-1); }

  //----------------------------------------------------
  // Get the TH2D
  // ** This step must be discarded at some point and the TH2D/obs must be normalized correctly ** 
  //----------------------------------------------------
  TH2D* h2 = (TH2D*)f->Get(fHistoName);
  if(!h2){ cout<<"Error: The histo "<< fHistoName << " cannot be found !"<<endl;  exit(-1); }

  //----------------------------------------------------
  // Compute rate and expected event count
  //----------------------------------------------------
  double rate = hE->Integral("width");
  fRates[contrib_name] = rate;
  
  double Ncount = rate * fExposure;

  cout<<"<Info::Builder> (" << fName <<  ") New componant rate = "<< rate<<" evt/y/t  -- > "<< Ncount <<" events" <<endl;


  //----------------------------------------------------
  //----------------------------------------------------
  RooRealVar N(Form("N_%s", contrib_name), Form("Number of %s bkg events", contrib_name), Ncount, 0, 5000);
  
  if(componantType == Background){
    
    fw->import(N);
      
    RooRealVar N0(Form("N0_%s", contrib_name),Form("Estimated number of %s bkg events from aux measurement", contrib_name), Ncount);
    fw->import(N0);

    RooRealVar sigmaN0(Form("sigmaN0_%s", contrib_name), Form("Absolute error on expected %s BG", contrib_name), error_factor*Ncount);
    fw->import(sigmaN0); 
    
    fw->factory(Form("Gaussian::%sConstraint(N_%s,N0_%s,sigmaN0_%s)", contrib_name, contrib_name, contrib_name, contrib_name));
    
    fNuisanceParameters.push_back(Form("N_%s", contrib_name));
    fNuisanceParametersConstrain.push_back(Form("%sConstraint", contrib_name));
    fNuisanceParametersGlobalObs.push_back(Form("N0_%s", contrib_name));
  }
  else{
    fInterestParameters.push_back(Form("N_%s", contrib_name));

    N.setRange(-1000,10000);
    fw->import(N);

  }

  
  RooDataHist sig(Form("sig%s", contrib_name),Form("sig%s", contrib_name), fw->argSet(fVariables), h2);
  RooHistPdf HistPdf(Form("%sHistPdf", contrib_name) ,Form("%sHistPdf",contrib_name) , fw->argSet(fVariables),sig);
  fw->import(HistPdf);

  fModelComponant.push_back(Form("N_%s*%sHistPdf",contrib_name,contrib_name));
  
}
//----------------------------------------------------
//----------------------------------------------------


//----------------------------------------------------
//----------------------------------------------------
void ModelBuilder::BuildModel(bool IsNPActivated){


  if(IsNPActivated)
    cout<<"<Info::Builder> (" << fName <<  ") Build the model with nuisance parameters" <<endl;
  else
    cout<<"<Info::Builder> (" << fName <<  ") Build the model without nuisance parameters" <<endl;
    
  
  //----------------------------------------------------
  // Create the model
  //----------------------------------------------------
  string fullModel = "SUM::EventModel(";
  fullModel += getStringFromArray(fModelComponant);
  fullModel += ")";

  //----------------------------------------------------
  // Register the model
  //----------------------------------------------------
  fw->factory(fullModel.c_str());


  //----------------------------------------------------
  // Create the nuisance parameters
  //----------------------------------------------------
  string nuisanceParam = "PROD::FullModel(EventModel";
  if(IsNPActivated) nuisanceParam += "," + getStringFromArray(fNuisanceParametersConstrain); 
  nuisanceParam += ")";
  
  //----------------------------------------------------
  // Register the nuisance parameters
  //----------------------------------------------------
  fw->factory(nuisanceParam.c_str());

 
  //----------------------------------------------------
  // Make the configuration for RooStats::ModelConfig
  //----------------------------------------------------

  //Set the observables
  fw->defineSet("obs",fVariables);

  //Set the parameters of interest
  string poi_str = getStringFromArray(fInterestParameters);  
  fw->defineSet("poi",poi_str.c_str());

  //Setup the nuisance parameteres
  if(IsNPActivated){

    //Global Observables: the auxiliary estimates of the BG rates
    string global_obs = getStringFromArray(fNuisanceParametersGlobalObs);
    fw->defineSet("gobs",global_obs.c_str());

    //Nuisance Paramenters: free parameters the value of which value will be estimated from the data
    string nuisance = getStringFromArray(fNuisanceParameters);
    fw->defineSet("nuis",nuisance.c_str());

    
  }else{
    
    fw->defineSet("gobs","");
    fw->defineSet("nuis","");

    for(size_t i=0; i< fNuisanceParameters.size(); ++i) fw->var(fNuisanceParameters[i].c_str())->setConstant(1);
  }


  ModelConfig* model = new ModelConfig("model", fw);
  model->SetPdf(*fw->pdf("FullModel"));
  model->SetObservables(*fw->set("obs"));
  model->SetParametersOfInterest(*fw->set("poi"));
  model->SetNuisanceParameters(*fw->set("nuis"));
  model->SetGlobalObservables(*fw->set("gobs"));
  fw->import(*model);
  //Note: you don't need to create a snapshop because the value of muSig will be varied during the run                                                                                                                                                                        

  ModelConfig* bModel = (ModelConfig*) model->Clone("bModel");
  bModel->SetName("bModel");
  RooRealVar* poi = (RooRealVar*) model->GetParametersOfInterest()->first();

  cout<<"poi->GetName = "<<poi->GetName()<<endl;;
  
  float oldval = poi->getVal();
  poi->setVal(0.);
  bModel->SetSnapshot(*poi);
  poi->setVal(oldval);
  fw->import(*bModel);


  cout<<"<Info::Builder> (" << fName <<  ") Model construction completed. " <<endl;
  cout<<"-----------------------------------------------------------------------" <<endl;
  cout<<"-----------------------------------------------------------------------" <<endl;
  fw->Print();
  cout<<"-----------------------------------------------------------------------" <<endl;
  cout<<"-----------------------------------------------------------------------" <<endl;
  


  
}
//----------------------------------------------------
//----------------------------------------------------



//----------------------------------------------------
//----------------------------------------------------
void ModelBuilder::BuildData(){

  RooAbsData* obsData;

  
  ModelConfig* mcB = (ModelConfig*) fw->obj("bModel");
  ModelConfig* mcSB = (ModelConfig*) fw->obj("model");
  
  double n0_obs = mcB->GetPdf()->expectedEvents(*fw->set("obs"));

  bool DISCOVERY_SENS_FLAG = true;
  
  if (bool(DISCOVERY_SENS_FLAG)){
    cout << "+++Generating random data with extended term from ALT model (S+B)+++" << endl;
    n0_obs = mcSB->GetPdf()->expectedEvents(*fw->set("obs"));
    obsData = (RooDataSet*) mcSB->GetPdf()->generate(*fw->set("obs"), n0_obs ,Extended(true),Verbose(1));
  }
  else{
    cout << "+++Generating random data with extended term from ALT model (B)+++" << endl;
    n0_obs = mcB->GetPdf()->expectedEvents(*fw->set("obs"));
    obsData = (RooDataSet*) mcB->GetPdf()->generate(*fw->set("obs"), n0_obs ,Extended(true),Verbose(1));
    
  }
  cout << "Number expected counts in observed data set = " << n0_obs << endl;

  char* extraOutputName = "toto";
  double mWimp = 100;
  int  livetime = 10;
  char* filename = Form("genData_Freq_%s_m%.0f_d%i",extraOutputName,mWimp,livetime);

  
  obsData->Print();
  obsData->SetName("obsData");
  string filename_with_ext =  (string)filename + (string)".root";
  TFile* fdata = new TFile(filename_with_ext.c_str(),"recreate");
  obsData->Write();
  cout << "Data file created: " << filename_with_ext << endl;
  fdata->Close();

  fw->import(*obsData);
  cout << "---------------------------------" << endl;
  
  fw->Print();
  //TString workspace_filename = Form("%s/Workspace_Dir/Workspace-m%.0f-d%i.root",outFolder,mWimp,livetime);
  //fw->SaveAs(workspace_filename);
  cout << "Workspace successfully created" << endl;
  
}
//----------------------------------------------------
//----------------------------------------------------

