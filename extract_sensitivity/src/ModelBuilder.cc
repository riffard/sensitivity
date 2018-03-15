#include "ModelBuilder.hh"

#include "TRandom.h"

ModelBuilder::ModelBuilder(string name, string analysis_variables , map<string, RooRealVar* >* variables_collection){


  fName = name;
  fAnalysisVariables = analysis_variables;
  fVariablesCollection = variables_collection;

  fPdfList = new RooArgList();;
  fFullPdf = NULL;
  fFakeDataSet = NULL;
}



ModelBuilder::~ModelBuilder(){


  
}


void ModelBuilder::AddComponant(string contrib_name, string file_name, string histo_name, RooArgList arglist){

  cout<<"<Info::Builder> ("<< fName<< ") "<<"Add new componant from "<<file_name<< endl;
  
  
  TFile *f = TFile::Open(file_name.c_str());
  if(!f){
    cout<<"Error: "<<file_name<<" cannot be open !"<<endl;
    exit(-1);
  }  
  
  TH1D* hE = (TH1D*)f->Get("E");
  if(!hE){ cout<<"Error: The histo E cannot be found !"<<endl;  exit(-1); }

  double rate = hE->Integral("width");

  fRates.push_back(rate);
  fTotalRate += rate;
  
  cout<<"<Info::Builder> (" << fName <<  ") New componant rate = "<< rate<<" evt/y/t" <<endl;
  
  
  TH2D* h2 = (TH2D*)f->Get(histo_name.c_str());
  if(!h2){ cout<<"Error: The histo "<< histo_name << " cannot be found !"<<endl;  exit(-1); }
   
  fComponants[contrib_name] = h2;
  
  RooDataHist* data = new RooDataHist(Form("data_%s",contrib_name.c_str()),"", arglist, h2);
  //RooDataHist* data = new RooDataHist(Form("data_%s",contrib_name.c_str()),"",RooArgList(fAnalysisVariables.c_str()), h2);
  fComponants_data[contrib_name] = data;

  RooHistPdf* dataPdf = new RooHistPdf(Form("pdf_%s",contrib_name.c_str()),"", arglist, *data);
  //RooHistPdf* dataPdf = new RooHistPdf(Form("pdf_%s",contrib_name.c_str()),"",RooArgList(fAnalysisVariables.c_str()), *data);
  fComponants_pdf[contrib_name] = dataPdf;

  
  RooRealVar* Amplitude = new RooRealVar(Form("A_%s",contrib_name.c_str()),"",rate,0.,10000000);
  (*fVariablesCollection)[Form("A_%s",contrib_name.c_str())] = Amplitude;
  
  
  RooExtendPdf* extendedPdf = new RooExtendPdf(Form("extented_pdf_%s",contrib_name.c_str()),"",*dataPdf, *Amplitude );
  fExtendedPdf[Form("extented_pdf_%s",contrib_name.c_str())] = extendedPdf;
  
  fPdfList->add(*extendedPdf);


  
  
}


void ModelBuilder::BuildPdf(){

  fFullPdf = new RooAddPdf("","",*fPdfList) ;


  cout<<"<Info::Builder> (" << fName <<  ") Build PDF : "<< fTotalRate<<" evt/y/t" <<endl;

}



RooDataSet* ModelBuilder::GenerateFakeData(double exposure, RooArgList arglist){

  if(!fFullPdf) BuildPdf();
  //if(fFakeDataSet) delete fFakeDataSet;

  fExposure = exposure;
  
  double mu_event = fExposure * fTotalRate;
  int N_event = gRandom->Poisson(mu_event);

  
  cout<<"<Info::Builder> (" << fName <<  ") Generate fake data. Exposure="<<exposure<<" t.y : Average = "<< mu_event<<" --> # event " << N_event<<endl;
  
  fFakeDataSet = fFullPdf->generate(arglist ,N_event);
  
  return fFakeDataSet;
}
