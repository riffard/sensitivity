#include "ModelBuilder.hh"


ModelBuilder::ModelBuilder(map<string, RooRealVar* >* variables){


  fVariables = variables;
  
}



ModelBuilder::~ModelBuilder(){


  
}


void ModelBuilder::AddComponant(string file_name, string histo_name, RooArgList arglist){


  TFile *f = TFile::Open(file_name.c_str());

  TH2D* h2 = (TH2D*)f->Get(histo_name.c_str());


  fComponants[histo_name] = h2;
    
  RooDataHist* data = new RooDataHist(Form("data_%s",histo_name.c_str()),"",arglist, h2);
  fComponants_data[histo_name] = data;
  
  RooHistPdf* dataPdf = new RooHistPdf(Form("pdf_%s",histo_name.c_str()),"",arglist, *data);
  fComponants_pdf[histo_name] = dataPdf;
  
  RooRealVar* Amplitude = new RooRealVar(Form("A_%s",histo_name.c_str()),"",200,0.,10000) ;

  
}


void BuildPdf(){
  
  
}
