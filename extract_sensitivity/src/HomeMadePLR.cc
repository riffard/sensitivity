#include "HomeMadePLR.hh"

#include "TRandom.h"

#include "RooRealVar.h"


#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"

#include "RooStats/RooStatsUtils.h"

HomeMadePLR::HomeMadePLR(RooWorkspace *w){

  fw = w;
 
  
}




HomeMadePLR::~HomeMadePLR(){



}

#include "TCanvas.h"

void HomeMadePLR::Process(int n_pseudo_exp, double multiplicator, HomeMadePLRData& plrData){

  ModelConfig* model = (ModelConfig*) fw->obj("model");
  RooRealVar* poi = (RooRealVar*) model->GetParametersOfInterest()->first();
  RooAbsPdf* pdf = model->GetPdf();

  double rate_init_poi = poi->getVal();
  
  double test_val = multiplicator*rate_init_poi;

  
  cout<<"<> PRL processing for: "<< n_pseudo_exp << " pseudo exp. and " << test_val << endl;
  
  for(int n=0; n<n_pseudo_exp; ++n){

    cout<<"Exp "<< n+1 << " / "<< n_pseudo_exp <<endl;

    
    poi->setVal(test_val);    
    double mu_SigBkg  = pdf->expectedEvents(*fw->set("obs"));
    int N_SigBkg  = gRandom->Poisson(mu_SigBkg);
    RooDataSet* dataSigBkg = pdf->generate(*fw->set("obs"), N_SigBkg );

    
    poi->setVal(0);
    double mu_bkgOnly = model->GetPdf()->expectedEvents(*fw->set("obs"));
    int N_bkgOnly = gRandom->Poisson(mu_bkgOnly);
    RooDataSet* dataBkgOnly = pdf->generate(*fw->set("obs"), N_bkgOnly ); 

    
    double qH0 = Getq(dataSigBkg, model, test_val);
    double qH1 = Getq(dataBkgOnly, model, test_val);
    
    
    plrData.AddData(qH0, qH1);

    poi->setVal(rate_init_poi);
    
  }

  poi->setVal(rate_init_poi);

}


double HomeMadePLR::Getq(RooDataSet* data, ModelConfig* model, double test_val){


  RooAbsPdf* pdf = model->GetPdf();
  RooRealVar* poi = (RooRealVar*) model->GetParametersOfInterest()->first();;
  RooArgSet* np =  (RooArgSet*)model->GetNuisanceParameters();

  int vervose_level = -1;

  int n_np = np->getSize();;
  
  poi->setVal(test_val);
  poi->setConstant(1);

  double lnL_constrain = 0;//fit_constrain->minNll();


  RooArgSet constrainParams;
  if (model->GetNuisanceParameters() ) constrainParams.add(*model->GetNuisanceParameters());
  RooStats::RemoveConstantParameters(&constrainParams);


      
  if(n_np == 0){
    RooAbsReal* fit_constrain = pdf->createNLL(*data);
    lnL_constrain = fit_constrain->getVal();    
  }else{
    RooFitResult* fit_constrain = pdf->fitTo(*data, Extended(kTRUE),
					     Save(1), PrintLevel(vervose_level),
					     //InitialHesse(false), Hesse(false),
					     SumW2Error(kFALSE), Constrain(constrainParams));
        
    lnL_constrain = fit_constrain->minNll();
  }

  
  /*new TCanvas;
  RooPlot* xframe = fw->var("s1")->frame();
  data->plotOn(xframe) ;
  pdf->plotOn(xframe);
  xframe->Draw();
  */
  
  cout<<"<Fit> fit Unconstrain!!"<<endl;
  poi->setConstant(0);
  poi->setVal(test_val);
  RooFitResult* fit_unconstrain = pdf->fitTo(*data, Extended(kTRUE),
					     Save(1), PrintLevel(vervose_level),
					     //InitialHesse(false), Hesse(false),
					     SumW2Error(kFALSE), Constrain(constrainParams)
					     );


  cout<<fit_unconstrain->status()<<endl;
  
  if(fit_unconstrain->status() != 0){
    cout<<"ResMFKhjgh !!!!!"<<endl;
    
  }
  
  double lnL_unconstrain = fit_unconstrain->minNll();
  double poi_bestfit = poi->getVal();
  
  if(poi_bestfit<0){
    
    /*    cout<<"Fit again !!"<<endl;
    fit_unconstrain = pdf->fitTo(*data, Extended(kTRUE),
				 Save(1), PrintLevel(vervose_level),
				 InitialHesse(false), Hesse(false),
				 SumW2Error(kFALSE), Constrain(constrainParams)
				 );
    poi_bestfit = poi->getVal();

    cout<<"hold: "<< poi_bestfit <<endl;
    getchar();*/
  }

  
  
  if(poi_bestfit<0){
    poi->setVal(0);
    RooAbsReal* nLL_unconstrain = pdf->createNLL(*data);
    lnL_unconstrain = nLL_unconstrain->getVal();

    
    
  }

  
  double q= 2 * (lnL_constrain - lnL_unconstrain); 

  return q;
  
}
