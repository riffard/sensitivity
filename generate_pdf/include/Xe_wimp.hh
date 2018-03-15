#ifndef Xe_wimp_hh
#define Xe_wimp_hh


//----------------------------------------------------
// dm_calc
//----------------------------------------------------
#include "Target.hh"
#include "FormFactorDataBase.hh"
#include "HaloModel.hh"
#include "wimp_calculator.hh"
#include "DetectorEfficiency.hh"

//----------------------------------------------------
//----------------------------------------------------
#include "PdfCollection.hh" 

class Xe_wimp{

public:

  Xe_wimp(Target* target){

    fTarget = target;
    
    //-------------------------------------------------------------------------
    // Halo model & parameters
    //-------------------------------------------------------------------------
    double v0 = 220; /*km.s-1*/
    double vesc = 554; /*km.s-1*/
    double sig0 = v0/sqrt(2.);
    
    HaloModel* halo = new HaloModel(v0,sig0 ,vesc, "SMH");


    calc = new wimp_calculator("wimp_calc", target, halo);

    
  };


  void GetRate(PdfCollection&pdfs, double mChi, double sigma0Si, DetectorEfficiency* Efficiency){

    
    for(int i=1; i<=pdfs.hEnergy->GetNbinsX(); ++i){

      double Er = pdfs.hEnergy->GetXaxis()->GetBinCenter(i);
      double rate = calc->Rate_SI( Er, mChi, sigma0Si);

      rate *= Efficiency->GetEfficiency(Er);
      
      pdfs.hEnergy->SetBinContent(i, rate);
      
    }

    
  }
  
  
private:
  Target* fTarget;
  wimp_calculator* calc;
  
  
};


#endif
