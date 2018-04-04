#ifndef HomeMadePLR_hh
#define HomeMadePLR_hh l

#include <iostream>

#include "TH1D.h"

#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooStats/ModelConfig.h"

#include "HomeMadePLRData.hh"


using namespace RooStats;
using namespace RooFit;



using namespace std;

class HomeMadePLR{

public:
   
  HomeMadePLR(RooWorkspace *w);

  ~HomeMadePLR();
  
  void Process(int n_pseudo_exp, double text_val, HomeMadePLRData& plrData);

private:
  RooWorkspace *fw;

  double Getq(RooDataSet* dataSigBkg, ModelConfig* model, double test_val);
  
};


#endif
