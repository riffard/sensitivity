#ifndef ReadFromFile_hh
#define ReadFromFile_hh 1

#include "Math/Interpolator.h"

#include "PdfCollection.hh"

void ReadFromFile(string filename, PdfCollection& pdf){

  cout<<"Read from the file: "<< filename  << endl;
  
  ifstream file(filename.c_str());
  
  vector<double> vec_E, vec_rate;
  
  string tmp;
  getline(file, tmp);
  
  double energy_keV, rate_dru;
  while(file >> energy_keV >> rate_dru){

    double rate_kyt = rate_dru * 1000. * 365.;
    
    vec_E.push_back(energy_keV);
    vec_rate.push_back(rate_kyt);
  }

  ROOT::Math::Interpolator* interp = new ROOT::Math::Interpolator(vec_E, vec_rate);

  for(int i = 0; i <= pdf.hEnergy->GetNbinsX(); ++i){

    double E = pdf.hEnergy->GetBinCenter(i);
    if(E < vec_E[0]) continue;
    if(E > vec_E.back()) break;
    
    pdf.hEnergy->SetBinContent(i, interp->Eval(E));
    
  }

}

#endif
