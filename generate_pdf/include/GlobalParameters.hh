#ifndef GLOBALPARAMETERS_hh
#define GLOBALPARAMETERS_hh 1

#include <iostream>

using namespace std;

class GlobalParameters{

public:

  static GlobalParameters* GetInstance(){ return (!instance) ? instance = new GlobalParameters() : instance; }
  
private:

  GlobalParameters();
  
  ~GlobalParameters(){}

  
  static  GlobalParameters* instance;

  
public:

  static const int    E_nbins;
  static const double E_min;
  static const double E_max;

  static const bool   is_E_log_scale;
  double* E_log_scale_bins;
  
    
  static const int    s1_nbins;
  static const double s1_max;
  
  static const int    s2_nbins;
  static const double s2_max;
    
  static const int    logs2s1_nbins;
  static const double logs2s1_min;
  static const double logs2s1_max;
 


};

#endif
