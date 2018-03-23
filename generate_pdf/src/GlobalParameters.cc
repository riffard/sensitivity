#include <math.h>
#include <cmath>

#include "GlobalParameters.hh"

GlobalParameters::GlobalParameters(){


    E_log_scale_bins = new double[E_nbins];
    
    double logxmin = log10(E_min);
    double logxmax = log10(E_max);
    double binwidth = (logxmax-logxmin)/(double)E_nbins;
    E_log_scale_bins[0] = E_min;
    for (int i = 1 ; i <= E_nbins ; i++) E_log_scale_bins[i] = E_min + pow(10, logxmin + i * binwidth);


}


GlobalParameters* GlobalParameters::instance = nullptr;

const bool   GlobalParameters::is_E_log_scale = true;
const int    GlobalParameters::E_nbins = 1000;
const double GlobalParameters::E_min = GlobalParameters::is_E_log_scale ?  1e-4: 0; // Only in log scale mode
const double GlobalParameters::E_max = 100;    

const int GlobalParameters::s1_nbins = 1000;
const double GlobalParameters::s1_max = 200;
  
const int GlobalParameters::s2_nbins = 1000;
const double GlobalParameters::s2_max = 16e3;

const int GlobalParameters::logs2s1_nbins = 1000;
const double GlobalParameters::logs2s1_min = 0;
const double GlobalParameters::logs2s1_max = 4;


  
