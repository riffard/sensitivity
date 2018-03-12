#ifndef GlobalParameters_hh
#define GlobalParameters_hh 1


class GlobalParameters{

public:

  static GlobalParameters* GetInstance(){ return (!instance) ? instance = new GlobalParameters() : instance; }
  
private:

  GlobalParameters(){

    E_log_scale_bins = new double[E_nbins];
    
    double logxmin = log10(E_min);
    double logxmax = log10(E_max);
    double binwidth = (logxmax-logxmin)/(double)E_nbins;
    E_log_scale_bins[0] = E_min;
    for (int i = 1 ; i <= E_nbins ; i++) E_log_scale_bins[i] = E_min + pow(10, logxmin + i * binwidth);
    
    
  }

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

GlobalParameters* GlobalParameters::instance = nullptr;

const bool   GlobalParameters::is_E_log_scale = false;
const int    GlobalParameters::E_nbins = 1000;
const double GlobalParameters::E_min = GlobalParameters::is_E_log_scale ?  1e-4: 0; // Only in log scale mode
const double GlobalParameters::E_max = 100;    

const int GlobalParameters::s1_nbins = 1000;
const double GlobalParameters::s1_max = 200;
  
const int GlobalParameters::s2_nbins = 1000;
const double GlobalParameters::s2_max = 100000;

const int GlobalParameters::logs2s1_nbins = 1000;
const double GlobalParameters::logs2s1_min = 1;
const double GlobalParameters::logs2s1_max = 4;


  
#endif
