#ifndef NESTER_hh
#define NESTER_hh 1

#include <iostream>

#include "NEST.hh"

#include "PdfCollection.hh"


using namespace std;
using namespace NEST;


class Nester: public NESTcalc{

public:
  Nester(string detector_name);
  ~Nester(){};
  

  void GetNestedPdf(PdfCollection& pdfs, int Nevent, string pt, double field);

    
private:

  void InitLuxRun3();
  void InitLuxRun4TB1();
  void InitLuxRun4TB2();
  void InitLuxRun4TB3();
  void InitLuxRun4TB4();

    
  string fDetectorName;
    
  int fSeed;
  double frho, fdvD, gasGap_cm;
  
  
private:
  
  double SetDriftVelocity ( double Kelvin, double eField ) {
  
    double speed = 0.0; // returns drift speed in mm/usec. based on Fig. 14 arXiv:1712.08607
    int i, j; double vi, vf, slope, Ti, Tf, offset;
  
    double polyExp[11][7] = { { -3.1046, 27.037, -2.1668, 193.27, -4.8024, 646.04, 9.2471 }, //100K
			      { -2.7394, 22.760, -1.7775, 222.72, -5.0836, 724.98, 8.7189 }, //120
			      { -2.3646, 164.91, -1.6984, 21.473, -4.4752, 1202.2, 7.9744 }, //140
			      { -1.8097, 235.65, -1.7621, 36.855, -3.5925, 1356.2, 6.7865 }, //155
			      { -1.5000, 37.021, -1.1430, 6.4590, -4.0337, 855.43, 5.4238 }, //157, merging Miller with Yoo
			      { -1.4939, 47.879, 0.12608, 8.9095, -1.3480, 1310.9, 2.7598 }, //163, merging Miller with Yoo
			      { -1.5389, 26.602, -.44589, 196.08, -1.1516, 1810.8, 2.8912 }, //165
			      { -1.5000, 28.510, -.21948, 183.49, -1.4320, 1652.9, 2.884 }, //167
			      { -1.1781, 49.072, -1.3008, 3438.4, -.14817, 312.12, 2.8049 }, //184
			      {  1.2466, 85.975, -.88005, 918.57, -3.0085, 27.568, 2.3823 }, //200
			      { 334.60 , 37.556, 0.92211, 345.27, -338.00, 37.346, 1.9834 } }; //230
  
    double Temperatures[11] = { 100., 120., 140., 155., 157., 163., 165., 167., 184., 200., 230. };

    i=0;
    if ( Kelvin >= Temperatures[0] && Kelvin < Temperatures[1] ) i = 0;
    else if ( Kelvin >= Temperatures[1] && Kelvin < Temperatures[2] ) i = 1;
    else if ( Kelvin >= Temperatures[2] && Kelvin < Temperatures[3] ) i = 2;
    else if ( Kelvin >= Temperatures[3] && Kelvin < Temperatures[4] ) i = 3;
    else if ( Kelvin >= Temperatures[4] && Kelvin < Temperatures[5] ) i = 4;
    else if ( Kelvin >= Temperatures[5] && Kelvin < Temperatures[6] ) i = 5;
    else if ( Kelvin >= Temperatures[6] && Kelvin < Temperatures[7] ) i = 6;
    else if ( Kelvin >= Temperatures[7] && Kelvin < Temperatures[8] ) i = 7;
    else if ( Kelvin >= Temperatures[8] && Kelvin < Temperatures[9] ) i = 8;
    else if ( Kelvin >= Temperatures[9] && Kelvin <= Temperatures[10] ) i = 9;
    else {
      cout << "\nERROR: TEMPERATURE OUT OF RANGE (100-230 K)\n";
    }
  
    j = i + 1;
    Ti = Temperatures[i];
    Tf = Temperatures[j];
    // functional form from http://zunzun.com
    vi = polyExp[i][0]*exp(-eField/polyExp[i][1])+polyExp[i][2]*exp(-eField/polyExp[i][3])+polyExp[i][4]*exp(-eField/polyExp[i][5])+polyExp[i][6];
    vf = polyExp[j][0]*exp(-eField/polyExp[j][1])+polyExp[j][2]*exp(-eField/polyExp[j][3])+polyExp[j][4]*exp(-eField/polyExp[j][5])+polyExp[j][6];
    if ( Kelvin == Ti ) return vi;
    if ( Kelvin == Tf ) return vf;
    if ( vf < vi ) {
      offset = (sqrt((Tf*(vf-vi)-Ti*(vf-vi)-4.)*(vf-vi))+sqrt(Tf-Ti)*(vf+vi))/(2.*sqrt(Tf-Ti));
      slope = -(sqrt(Tf-Ti)*sqrt((Tf*(vf-vi)-Ti*(vf-vi)-4.)*(vf-vi))-(Tf+Ti)*(vf-vi))/(2.*(vf-vi));
      speed = 1. / ( Kelvin - slope ) + offset;
    }
    else {
      slope = ( vf - vi ) / ( Tf - Ti );
      speed = slope * ( Kelvin - Ti ) + vi;
    }
  
    return speed;
  
  }

  double SetDensity ( double Kelvin ) { // currently only for fixed pressure (saturated vapor pressure); will add pressure dependence later
  
    if ( Kelvin < 161.40 ) // solid Xenon
      return 3.41; // from Yoo at 157K; other sources say 3.100 (Wikipedia, 'max') and 3.64 g/mL at unknown T's
  
    return 
      2.9970938084691329E+02 * exp ( -8.2598864714323525E-02 * Kelvin ) - 1.8801286589442915E+06 * exp ( - pow ( ( Kelvin - 4.0820251276172212E+02 ) / 2.7863170223154846E+01, 2. ) )
      - 5.4964506351743057E+03 * exp ( - pow ( ( Kelvin - 6.3688597345042672E+02 ) / 1.1225818853661815E+02, 2. ) )
      + 8.3450538370682614E+02 * exp ( - pow ( ( Kelvin + 4.8840568924597342E+01 ) / 7.3804147172071107E+03, 2. ) )
      - 8.3086310405942265E+02; // in grams per cubic centimeter based on zunzun fit to NIST data; will add gas later
  
  }


  
  
public:

  int coinLevel, numPMTs;
  double g1, sPEres, sPEthr = 0.3, sPEeff, noise[2], P_dphe, s1poly[5], efpoly[6], g1_gas, s2Fano, s2_thr, S2botTotRatio, E_gas, eLife_us, T_Kelvin, p_bar, dtCntr, dt_min, dt_max, liquidBorder, gasGap_mm;
    
  

};



#endif
