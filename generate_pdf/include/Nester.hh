#ifndef NESTER_hh
#define NESTER_hh 1

#include <iostream>

#include "NEST.hh"

#include "PdfCollection.hh"


using namespace std;
using namespace NEST;


class Nester{

public:
  Nester(string detector_name);
  ~Nester(){};
  

  void GetNestedPdf(PdfCollection& pdfs, int Nevent, string pt, double field);

  
  NEST::NESTcalc* GetNEST(){return &fNEST;};
  
private:

  void InitLuxRun3();
  void InitLuxRun4TB1();
  void InitLuxRun4TB2();
  void InitLuxRun4TB3();
  void InitLuxRun4TB4();

    
  string fDetectorName;
  
  NEST::NESTcalc fNEST;
  
  int fSeed;
  double frho, fdvD, gasGap_cm;
  
  //----------------------------------------------------
  //----------------------------------------------------
  // NEST function
  //----------------------------------------------------
  //----------------------------------------------------
  vector<double> GetS1 ( int Nph, double dz ) {
  
    vector<double> scintillation(8);  // return vector
  
    // Add some variability in g1 drawn from a polynomial spline fit
    double posDep = s1poly[0] + s1poly[1] * dz +
      s1poly[2] * pow(dz,2.)+
      s1poly[3] * pow(dz,3.)+
      s1poly[4] * pow(dz,4.);
    double dz_center = liquidBorder - fdvD * dtCntr; //go from t to z
    posDep /= s1poly[0]+
      s1poly[1] * pow(dz_center/10.,1.)+
      s1poly[2] * pow(dz_center/10.,2.)+
      s1poly[3] * pow(dz_center/10.,3.)+
      s1poly[4] * pow(dz_center/10.,4.); //factor 10 for mm to cm conversion
  
    // generate a number of PMT hits drawn from a binomial distribution. Initialize number of photo-electrons
    int nHits=fNEST.BinomFluct(Nph,g1*posDep), Nphe = 0; if ( nHits >= numPMTs ) { sPEeff = 1.0; sPEthr = 0.0; }
  
    // Initialize the pulse area and spike count variables
    double pulseArea = 0., spike = 0., prob;
  
    // If single photo-electron efficiency is under 1 and the threshold is above 0 (some phe will be below threshold)
    if ( sPEthr > 0. ) {
      // Step through the pmt hits
      for ( int i = 0; i < nHits; i++ ) {
	// generate photo electron, integer count and area
	double phe1 = fNEST.rand_gauss(1.,sPEres) + fNEST.rand_gauss(noise[0],noise[1]); Nphe++;
	prob = fNEST.rand_uniform();
	// zero the area if random draw determines it wouldn't have been observed.
	if ( prob > sPEeff ) { phe1 = 0.; } //add an else with Nphe++ if not doing mc truth
	// Generate a double photo electron if random draw allows it
	double phe2 = 0.;
	if ( fNEST.rand_uniform() < P_dphe ) {
	  // generate area and increment the photo-electron counter
	  phe2 = fNEST.rand_gauss(1.,sPEres) + fNEST.rand_gauss(noise[0],noise[1]); Nphe++;
	  // zero the area if phe wouldn't have been observed
	  if ( fNEST.rand_uniform() > sPEeff && prob > sPEeff ) { phe2 = 0.; } //add an else with Nphe++ if not doing mc truth
	  // The dphe occurs simultaneously to the first one from the same source photon. If the first one is seen, so should be the second one
	}
	// Save the phe area and increment the spike count (very perfect spike count) if area is above threshold
	if ( (phe1+phe2) > sPEthr ) { spike++; pulseArea += phe1 + phe2; }
      }
    }
    else { // apply just an empirical efficiency by itself, without direct area threshold
      Nphe = nHits + fNEST.BinomFluct(nHits,P_dphe);
      pulseArea = fNEST.rand_gauss(fNEST.BinomFluct(Nphe,1.-(1.-sPEeff)/(1.+P_dphe)),sPEres*sqrt(Nphe));
      spike = (double)nHits;
    }
    if ( pulseArea < 0. ) pulseArea = 0.;
    double pulseAreaC= pulseArea / posDep;
    double Nphd = pulseArea / (1.+P_dphe);
    double NphdC= pulseAreaC/ (1.+P_dphe);
    double spikeC = spike / posDep;
  
    scintillation[0] = nHits; scintillation[1] = Nphe;
    scintillation[2] = pulseArea; scintillation[3] = pulseAreaC;
    scintillation[4] = Nphd; scintillation[5] = NphdC;
    scintillation[6] = spike; scintillation[7] = spikeC;
  
    if ( spike > 10 )
      prob = 1.;
    else {
      if ( spike >= coinLevel ) {

	double numer = 0;
	double denom = 0;

	for ( int i = spike; i > 0; i-- ) {
	  denom += nCr ( numPMTs, i );
	  if ( i >= coinLevel ) numer += nCr ( numPMTs, i );
	}
	prob = numer / denom;
      }
      else
	prob = 0.;
    }
  
    if ( fNEST.rand_uniform() < prob ) // coincidence has to happen in different PMTs
      { ; }
    else { // some of these are set to -1 to flag them as having been below threshold
      //scintillation[0] *= -1.;
      //scintillation[1] *= -1.;
      scintillation[2] *= -1.;
      scintillation[3] *= -1.;
      //scintillation[4] *= -1.;
      //scintillation[5] *= -1.;
      scintillation[6] *= -1.;
      scintillation[7] *= -1.;
    }
  
    return scintillation;
  
  }

  vector<double> GetS2 ( int Ne, double dt ) {
  
    vector<double> ionization(8);
    double alpha = 0.137, beta = 177., gamma = 45.7;
    double epsilon = 1.85 / 1.00126;
  
    double E_liq = E_gas / epsilon; //kV per cm
    double ExtEff = -0.03754*pow(E_liq,2.)+0.52660*E_liq-0.84645; // arXiv:1710.11032
    if ( ExtEff > 1. ) ExtEff = 1.;
    if ( ExtEff < 0. ) ExtEff = 0.;
    int Nee = fNEST.BinomFluct(Ne,ExtEff*exp(-dt/eLife_us));
  
    double elYield = Nee* (alpha*E_gas*1000.-beta*p_bar-gamma)* gasGap_cm; // arXiv:1207.2292
    int Nph = int(floor(fNEST.rand_gauss(elYield,sqrt(s2Fano*elYield))+0.5));
    int nHits = fNEST.BinomFluct(Nph,g1_gas);
    int Nphe = nHits + fNEST.BinomFluct(nHits,P_dphe);
    double pulseArea=fNEST.rand_gauss(Nphe,sPEres*sqrt(Nphe));
    double pulseAreaC= pulseArea/exp(-dt/eLife_us);
    double Nphd = pulseArea / (1.+P_dphe);
    double NphdC= pulseAreaC/ (1.+P_dphe);
  
    double S2b = fNEST.rand_gauss(S2botTotRatio*pulseArea,sqrt(S2botTotRatio*pulseArea*(1.-S2botTotRatio)));
    double S2bc= S2b / exp(-dt/eLife_us); // for detectors using S2 bottom-only in their analyses
  
    ionization[0] = Nee; ionization[1] = Nph;
    ionization[2] = nHits; ionization[3] = Nphe;
    if ( s2_thr >= 0 ) {
      ionization[4] = pulseArea; ionization[5] = pulseAreaC;
      ionization[6] = Nphd; ionization[7] = NphdC;
    }
    else {
      ionization[4] = S2b; ionization[5] = S2bc;
      ionization[6] = S2b / (1.+P_dphe); ionization[7] = S2bc / (1.+P_dphe);
    }
  
    return ionization;
  
  }


  
private:

  long double Factorial ( double x ) {
  
    return tgammal ( x + 1. );
  
  }

  double nCr ( double n, double r ) {
  
    return Factorial(n) /
      ( Factorial(r) * Factorial(n-r) );
  
  }

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
