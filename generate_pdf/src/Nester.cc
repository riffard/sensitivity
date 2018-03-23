#include "Nester.hh"
#include "TestSpectra.hh"




//----------------------------------------------------
//----------------------------------------------------
Nester::Nester(string detector_name){

  fDetectorName = detector_name;


  if(fDetectorName ==  "LuxRun3")  InitLuxRun3();
  else if(fDetectorName == "LuxRun4TB1")  InitLuxRun4TB1();
  else if(fDetectorName == "LuxRun4TB2")  InitLuxRun4TB2();
  else if(fDetectorName == "LuxRun4TB3")  InitLuxRun4TB3();
  else if(fDetectorName == "LuxRun4TB4")  InitLuxRun4TB4();
  else {
    cout<<"VMDMVLAMSV"<<endl;
    exit(-1);    
  }


  fNEST = NEST::NESTcalc();

  fSeed = 0;
  fNEST.SetRandomSeed(fSeed);

  frho = -1;
  fdvD = -1;
  
  gasGap_cm = gasGap_mm / 10.;


  
}
//----------------------------------------------------
//----------------------------------------------------


//----------------------------------------------------
//----------------------------------------------------
void Nester::InitLuxRun3(){

  cout<<"<Info::Nester> Load LUX run 3 parameters"<<endl;
  
  // Primary Scintillation (S1) parameters
  g1 = 0.117; //phd per S1 photon in liquid at dtCntr (not phe)
  sPEres = 0.37; //single phe resolution (Gaussian assumed)
  sPEthr = 0.3; //POD threshold in phe
  sPEeff = 1.0; //actual efficiency, can be used in lieu of POD threshold
  //noise[2] = {-0.01,0.08}; //baseline noise mean and width in PE (Gaussian)
  noise[0] = -0.01; //baseline noise mean and width in PE (Gaussian)
  noise[1] = 0.08; //baseline noise mean and width in PE (Gaussian)
  P_dphe = 0.173; //chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT
    
  coinLevel = 2; //how many PMTs have to fire for an S1 to count
  numPMTs = 119; //For coincidence calculation
    
  //S1 PDE quartic polynomial for function of z
  //s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binomial distribution)
  //s1poly[5] = {0.84637, 0.0016222, -9.9224e-6, 4.7508e-8, -7.1692e-11}; // unitless, 1.000 at detector center
  s1poly[0] = 0.84637;
  s1poly[1] = 0.0016222;
  s1poly[2] = -9.9224e-6;
  s1poly[3] = 4.7508e-8;
  s1poly[4] = -7.1692e-11;
  
  //Drift electric field as function of Z in mm
  //The coefficients for a quintic poly, in rising order
  //efpoly[6] = {158.92, -0.22090, 0.0024485, -8.7098e-6, 1.5049e-8, -1.0110e-11}; // in V/cm
  efpoly[0] = 158.92;
  efpoly[1] = -0.22090;
  efpoly[2] = 0.0024485;
  efpoly[3] = -8.7098e-6;
  efpoly[4] = 1.5049e-8;
  efpoly[5] = -1.0110e-11;
  
  // Ionization and Secondary Scintillation (S2) parameters
  g1_gas = 0.102; //phd per S2 photon in gas, used to get SE size
  s2Fano = 1.4; //Fano-like fudge factor for SE width
  s2_thr = 150.; //the S2 threshold in phd. Effects NR most
  S2botTotRatio = 0.449; //S2 bottom-to-total ratio, not really used anymore
  E_gas = 6.148; //field in kV/cm between liquid/gas border and anode
  eLife_us = 800.; //the drift electron mean lifetime in micro-seconds
    
  // Thermodynamic Properties
  T_Kelvin = 176.; //for liquid drift speed calculation
  p_bar = 1.8; //gas pressure in units of bars, it controls S2 size
    
  // Data Analysis Parameters and Geometry
  dtCntr = 160.; //center of detector for S1 corrections, in usec.
  dt_min = 38.; //minimum. Top of detector fiducial volume
  dt_max = 305.; //maximum. Bottom of detector fiducial volume
  liquidBorder = 544.2198; // mm
  gasGap_mm = 5.6; //EL gap in mm, affecting both field and linear S2 term
    
}
//----------------------------------------------------
//----------------------------------------------------

//----------------------------------------------------
//----------------------------------------------------
void Nester::InitLuxRun4TB1(){

  cout<<"<Info::Nester> Load LUX run 4 time-bin 1 parameters"<<endl;

  // Primary Scintillation (S1) parameters
  g1 = 0.1043; //phd per S1 photon in liquid at dtCntr (not phe)
  sPEres = 0.37; //single phe resolution (Gaussian assumed)
  sPEthr = 0.3; //POD threshold in phe
  sPEeff = 1.0; //actual efficiency, can be used in lieu of POD threshold
  //noise[2] = {-0.01,0.08}; //baseline noise mean and width in PE (Gaussian)
  noise[0] = -0.01;
  noise[1] = 0.08;
    
  P_dphe = 0.173; //chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

  coinLevel = 2; //how many PMTs have to fire for an S1 to count
  numPMTs = 118; //For coincidence calculation

  //S1 PDE quartic polynomial for function of z
  //s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binomial distribution)
  //s1poly[5] = {1.1775, -9.9623e-4, 7.7049e-7, 0, 0}; // unitless, 1.000 at detector center
  s1poly[0] = 1.1775;
  s1poly[1] = -9.9623e-4;
  s1poly[2] = 7.7049e-7;
  s1poly[3] = 0;
  s1poly[4] = 0;
    
  //Drift electric field as function of Z in mm
  //The coefficients for a quintic poly, in rising order
  //efpoly[6] = {0, 0, 0, 0, 0, 0}; // in V/cm
  efpoly[0] = 126.120812197;
  efpoly[1] = -1.90341778291;
  efpoly[2] = 0.0221484569642;
  efpoly[3] = -9.87866796872e-5;
  efpoly[4] = 2.0828842974e-7;
  efpoly[5] = -1.57249357722e-10;

  // Ionization and Secondary Scintillation (S2) parameters
  g1_gas = 0.08845; //phd per S2 photon in gas, used to get SE size
  s2Fano = 0.; //Fano-like fudge factor for SE width
  s2_thr = 0.; //the S2 threshold in phd. Effects NR most
  S2botTotRatio = 0.449; //S2 bottom-to-total ratio, not really used anymore
  E_gas = 7.5; //field in kV/cm between liquid/gas border and anode
  eLife_us = 735.; //the drift electron mean lifetime in micro-seconds

  // Thermodynamic Properties
  T_Kelvin = 178.; //for liquid drift speed calculation
  p_bar = 1.95; //gas pressure in units of bars, it controls S2 size

  // Data Analysis Parameters and Geometry
  dtCntr = 162.535; //center of detector for S1 corrections, in usec.
  dt_min = 40.; //minimum. Top of detector fiducial volume
  dt_max = 300.; //maximum. Bottom of detector fiducial volume
  liquidBorder = 544.2198; // mm
  gasGap_mm = 5.027; //EL gap in mm, affecting both field and linear S2 term
 
}
//----------------------------------------------------
//----------------------------------------------------

//----------------------------------------------------
//----------------------------------------------------
void Nester::InitLuxRun4TB2(){

  cout<<"<Info::Nester> Load LUX run 4 time-bin 2 parameters"<<endl;

  // Primary Scintillation (S1) parameters
  g1 = 0.1013; //phd per S1 photon in liquid at dtCntr (not phe)
  sPEres = 0.37; //single phe resolution (Gaussian assumed)
  sPEthr = 0.3; //POD threshold in phe
  sPEeff = 1.0; //actual efficiency, can be used in lieu of POD threshold
  //noise[2] = {-0.01,0.08}; //baseline noise mean and width in PE (Gaussian)
  noise[0] = -0.01;
  noise[1] = 0.08;
  P_dphe = 0.173; //chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

  coinLevel = 2; //how many PMTs have to fire for an S1 to count
  numPMTs = 118; //For coincidence calculation

  //S1 PDE quartic polynomial for function of z
  //s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binomial distribution)
  //s1poly[5] = {1.1775, -9.9623e-4, 7.7049e-7, 0, 0}; // unitless, 1.000 at detector center
  s1poly[0] = 1.1775;
  s1poly[1] = -9.9623e-4;
  s1poly[2] = 7.7049e-7;
  s1poly[3] = 0;
  s1poly[4] = 0;

  //Drift electric field as function of Z in mm
  //The coefficients for a quintic poly, in rising order
  //efpoly[6] = {0, 0, 0, 0, 0, 0}; // in V/cm
  efpoly[0] = 147.431375761;
  efpoly[1] = -3.03129014586;
  efpoly[2] = 0.0343192819835;
  efpoly[3] = -0.000156379508552;
  efpoly[4] = 3.31145682998e-7;
  efpoly[5] = -2.50628583785e-10;

  // Ionization and Secondary Scintillation (S2) parameters
  g1_gas = 0.08708; //phd per S2 photon in gas, used to get SE size
  s2Fano = 0.2; //Fano-like fudge factor for SE width
  s2_thr = 0.; //the S2 threshold in phd. Effects NR most
  S2botTotRatio = 0.449; //S2 bottom-to-total ratio, not really used anymore
  E_gas = 7.73; //field in kV/cm between liquid/gas border and anode
  eLife_us = 947.; //the drift electron mean lifetime in micro-seconds

  // Thermodynamic Properties
  T_Kelvin = 178.; //for liquid drift speed calculation
  p_bar = 1.95; //gas pressure in units of bars, it controls S2 size

  // Data Analysis Parameters and Geometry
  dtCntr = 162.535; //center of detector for S1 corrections, in usec.
  dt_min = 40.; //minimum. Top of detector fiducial volume
  dt_max = 300.; //maximum. Bottom of detector fiducial volume
  liquidBorder = 544.2198; // mm
  gasGap_mm = 5.027; //EL gap in mm, affecting both field and linear S2 term
}
//----------------------------------------------------
//----------------------------------------------------


//----------------------------------------------------
//----------------------------------------------------
void Nester::InitLuxRun4TB3(){

  cout<<"<Info::Nester> Load LUX run 4 time-bin 3 parameters"<<endl;

  // Primary Scintillation (S1) parameters
  g1 = 0.1013; //phd per S1 photon in liquid at dtCntr (not phe)
  sPEres = 0.37; //single phe resolution (Gaussian assumed)
  sPEthr = 0.3; //POD threshold in phe
  sPEeff = 1.0; //actual efficiency, can be used in lieu of POD threshold
  //noise[2] = {-0.01,0.08}; //baseline noise mean and width in PE (Gaussian)
  noise[0] = -0.01;
  noise[1] = 0.08;
  P_dphe = 0.173; //chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

  coinLevel = 2; //how many PMTs have to fire for an S1 to count
  numPMTs = 118; //For coincidence calculation

  //S1 PDE quartic polynomial for function of z
  //s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binomial distribution)
  //s1poly[5] = {1.1775, -9.9623e-4, 7.7049e-7, 0, 0}; // unitless, 1.000 at detector center
  s1poly[0] = 1.1775;
  s1poly[1] = -9.9623e-4;
  s1poly[2] = 7.7049e-7;
  s1poly[3] = 0;
  s1poly[4] = 0;

  //Drift electric field as function of Z in mm
  //The coefficients for a quintic poly, in rising order
  //efpoly[6] = {0, 0, 0, 0, 0, 0}; // in V/cm
  efpoly[0] = 103.47236631;
  efpoly[1] = -1.82104229385;
  efpoly[2] = 0.0222141312792;
  efpoly[3] = -0.000105722325686;
  efpoly[4] = 2.38267089834e-7;
  efpoly[5] = -1.87687319572e-10;

  // Ionization and Secondary Scintillation (S2) parameters
  g1_gas = 0.086; //phd per S2 photon in gas, used to get SE size
  s2Fano = 1.14; //Fano-like fudge factor for SE width
  s2_thr = 0.; //the S2 threshold in phd. Effects NR most
  S2botTotRatio = 0.449; //S2 bottom-to-total ratio, not really used anymore
  E_gas = 7.8; //field in kV/cm between liquid/gas border and anode
  eLife_us = 871.; //the drift electron mean lifetime in micro-seconds

  // Thermodynamic Properties
  T_Kelvin = 178.; //for liquid drift speed calculation
  p_bar = 1.95; //gas pressure in units of bars, it controls S2 size

  // Data Analysis Parameters and Geometry
  dtCntr = 162.535; //center of detector for S1 corrections, in usec.
  dt_min = 40.; //minimum. Top of detector fiducial volume
  dt_max = 300.; //maximum. Bottom of detector fiducial volume
  liquidBorder = 544.2198; // mm
  gasGap_mm = 5.027; //EL gap in mm, affecting both field and linear S2 term


}
//----------------------------------------------------
//----------------------------------------------------


//----------------------------------------------------
//----------------------------------------------------
void Nester::InitLuxRun4TB4(){

  cout<<"<Info::Nester> Load LUX run 4 time-bin 4 parameters"<<endl;

  // Primary Scintillation (S1) parameters
  g1 = 0.099274; //phd per S1 photon in liquid at dtCntr (not phe)
  sPEres = 0.37; //single phe resolution (Gaussian assumed)
  sPEthr = 0.3; //POD threshold in phe
  sPEeff = 1.0; //actual efficiency, can be used in lieu of POD threshold
  //noise[2] = {-0.01,0.08}; //baseline noise mean and width in PE (Gaussian)
  noise[0] = -0.01;
  noise[1] = 0.08;
  P_dphe = 0.173; //chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

  coinLevel = 2; //how many PMTs have to fire for an S1 to count
  numPMTs = 118; //For coincidence calculation

  //S1 PDE quartic polynomial for function of z
  //s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binomial distribution)
  //s1poly[5] = {1.1775, -9.9623e-4, 7.7049e-7, 0, 0}; // unitless, 1.000 at detector center
  s1poly[0] = 1.1775;
  s1poly[1] = -9.9623e-4;
  s1poly[2] = 7.7049e-7;
  s1poly[3] = 0;
  s1poly[4] = 0;

  //Drift electric field as function of Z in mm
  //The coefficients for a quintic poly, in rising order
  //efpoly[6] = {0, 0, 0, 0, 0, 0}; // in V/cm
  efpoly[0] = 212.716223709;
  efpoly[1] = -5.0145653897;
  efpoly[2] = 0.0514514976033;
  efpoly[3] = -0.000224872819366;
  efpoly[4] = 4.6123475598e-7;
  efpoly[5] = -3.42069868433e-10;

  // Ionization and Secondary Scintillation (S2) parameters
  g1_gas = 0.08473; //phd per S2 photon in gas, used to get SE size
  s2Fano = 1.15; //Fano-like fudge factor for SE width
  s2_thr = 0.; //the S2 threshold in phd. Effects NR most
  S2botTotRatio = 0.449; //S2 bottom-to-total ratio, not really used anymore
  E_gas = 7.768; //field in kV/cm between liquid/gas border and anode
  eLife_us = 871.; //the drift electron mean lifetime in micro-seconds

  // Thermodynamic Properties
  T_Kelvin = 178.; //for liquid drift speed calculation
  p_bar = 1.95; //gas pressure in units of bars, it controls S2 size

  // Data Analysis Parameters and Geometry
  dtCntr = 162.535; //center of detector for S1 corrections, in usec.
  dt_min = 40.; //minimum. Top of detector fiducial volume
  dt_max = 300.; //maximum. Bottom of detector fiducial volume
  liquidBorder = 544.2198; // mm
  gasGap_mm = 5.027; //EL gap in mm, affecting both field and linear S2 term

}
//----------------------------------------------------
//----------------------------------------------------


//----------------------------------------------------
//----------------------------------------------------
void Nester::GetNestedPdf(PdfCollection& pdfs, int Nevent, string pt, double field){
  
  INTERACTION_TYPE type_num;
  if (pt =="NR"){
    type_num = NR;
  }else if(pt =="DD"){
    type_num = DD;
  }else if (pt =="CH3T" || pt == "ER"){
    type_num = beta;
  }
  else{
    std::cout << "WARNING YOU HAVE SET AND ILLEGAL PARTICLE TYPE!" << std::endl;
    std::cout << "Please use DD or NR" << std::endl;
    exit(-1);
    
  }

  frho = SetDensity(T_Kelvin);
  fdvD = SetDriftVelocity(T_Kelvin, field);
  //double driftTime = ( liquidBorder - pos_z*10. ) / vD;
  double driftTime = 10; 

  double dt_min = 5;
  double dt_max = 300;
  
  double dz_max = liquidBorder - fdvD * dt_min; // mm - (mm/us)*us = mm
  double dz_min = liquidBorder - fdvD * dt_max; // ditto

  for(int i=0; i<Nevent; i++){

    if(i % 100000 == 0){
      fSeed ++ ;
      fNEST.SetRandomSeed(fSeed);
    }
    
    double energy = 0;
    if(pt =="DD"){

      double eMin = 0;
      double eMax = 100;
      energy = DD_spectrum(eMin, eMax, fNEST);

      pdfs.hEnergy->Fill(energy);
      
    }else energy = pdfs.hEnergy->GetRandom();

    double pos_z = ( dz_min + ( dz_max - dz_min ) * fNEST.rand_uniform() ) * 0.1; //cm
    
    NEST::YieldResult yields = fNEST.GetYields(type_num, energy, frho, field);
    NEST::QuantaResult quanta = fNEST.GetQuanta(yields, frho);
    vector<double> scint_s1 = GetS1(quanta.photons,pos_z);
    vector<double> scint_s2 = GetS2(quanta.electrons,driftTime);
    
    double drift = pos_z/fdvD;
    double S1 = scint_s1[2];
    double S1c = scint_s1[5];
    
    double S2 = scint_s2[4];
    double S2c = scint_s2[7];
    
    pdfs.hs1->Fill(S1c);
    pdfs.hs2->Fill(S2c);
    pdfs.htdrift->Fill(drift);
    pdfs.h2_s1_s2->Fill(S1c, S2c);
    pdfs.h2_s1_logs2s1->Fill(S1c, log10(S2c/S1c));
    
  }

  //Spectra normalization
  
  double Total_rate = pdfs.hEnergy->Integral("width");
  
  pdfs.hs1->Scale(Total_rate / pdfs.hs1->Integral("width"));
  pdfs.hs2->Scale(Total_rate / pdfs.hs2->Integral("width"));
  pdfs.htdrift->Scale(Total_rate / pdfs.htdrift->Integral("width"));
  pdfs.h2_s1_s2->Scale(Total_rate / pdfs.h2_s1_s2->Integral("width"));
  pdfs.h2_s1_logs2s1->Scale(Total_rate / pdfs.h2_s1_logs2s1->Integral("width"));
  
  
}

//----------------------------------------------------
//----------------------------------------------------
