
// Primary Scintillation (S1) parameters
double g1 = 0.1013; //phd per S1 photon in liquid at dtCntr (not phe)
double sPEres = 0.37; //single phe resolution (Gaussian assumed)
double sPEthr = 0.3; //POD threshold in phe
double sPEeff = 1.0; //actual efficiency, can be used in lieu of POD threshold
double noise[2] = {-0.01,0.08}; //baseline noise mean and width in PE (Gaussian)
double P_dphe = 0.173; //chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

int coinLevel = 2; //how many PMTs have to fire for an S1 to count
int numPMTs = 118; //For coincidence calculation

//S1 PDE quartic polynomial for function of z
//s1polA + s1polB*z[mm] + s1polC*z^2+... (QE included, for binomial distribution)
double s1poly[5] = {1.1775, -9.9623e-4, 7.7049e-7, 0, 0}; // unitless, 1.000 at detector center

//Drift electric field as function of Z in mm
//The coefficients for a quintic poly, in rising order
double efpoly[6] = {0, 0, 0, 0, 0, 0}; // in V/cm

 // Ionization and Secondary Scintillation (S2) parameters
double g1_gas = 0.086; //phd per S2 photon in gas, used to get SE size
double s2Fano = 1.14; //Fano-like fudge factor for SE width
double s2_thr = 0.; //the S2 threshold in phd. Effects NR most
double S2botTotRatio = 0.449; //S2 bottom-to-total ratio, not really used anymore
double E_gas = 7.8; //field in kV/cm between liquid/gas border and anode
double eLife_us = 871.; //the drift electron mean lifetime in micro-seconds

 // Thermodynamic Properties
double T_Kelvin = 178.; //for liquid drift speed calculation
double p_bar = 1.95; //gas pressure in units of bars, it controls S2 size

 // Data Analysis Parameters and Geometry
double dtCntr = 162.535; //center of detector for S1 corrections, in usec.
double dt_min = 40.; //minimum. Top of detector fiducial volume
double dt_max = 300.; //maximum. Bottom of detector fiducial volume
double liquidBorder = 544.2198; // mm
double gasGap_mm = 5.027; //EL gap in mm, affecting both field and linear S2 term
