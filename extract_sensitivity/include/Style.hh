#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"


void LoadStyle(){
  
  gROOT->SetStyle("Modern");
  // color palette
  const Int_t _NRGBs = 5;
  const Int_t _NCont = 255;
  Double_t stops[ _NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[ _NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[ _NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[ _NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(_NRGBs, stops, red, green, blue, _NCont);
  gStyle->SetNumberContours(_NCont);

  
  gROOT->GetColor(10)->SetRGB(216/255.,0/255.,23/255.); // Red
  gROOT->GetColor(11)->SetRGB(42/255.,103/255.,168/255.); // blue
  gROOT->GetColor(12)->SetRGB(64/255.,160/255.,58/255.); // green 
  gROOT->GetColor(13)->SetRGB(131/255.,52/255.,146/255.); // purpel
  gROOT->GetColor(14)->SetRGB(247/255.,106/255.,10/255.); // orange
  gROOT->GetColor(15)->SetRGB(251/255.,195/255.,13/255.);  // yellow
  gROOT->GetColor(16)->SetRGB(147/255.,67/255.,31/255.); // brown
  gROOT->GetColor(17)->SetRGB(238/255.,104/255.,177/255.); // pink
  gROOT->GetColor(18)->SetRGB(40/255.,172/255.,142/255.);
  gROOT->GetColor(19)->SetRGB(77/255.,77/255.,77/255.);

  gROOT->GetColor(20)->SetRGB( 67/255., 148/255., 246/255.); // blue block
  gROOT->GetColor(21)->SetRGB( 17/255., 78/255., 179/255.);
  gROOT->GetColor(22)->SetRGB( 19/255., 60/255., 115/255.);
  gROOT->GetColor(23)->SetRGB( 4/255., 25/255., 65/255.);

  gROOT->GetColor(24)->SetRGB( 95/255., 182/255., 50/255.); // green block
  gROOT->GetColor(25)->SetRGB( 28/255., 121/255., 32/255.);
  gROOT->GetColor(26)->SetRGB( 18/255., 77/255., 19/255.);
  gROOT->GetColor(27)->SetRGB( 12/255., 51/255., 10/255.);

  gROOT->GetColor(28)->SetRGB( 241/255., 204/255., 32/255.); // yellow block
  gROOT->GetColor(29)->SetRGB( 211/255., 178/255., 28/255.);
  gROOT->GetColor(30)->SetRGB( 181/255., 135/255., 21/255.);
  gROOT->GetColor(31)->SetRGB( 145/255., 98/255., 15/255.);

  gROOT->GetColor(32)->SetRGB( 238/255., 124/255., 21/255.); // orange block
  gROOT->GetColor(33)->SetRGB( 212/255., 84/255., 16/255.);
  gROOT->GetColor(34)->SetRGB( 174/255., 71/255., 12/255.);
  gROOT->GetColor(35)->SetRGB( 127/255., 53/255., 10/255.);

  gROOT->GetColor(36)->SetRGB( 229/255., 69/255., 69/255.); // red
  gROOT->GetColor(37)->SetRGB( 186/255., 18/255., 9/255.);
  gROOT->GetColor(38)->SetRGB( 112/255., 4/255., 3/255.);
  gROOT->GetColor(39)->SetRGB( 68/255., 1/255., 8/255.);

  gROOT->GetColor(40)->SetRGB( 162/255., 77/255., 218/255.); // violet
  gROOT->GetColor(41)->SetRGB( 98/255., 41/255., 138/255.);
  gROOT->GetColor(42)->SetRGB( 75/255., 32/255., 105/255.);
  gROOT->GetColor(43)->SetRGB( 44/255., 21/255., 61/255.);

  // Extra colors
  gROOT->GetColor(44)->SetRGB( 232/255., 54/255., 56/255.);
  gROOT->GetColor(45)->SetRGB( 52/255., 132/255., 136/255.);
  gROOT->GetColor(46)->SetRGB( 150/255., 208/255., 52/255.);
  gROOT->GetColor(47)->SetRGB( 112/255., 166/255., 22/255.);
  gROOT->GetColor(48)->SetRGB( 245/255., 82/255., 31/255.);
  gROOT->GetColor(49)->SetRGB( 50/255., 212/255., 214/255.);
  
  
  //titles
  gStyle->SetTitleSize(0.04, "xyz") ; 
  gStyle->SetTitleOffset(1.05, "xz") ; 
  gStyle->SetTitleOffset(1.1, "y") ; 
  gStyle->SetTextFont(62) ; 
  
  //stat
  gStyle->SetOptStat("eIMR") ;
  
  //legend
  gStyle->SetLegendFillColor(0) ; 
  gStyle->SetLegendFillColor(0) ; 
  
	
  gROOT->ForceStyle() ;
  
  
}

