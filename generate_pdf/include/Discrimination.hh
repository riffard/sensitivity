#ifndef Discrimination_hh
#define Discrimination_hh 1


#include "PdfCollection.hh"


void AddFlatDiscrimination(double discrimination_rate, PdfCollection& pdf){

  pdf.hEnergy->Scale(1-discrimination_rate, "nosw2");
}




#endif
