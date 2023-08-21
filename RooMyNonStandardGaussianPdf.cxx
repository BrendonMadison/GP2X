/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooMyNonStandardGaussianPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooMyNonStandardGaussianPdf); 

 RooMyNonStandardGaussianPdf::RooMyNonStandardGaussianPdf(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mean,
                        RooAbsReal& _sigmapc) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mean("mean","mean",this,_mean),
   sigmapc("sigmapc","sigmapc",this,_sigmapc)
 { 
 } 


 RooMyNonStandardGaussianPdf::RooMyNonStandardGaussianPdf(const RooMyNonStandardGaussianPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mean("mean",this,other.mean),
   sigmapc("sigmapc",this,other.sigmapc)
 { 
 } 



 Double_t RooMyNonStandardGaussianPdf::evaluate() const 
 { 
   // ENTER EXPRESSION IN TERMS OF VARIABLE ARGUMENTS HERE 
   double rms = 0.01*sigmapc*mean;
   double z = (x-mean)/rms;
   double value = (1.0/rms)*TMath::Exp(-0.5*z*z);
   return value;
 } 
 
