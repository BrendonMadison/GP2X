/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooMyConvolvedBetaPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooMyConvolvedBetaPdf); 

RooMyConvolvedBetaPdf::RooMyConvolvedBetaPdf(const char *name, const char *title, 
                       RooAbsReal& _x,
                       RooAbsReal& _Escale,
                       RooAbsReal& _alpha,
                       RooAbsReal& _beta,
                       RooAbsReal& _sigma) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   Escale("Escale","Escale",this,_Escale),
   alpha("alpha","alpha",this,_alpha),
   beta("beta","beta",this,_beta),
   sigma("sigma","sigma",this,_sigma)
{ 
} 


RooMyConvolvedBetaPdf::RooMyConvolvedBetaPdf(const RooMyConvolvedBetaPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   Escale("Escale",this,other.Escale),
   alpha("alpha",this,other.alpha),
   beta("beta",this,other.beta),
   sigma("sigma",this,other.sigma)
{ 
}


Double_t RooMyConvolvedBetaPdf::prefactor() const
{
  
    double bnorm = TMath::Beta(alpha,beta);
    double prefactden = TMath::Sqrt(2.0*TMath::Pi())*sigma*bnorm;
    return prefactden;
     
}


Double_t RooMyConvolvedBetaPdf::evaluate1() const 
{ 

//    const double eta = 5.0;
    double bnorm = TMath::Beta(alpha,beta);
    double prefact = eta/(TMath::Sqrt(2.0*TMath::Pi())*sigma*bnorm);
   
// Let's do simple Simpson's rule for now in place for the t integral 
// from t = 0.0 to t = 1.0.
    const int N=500;             // make sure this number is even!
    double dt = 1.0/double(N);
    double t;
    double integral = 0.0;
    double xscaled = x/Escale;
    double sigmafrac = 0.01*sigma;
    
    for (int i = 0; i <=N; i++ ) {
        t = double(i)*dt;
        double poly = TMath::Power(t, eta*beta - 1.0) * TMath::Power(1.0 - TMath::Power(t, eta), alpha - 1.0);
        double z = (xscaled + TMath::Power(t, eta) - 1.0)/sigmafrac;
        double gaus = TMath::Exp(-0.5*z*z);
        double tintegrand = prefact*poly*gaus;
        if(i%2==1){
            integral += 4.0*tintegrand; // interior odd abscissa
        }
        else if(i==0 || i== N){         // exterior abscissa
            integral += tintegrand;
        }
        else{
            integral += 2.0*tintegrand; // interior even abscissa      
        }
    }
    
    integral *= dt/3.0;
    return integral;
     
} 


Double_t RooMyConvolvedBetaPdf::evaluate2() const 
{ 

    bool DEBUG = false;   //BEWARE of your disk filling up ...
    double prefactden = RooMyConvolvedBetaPdf::prefactor();
//    const double eta = 5.0;
    double prefact = eta/prefactden;
   
// Let's do simple Simpson's rule for now in place for the t integral 
// from t = 0.0 to t = 1.0.
    const int N=250;             // make sure this number is even!
    double dt = 1.0/double(N);
    double t;
    double integral = 0.0;
    double xscaled = x/Escale;
    double sigmafrac = 0.01*sigma;
    
    for (int i = 0; i <=N; i++ ) {
        t = double(i)*dt;
        double poly = TMath::Power(t, eta*beta - 1.0) * TMath::Power(1.0 - TMath::Power(t, eta), alpha - 1.0);
        double z = (xscaled + TMath::Power(t, eta) - 1.0)/sigmafrac;
        double gaus = TMath::Exp(-0.5*z*z);
        double tintegrand = prefact*poly*gaus;
        if(i%2==1){
            integral += 4.0*tintegrand; // interior odd abscissa (counting from 0: So i=1,3,5...)
        }
        else if(i==0 || i== N){         // exterior abscissa
            integral += tintegrand;
        }
        else{
            integral += 2.0*tintegrand; // interior even abscissa (counting from 0: So i=2,4,6,..))     
        }
    }
    
    integral *= dt/3.0;
    
    if(DEBUG){
        cout << "evaluate2 (Simpson's) N = " << N << " x = " << x << " integral = " << integral << endl;
    }
    
    return integral;
    
}

Double_t RooMyConvolvedBetaPdf::evaluate3() const 
{ 

    double prefactden = RooMyConvolvedBetaPdf::prefactor();
//    const double eta = 5.0;
    double prefact = eta/prefactden;
   
// Let's do Boole's rule for now in place for the t integral 
// from t = 0.0 to t = 1.0.
    const int N=500;             // make sure this number is a multiple of 4.
    double dt = 1.0/double(N);
    double t;
    double integral = 0.0;
    double xscaled = x/Escale;
    double sigmafrac = 0.01*sigma;
    
    for (int i = 0; i <=N; i++ ) {
        t = double(i)*dt;
        double poly = TMath::Power(t, eta*beta - 1.0) * TMath::Power(1.0 - TMath::Power(t, eta), alpha - 1.0);
        double z = (xscaled + TMath::Power(t, eta) - 1.0)/sigmafrac;
        double gaus = TMath::Exp(-0.5*z*z);
        double tintegrand = prefact*poly*gaus;

        if ( i==0 || i==N){
            integral += 7.0*tintegrand;
        }
        else if (i%2==1){            
            integral += 32.0*tintegrand;        
        }    
        else if (i%4==2){
            integral += 12.0*tintegrand;
        }
        else{
            integral += 14.0*tintegrand; // interior even abscissa (counting from 0: So i=2,4,6,..))     
        }
        
    }
    
    integral *= 2.0*dt/45.0;
    return integral;
    
}

Double_t RooMyConvolvedBetaPdf::evaluate4() const 
{ 

    bool DEBUG = false;   //BEWARE of your disk filling up ...
    double prefactden = RooMyConvolvedBetaPdf::prefactor();
//    const double eta = 5.0;
    double prefact = eta/prefactden;
   
// Let's do simple Simpson's rule for now in place for the t integral 
// from t = 0.0 to t = 1.0.
    const int N=250;             // make sure this number is even!
    double t;
    double integral = 0.0;
    double xscaled = x/Escale;
    double sigmafrac = 0.01*sigma;
    
    double tmin = 0.0;
    double tmax = 1.0;
    double dt = (tmax - tmin)/double(N);
    
    for (int i = 0; i <=N; i++ ) {
        t = tmin + double(i)*dt;
        double poly = TMath::Power(t, eta*beta - 1.0) * TMath::Power(1.0 - TMath::Power(t, eta), alpha - 1.0);
        double z = (xscaled + TMath::Power(t, eta) - 1.0)/sigmafrac;
        double gaus = TMath::Exp(-0.5*z*z);
        double tintegrand = prefact*poly*gaus;
        if(i%2==1){
            integral += 4.0*tintegrand; // interior odd abscissa (counting from 0: So i=1,3,5...)
        }
        else if(i==0 || i== N){         // exterior abscissa
            integral += tintegrand;
        }
        else{
            integral += 2.0*tintegrand; // interior even abscissa (counting from 0: So i=2,4,6,..))     
        }
    }
    
    integral *= dt/3.0;
    
    if(DEBUG){
        cout << "evaluate4 (Simpson's) N = " << N << " x = " << x << " integral = " << integral << endl;
    }
    
    return integral;
    
}

Double_t RooMyConvolvedBetaPdf::evaluate5() const 
{ 

    bool DEBUG = false;   //BEWARE of your disk filling up ...
    double prefactden = RooMyConvolvedBetaPdf::prefactor();
//    const double eta = 5.0;
    double prefact = eta/prefactden;
   
// Let's do simple Simpson's rule on range of t including Gaussian +- 6 sigma.

    const int N=100;             // make sure this number is even!
    double t;
    double integral = 0.0;
    double xscaled = x/Escale;
    double sigmafrac = 0.01*sigma;
    
    double tmin = 0.0;
    double tmax = 1.0;
    
    const double SMAX = 6.0;     // Integration range for Gaussian in standard deviation
    double tmineta = (1.0 - xscaled) - SMAX*sigmafrac;
    double tmaxeta = (1.0 - xscaled) + SMAX*sigmafrac;
    if (tmineta > 0.0) tmin = TMath::Power(tmineta, 1.0/eta);
    if (tmaxeta < 1.0) tmax = TMath::Power(tmaxeta, 1.0/eta);
    
    double dt = (tmax-tmin)/double(N);
    
    for (int i = 0; i <=N; i++ ) {
        t = tmin + double(i)*dt;
        double poly = TMath::Power(t, eta*beta - 1.0) * TMath::Power(1.0 - TMath::Power(t, eta), alpha - 1.0);
        double z = (xscaled + TMath::Power(t, eta) - 1.0)/sigmafrac;
        double gaus = TMath::Exp(-0.5*z*z);
        double tintegrand = prefact*poly*gaus;
        if(i%2==1){
            integral += 4.0*tintegrand; // interior odd abscissa (counting from 0: So i=1,3,5...)
        }
        else if(i==0 || i== N){         // exterior abscissa
            integral += tintegrand;
        }
        else{
            integral += 2.0*tintegrand; // interior even abscissa (counting from 0: So i=2,4,6,..))     
        }
    }
    
    integral *= dt/3.0;
    
    if(DEBUG){
        cout << "evaluate5 (Simpson's with N = " << N << "Using t range " << tmin << " " << tmax << " x = " << x << " integral = " << integral << endl;
    }
    
    return integral;
    
}

Double_t RooMyConvolvedBetaPdf::evaluate6() const 
{ 

    bool DEBUG = false;   //BEWARE of your disk filling up ...
    double prefactden = RooMyConvolvedBetaPdf::prefactor();
//    const double eta = 5.0;
    double prefact = eta/prefactden;
   
// Let's do Boole's rule on range of t including Gaussian +- 6 sigma.

    const int N=180;             // make sure this number is divisible by 4!
    double t;
    double integral = 0.0;
    double xscaled = x/Escale;
    double sigmafrac = 0.01*sigma;
    
    double tmin = 0.0;
    double tmax = 1.0;
    
    const double SMAX = 6.0;     // Integration range for Gaussian in standard deviation
    double tmineta = (1.0 - xscaled) - SMAX*sigmafrac;
    double tmaxeta = (1.0 - xscaled) + SMAX*sigmafrac;
    if (tmineta > 0.0) tmin = TMath::Power(tmineta, 1.0/eta);
    if (tmaxeta < 1.0) tmax = TMath::Power(tmaxeta, 1.0/eta);
    
    double dt = (tmax-tmin)/double(N);
    
    for (int i = 0; i <=N; i++ ) {
        t = tmin + double(i)*dt;
        double poly = TMath::Power(t, eta*beta - 1.0) * TMath::Power(1.0 - TMath::Power(t, eta), alpha - 1.0);
        double z = (xscaled + TMath::Power(t, eta) - 1.0)/sigmafrac;
        double gaus = TMath::Exp(-0.5*z*z);
        double tintegrand = prefact*poly*gaus;
        
        if ( i==0 || i==N){
            integral += 7.0*tintegrand;
        }
        else if (i%2==1){            
            integral += 32.0*tintegrand;        
        }    
        else if (i%4==2){
            integral += 12.0*tintegrand;
        }
        else{
            integral += 14.0*tintegrand; // interior even abscissa (counting from 0: So i=2,4,6,..))     
        }                
    }
    
    integral *= 2.0*dt/45.0;
    
    if(DEBUG){
        cout << "evaluate6 (Boole's with N = " << N << "Using t range " << tmin << " " << tmax << " x = " << x << " integral = " << integral << endl;
    }
    
    return integral;
    
}

Double_t RooMyConvolvedBetaPdf::integrand(double t) const{

// Assume for now that t is in range as computed from calling routine

//    double prefactden = RooMyConvolvedBetaPdf::prefactor();
//    double prefact = eta/prefactden;
    
    double xscaled = x/Escale;
    double sigmafrac = 0.01*sigma;

    double poly = TMath::Power(t, eta*beta - 1.0) * TMath::Power(1.0 - TMath::Power(t, eta), alpha - 1.0);
    double z = (xscaled + TMath::Power(t, eta) - 1.0)/sigmafrac;
    double tintegrand = poly*TMath::Exp(-0.5*z*z);

    return tintegrand;
    
}

Double_t RooMyConvolvedBetaPdf::evaluate7() const 
{ 

    bool DEBUG = false;   //BEWARE of your disk filling up if this is set true ...
    double prefactden = RooMyConvolvedBetaPdf::prefactor();
    double prefact = eta/prefactden;
   
    double t;
    double integral = 0.0;
    double xscaled = x/Escale;
    double sigmafrac = 0.01*sigma;
    
    double tmin = 0.0;
    double tmax = 1.0;
    
    const double SMAX = 6.0;     // Integration range for Gaussian in standard deviation
    double tmineta = (1.0 - xscaled) - SMAX*sigmafrac;
    double tmaxeta = (1.0 - xscaled) + SMAX*sigmafrac;
    if (tmineta > 0.0) tmin = TMath::Power(tmineta, 1.0/eta);
    if (tmaxeta < 1.0) tmax = TMath::Power(tmaxeta, 1.0/eta);
    
 // Integrate from [tmin, tmax] using Gaussian or Gauss-Kronrod quadrature.
 
 // First implement 30-point Gaussian quadrature (F). The xgkABCDEF-variable is in [-1,1].
    double c = (tmax-tmin)/2.0;
    double d = (tmin+tmax)/2.0;
    

    bool GAUSS = true;
    if(GAUSS){    
// Gaussian quadrature with 30 points    
        for (int i=0; i<15; ++i){
            int k=2*i+1;
            double pintegral = wgF[i]*(RooMyConvolvedBetaPdf::integrand(-c*xgkF[k] + d) + RooMyConvolvedBetaPdf::integrand(c*xgkF[k] + d));
            integral += pintegral;
        }
/*        
// Gaussian quadrature with 20 points        
        for (int i=0; i<10; ++i){
            int k=2*i+1;
            double pintegral = wgD[i]*(RooMyConvolvedBetaPdf::integrand(-c*xgkD[k] + d) + RooMyConvolvedBetaPdf::integrand(c*xgkD[k] + d));
            integral += pintegral;
        }
*/        
    }    
    else{
// Kronrod quadrature with 41 points
        integral += wgkD[20]*RooMyConvolvedBetaPdf::integrand(c*xgkD[20] + d);
        for (int i=0; i<20; ++i){
             double pintegral = wgkD[i]*(RooMyConvolvedBetaPdf::integrand(-c*xgkD[i] + d) + RooMyConvolvedBetaPdf::integrand(c*xgkD[i] + d));
             integral += pintegral;
        }    
/*
// Kronrod quadrature with 61 points
        integral += wgkF[30]*RooMyConvolvedBetaPdf::integrand(c*xgkF[30] + d);
        for (int i=0; i<30; ++i){
             double pintegral = wgkF[i]*(RooMyConvolvedBetaPdf::integrand(-c*xgkF[i] + d) + RooMyConvolvedBetaPdf::integrand(c*xgkF[i] + d));
             integral += pintegral;
        }
*/        
    }    
    
    integral = c*prefact*integral;
       
// FIX ME. For debugging to be useful for precision need more significant figures in the ouput
// Could also print-out the corresponding values of x
    if(DEBUG){
        if(GAUSS){
            cout << "evaluate7 (Gaussian quadrature with N = 30 Using t range " << tmin << " " << tmax << " x = " << x << " integral = " << integral << endl;
        }
        else{
            cout << "evaluate7 (Kronrod quadrature with N = 41 Using t range " << tmin << " " << tmax << " x = " << x << " integral = " << integral << endl;        
        }
    }
    
    return integral;
    
}

 
Double_t RooMyConvolvedBetaPdf::evaluate() const 
{ 

//    double value = RooMyConvolvedBetaPdf::evaluate1();   //N=500 Simpson's
//    double value = RooMyConvolvedBetaPdf::evaluate2();     //N=250 Simpson's 
//    double value = RooMyConvolvedBetaPdf::evaluate3();   //N=500 Boole's 
//    double value = RooMyConvolvedBetaPdf::evaluate4();   // A better implementation?
    double value = RooMyConvolvedBetaPdf::evaluate7();   // A better implementation?
// I had played previously with different weights
//    value /= (x*x);    // Apply 1/s cross-section weight
//    value *=Escale;    // include Jacobian term?

    return value;
    
} 
