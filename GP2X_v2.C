#include <cmath>
//So that we can use round else we get some weird results later on

void vBoost(TLorentzVector *FV, Double_t Px, Double_t Py, Double_t Pz, Double_t Ene)
{
  Double_t rx = atanh(Px/Ene);
  Double_t ry = atanh(Py/Ene);
  Double_t rz = atanh(Pz/Ene);
  Double_t nPx = sinh(rx)*FV->E() + cosh(rx)*FV->Px();
  Double_t nPy = sinh(ry)*FV->E() + cosh(ry)*FV->Py();
  Double_t nPz = sinh(rz)*FV->E() + cosh(rz)*FV->Pz();
  Double_t nEn = cosh(rx+ry+rz)*FV->E() + sinh(rx)*FV->Px() + sinh(ry)*FV->Py() + sinh(rz)*FV->Pz();
  FV->SetXYZT(nPx,nPy,nPz,nEn);
}

//Returns a gaussian random variate using the analytical solution
//I've put this in here because I've learnt that the standard
//ROOT implementation using TRandom->Gaus() is actually quite bad
Double_t RanGaus(Double_t MU, Double_t VAR, Double_t RUNI)
{
  return TMath::ErfInverse(RUNI * 2.0 - 1.0) * sqrt(2.0)*VAR + MU;
}

//Struct  is used to call and pass the KKMC (or BHWIDE) 
//fit and simulation values.
struct FinalFourVectors {
  TLorentzVector *Mu;
  TLorentzVector *AMu;
  TLorentzVector *Gam;
  TLorentzVector *DiMu;
  TLorentzVector *GetDiMu()
  {
    TLorentzVector *OutDiMu = new TLorentzVector();
    OutDiMu->SetXYZM(Mu->Px()+AMu->Px(),Mu->Py()+AMu->Py(),Mu->Pz()+AMu->Pz(),sqrt((Mu->E() + AMu->E()) * (Mu->E() + AMu->E()) - (Mu->Px()+AMu->Px())*(Mu->Px()+AMu->Px()) - (Mu->Py()+AMu->Py())*(Mu->Py()+AMu->Py()) - (Mu->Pz()+AMu->Pz())*(Mu->Pz()+AMu->Pz())));
    return OutDiMu;
  }
};

//Sets that the name for the variable is FFV
typedef struct FinalFourVectors FinalFourVec;

//Struct FileValues is used to call and pass the KKMC (or BHWIDE) 
//fit and simulation values.
struct FileValues {
  Double_t DELENE, SX0, SX1, SX2, SE0, SE1, ELO, EHI, ENOM1, ENOM2, NUM, PXMEAN, PXSTD, PYMEAN, PYSTD, PZMEAN, PZSTD;
  void PrintValues()
  {
    cout << "\n----------PRINTING FIT VALUES----------" << endl;
    cout << "Energy bins of KKMC(BHWIDE) [GeV] : " << DELENE << endl;
    cout << "X-section fit parameters [0]+[1]*x+[2]*x^2 : \n\t" << SX0 << "\n\t" << SX1 << "\n\t" << SX2 << endl;
    cout << "Entry$ fit parameters [0]+[1]*x : \n\t" << SE0 << "\n\t" << SE1 << endl;
    cout << "Lowest energy bin [GeV] : " << ELO << endl;
    cout << "Highest energy bin [GeV] : " << EHI << endl;
    cout << "Mode energy of E1 [GeV] : " << ENOM1 << endl;
    cout << "Mode energy of E2 [GeV] : " << ENOM2 << endl;
    cout << "Number of events per KKMC(BHWIDE) run : " << NUM << endl;
    cout << "GuineaPig Position Gaussian Fit Values : " << "\nX-MEAN = " << PXMEAN << "\tX-STD = " << PXSTD << endl;
    cout << "Y-MEAN = " << PYMEAN << "\tY-STD = " << PYSTD << "\nZ-MEAN = " << PZMEAN << "\tZ-STD = " << PZSTD << endl;
  }
};

//Sets that the name for the variable is "FtVal"
typedef struct FileValues FtVal;

//Struct to use for getting smeared tracker Four Vector
//The parameterized tracker resolution and smearing code is based off of 
//Graham Wilson's https://github.com/grahamwwilson/MyLHEReader readGenericWLHE.py as of Feb. 21st 2023
struct ParticleTrack {
  TLorentzVector *PartTrack;
  TLorentzVector *OutTrack = new TLorentzVector();
  TRandom3 *randSmear = new TRandom3();
  TTimeStamp *smearTime = new TTimeStamp();
  ParticleTrack(){smearTime->Set(); randSmear->SetSeed(smearTime->GetNanoSec());}

  TLorentzVector *SmearECAL()
  {
    //# Electromagnetic calorimeter smearing of energy 
    OutTrack->SetXYZM(PartTrack->Px(),PartTrack->Py(),PartTrack->Pz(),PartTrack->M());

    //Here we establish a minimum energy of 0.1 GeV.
    Double_t MinEne = 0.1;
    //z  =  rg.Gaus(0.0,1.0)
    //E = v.E()
    //c1 = 0.18
    //c2 = 0.01
    //varE = c1*c1*E + c2*c2*E*E
    //sigE = math.sqrt(varE)
    //Enew = E + sigE*z
    Double_t smearIt = 0.0;
    Double_t Enew = -1.0;
    //Smear the new energy
    //smearTime->Set();
    //randSmear->SetSeed(smearTime->GetNanoSec());
    //smearIt = RanGaus(0.0,1.0,randSmear->Uniform());
    smearIt = randSmear->Gaus(0.0,1.0);
    Enew = OutTrack->E() + smearIt*(sqrt(0.18*0.18*OutTrack->E() + 0.01*0.01*OutTrack->E() * OutTrack->E()));
    //If we fail the minum energy threshold then it is undetected, "0.0"
    if (Enew < MinEne)
      {
	Enew = 0.0;
      }

    //Computes the energy rescaling in the massless case
    //And if it doesn't pass that check does it in the massive case
    //
    Double_t KE = Enew/(OutTrack->E());
    if (OutTrack->M() > 0.0005)
      {
	KE = sqrt((Enew*Enew - OutTrack->M() * OutTrack->M()) /(OutTrack->E() * OutTrack->E() - OutTrack->M() * OutTrack->M()));
      }
    //k = Enew/E
    OutTrack->SetXYZT(KE * OutTrack->Px(),KE * OutTrack->Py(),KE * OutTrack->Pz(),Enew);
    //vsmear = rt.TLorentzVector(v.Px()*k,v.Py()*k,v.Pz()*k,Enew)
    //#    vsmear +=v
    //#    vsmear +=v
    return OutTrack;
  }
  TLorentzVector *SmearTrack()
  {
    //TLorentzVector *OutTrack = new TLorentzVector();
    OutTrack->SetXYZM(PartTrack->Px(),PartTrack->Py(),PartTrack->Pz(),PartTrack->M());
    //OutDiMu->SetXYZM(Mu->Px()+AMu->Px(),Mu->Py()+AMu->Py(),Mu->Pz()+AMu->Pz(),sqrt((Mu->E() + AMu->E()) * (Mu->E() + AMu->E()) - (Mu->Px()+AMu->Px())*(Mu->Px()+AMu->Px()) - (Mu->Py()+AMu->Py())*(Mu->Py()+AMu->Py()) - (Mu->Pz()+AMu->Pz())*(Mu->Pz()+AMu->Pz())));
    //Get the angle in degrees
    Double_t thetad = 180.0*OutTrack->Theta()/TMath::Pi();
    if (thetad > 90.0)
      {
	thetad = 180.0 - thetad;
      }
    //Get the a and b smearing parameters
    Double_t a = 1.0;
    Double_t b = 1.0;

    if (thetad>85.0)
      {
	a = 2.09e-5;
	b = 0.978e-3;
      }

    if (thetad>40.0 and thetad<=85.0)
      {
	// linear interpolation between 40 and 85 degrees
	//aparh = 2.09e-5;
	//bparh = 0.978e-3;
	//aparl = 1.96e-5;
	//bparl = 0.755e-3;
	//thdmin = 40.0;
	//dtheta = 85.0-40.0;
	//a = aparl + (aparh-aparl)*(thetad-thdmin)/dtheta;
	//b = bparl + (bparh-bparl)*(thetad-thdmin)/dtheta;
	a = 1.96e-5 + (2.09e-5 - 1.96e-5)*(thetad-40.0)/(85.0-40.0);
	b = 0.755e-3 + (0.978e-3 - 0.755e-3)*(thetad-40.0)/(85.0-40.0);  
      }
  
    if (thetad>37.573305 and thetad<=40.0)
      {
	// Use 40 deg parameters - as still in barrel part of TPC
	a = 1.96e-5;
	b = 0.755e-3;
      }

    if (thetad>20.0 and thetad<=37.573305)      
      {  
	//  linear interpolation between 20 and 37.57 degrees
	//aparh = 1.96e-5
	//bparh = 0.755e-3
	//aparl = 1.54e-4
	//bparl = 1.48e-3
	//thdmin = 20.0
	//dtheta = 37.573305-20.0
	//rh = 1.808
	// rl (=0.8553)
	//rl = 2.35*tan(20.0*pie/180.0);
	Double_t r = 2.35*tan(thetad*TMath::Pi()/180.0);
	//drsq = rh*rh - rl*rl;
	//a = aparl + (aparh-aparl)*(r*r-rl*rl)/drsq;
	//b = bparl + (bparh-bparl)*(thetad-thdmin)/dtheta;
	a = 1.54e-4 + (1.96e-5 - 1.54e-4)*(r*r - 0.8553*0.8553)/(1.808*1.808 - 0.8553*0.8553);
	b = 1.48e-3 + (0.755e-3 - 1.48e-3)*(thetad - 20.0)/(37.573305-20.0);
      }
  
    if(thetad>7.0 and thetad<=20.0)
      {
	//linear interpolation between 7 and 20 degrees
	//aparh = 1.54e-4
	//bparh = 1.48e-3
	//aparl = 1.61e-3
	//bparl = 1.29e-3
	//thdmin = 7.0
	//dtheta = 20.0-7.0
	// rh (=0.8553)
	//rh = 2.35*tan(20.0*pie/180.0);
	//rl = 2.35*tan(7.0*pie/180.0);
	//rl = 0.2885
	Double_t r = 2.35*tan(thetad*TMath::Pi()/180.0);
	//drsq = rh*rh - rl*rl;
	a = 1.61e-3 + (1.54e-4 - 1.61e-3)*(r*r - 0.2885*0.2885)/(0.8553*0.8553 - 0.2885*0.2885);
	b = 1.29e-3 + (1.48e-3 - 1.29e-3)*(thetad-7.0)/(13.0);
      }
    //You now have the a and b values so you can smear things!
    //Again, these are originally in Graham Wilson's code preserved to best ability as it was originally python
    // Tracker smearing of inverse pT following ILD LOI. Using theta dependent parametrization.
    // See invptvariance in zgamma.f for details.
    
    //Double_t smearIt = RanGaus(0.0,1.0,randSmear->Uniform());
    Double_t smearIt = randSmear->Gaus(0.0,1.0);
    //mass = v.M()
    //Theta = v.Theta()
    //Phi = v.Phi()
    //pTinv = 1.0/v.Perp()
    //a = 2.0e-5
    //b = 0.5e-3
    //a, b = TrackerResolutionParameters(Theta)
    //Double_t varpTinv = a*a + b*b/(OutTrack->Perp() * OutTrack->Perp() * sin(OutTrack->Theta())*sin(OutTrack->Theta()));
    //Double_t sigpTinv = sqrt(varpTinv);
    //Double_t pTinvnew = pTinv + sigpTinv*z
    //Double_t pT = 1.0/pTinvnew
    //
    //Compute the new smeared values!
    //
    //cout << "Smear values: " << endl;
    //cout << thetad << endl;
    //cout << smearIt << endl;
    Double_t pT = 1.0/(1.0/OutTrack->Perp() + sqrt(a*a + b*b/(OutTrack->Perp() * OutTrack->Perp() * sin(OutTrack->Theta())*sin(OutTrack->Theta())))*smearIt);
    Double_t px = pT*cos(OutTrack->Phi());
    Double_t py = pT*sin(OutTrack->Phi());
    Double_t pz = pT/tan(OutTrack->Theta());
    Double_t E = sqrt(px*px + py*py + pz*pz + OutTrack->M()*OutTrack->M());
    //cout << pT << endl;
    //Write them to the output track and then return it
    OutTrack->SetXYZT(px,py,pz,E);
    //Delete these or else memory overflow!
    return OutTrack;
  }
};

//Sets that the name for the variable is FFV
typedef struct ParticleTrack PartTrk;


FtVal GetGFitVals(char *gname)
{
  //The FtVal variable used to return at the end of this
  FtVal fvg;
  // Read in the KKMC (or BHWIDE) file and then get the cross-section and CoM
  TFile *f = new TFile(gname);
  TTree *t = (TTree*)f->Get("t");
  Double_t Px,Py,Pz;
  t->SetBranchAddress("PosX",&Px);
  t->SetBranchAddress("PosY",&Py);
  t->SetBranchAddress("PosZ",&Pz);

  Int_t n = t->Draw("PosX>>hPX","Entry$ < 100001","goff");
  TH1 *hPX = (TH1*)gDirectory->Get("hPX");
  hPX->Fit("gaus");
  TF1 *fPX = (TF1*)hPX->GetFunction("gaus");
  fvg.PXMEAN = fPX->GetParameter(1);
  fvg.PXSTD = fPX->GetParameter(2);

  Int_t m = t->Draw("PosY>>hPY","Entry$ < 100001","goff");
  TH1 *hPY = (TH1*)gDirectory->Get("hPY");
  hPY->Fit("gaus");
  TF1 *fPY = (TF1*)hPY->GetFunction("gaus");
  fvg.PYMEAN = fPY->GetParameter(1);
  fvg.PYSTD = fPY->GetParameter(2);

  Int_t k = t->Draw("PosZ>>hPZ","Entry$ < 100001","goff");
  TH1 *hPZ = (TH1*)gDirectory->Get("hPZ");
  hPZ->Fit("gaus");
  TF1 *fPZ = (TF1*)hPZ->GetFunction("gaus");
  fvg.PZMEAN = fPZ->GetParameter(1);
  fvg.PZSTD = fPZ->GetParameter(2);
  f->Close();
  return fvg;
}

//Function for filling FtVal with values from kname (KKMC or BHWIDE file)
FtVal GetKFitVals(char *kname)
{
  //The FtVal variable used to return at the end of this
  FtVal fv;
  // Read in the KKMC (or BHWIDE) file and then get the cross-section and CoM
  TFile *f = new TFile(kname);
  TTree *t = (TTree*)f->Get("t");
  Float_t sqs, xsec;
  t->SetBranchAddress("sqs",&sqs);
  t->SetBranchAddress("xsec",&xsec);

  //Get the low and high energies
  //We start with the low energies and essentially "zoom in" on it till we get to it at high resolution
  t->Draw("sqs>>hs(1000,0,10000)","Entry$ % 100 == 0","goff");
  TH1 *hs = (TH1*) gDirectory->Get("hs");
  fv.ELO = round(1000.0*hs->GetBinLowEdge(hs->FindFirstBinAbove(2.0)))/1000.0;
  //zoom to resolution of 0.04 GeV
  t->Draw(Form("sqs>>ht(1000,%f,%f)",fv.ELO-20.0,fv.ELO+20.0),"Entry$ % 100 == 0","goff");
  TH1 *ht = (TH1*) gDirectory->Get("ht");
  fv.ELO = round(1000.0*ht->GetBinLowEdge(ht->FindFirstBinAbove(2.0)))/1000.0;
  //zoom to resolution of 0.0005 GeV
  t->Draw(Form("sqs>>hu(1000,%f,%f)",fv.ELO-0.25,fv.ELO+0.25),"Entry$ % 100 == 0","goff");
  TH1 *hu = (TH1*) gDirectory->Get("hu");
  fv.ELO = round(1000.0*hu->GetBinCenter(hu->FindFirstBinAbove(2.0)))/1000.0;

  //We repeat for the high energies
  t->Draw("sqs>>hv(1000,0,10000)","Entry$ % 100 == 0","goff");
  TH1 *hv = (TH1*) gDirectory->Get("hv");
  fv.EHI = round(1000.0*hv->GetBinCenter(hv->FindLastBinAbove(2.0)))/1000.0;
  //zoom to resolution of 0.04 GeV
  t->Draw(Form("sqs>>hw(1000,%f,%f)",fv.EHI-20.0,fv.EHI+20.0),"Entry$ % 100 == 0","goff");
  TH1 *hw = (TH1*) gDirectory->Get("hw");
  fv.EHI = round(1000.0*hw->GetBinCenter(hw->FindLastBinAbove(2.0)))/1000.0;
  //zoom to resolution of 0.0005 GeV
  t->Draw(Form("sqs>>hx(1000,%f,%f)",fv.EHI-0.25,fv.EHI+0.25),"","goff");
  TH1 *hx = (TH1*) gDirectory->Get("hx");
  fv.EHI = round(1000.0*hx->GetBinCenter(hx->FindLastBinAbove(2.0)))/1000.0;
  //Get number of events per bin (aka KKMC or BHWIDE run)
  fv.NUM = hx->GetBinContent(hx->FindFirstBinAbove(2.0));

  //Get the energy increment by taking the remainder with a given energy in the range
  //For this to work you need a long enough range to cover remainder beyond -0.2
  //And for more than 10 events mod 100 (1000 events)
  //Or to adjust this code respectively
  //t->Draw(Form("remainder(sqs,%f)>>ha(1000,-0.2,-0.1001)",0.5*(fv.EHI+fv.ELO)),"Entry$ % 100 == 0","goff");
  //TH1 *ha = (TH1*) gDirectory->Get("ha");
  //t->Draw(Form("remainder(sqs,%f)>>hl(1000,-0.1,-0.0001)",0.5*(fv.EHI+fv.ELO)),"Entry$ % 100 == 0","goff");
  //TH1 *hl = (TH1*) gDirectory->Get("hl");
  //New method as the old method failed when sqs was rounded in odd ways by the event generators
  //t->Draw(Form("sqs>>hb(1000,%f,%f)",fv.ELO,fv.EHI),"Entry$ % 100 == 0","goff");
  //TH1D *hb = (TH1D*) gDirectory->Get("hb");
  //Int_t firstbin = 0;
  //firstbin = hb->FindFirstBinAbove(10,1,0,-1);
  //cout << "Energy Increment: " << setprecision(3) << 0.5*abs(hl->GetBinCenter(hl->FindFirstBinAbove(10.0)) - ha->GetBinCenter(ha->FindLastBinAbove(10.0))) << endl;
  //fv.DELENE = round(1000.0*0.5*abs(hl->GetBinCenter(hl->FindFirstBinAbove(10.0)) - ha->GetBinCenter(ha->FindLastBinAbove(10.0))))/1000.0;
  //fv.DELENE = hb->GetBinCenter(hb->FindFirstBinAbove(10,1,firstbin+1,-1)) - hb->GetBinCenter(firstbin);
  fv.DELENE = 0.1;

  //The polynomial fit for the cross-section as it is approximately parabolic (more so asymptotic)
  //Plot and save the plot for diagnostic purposes
  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  Int_t k = t->Draw("sqs:xsec>>SXGraph","Entry$ % 100 == 0","");
  TGraph *gSX = new TGraph(k,t->GetV1(),t->GetV2());
  //gSX->Draw("ap");
  //cout << "\n---------------FIT VALUES FOR CROSS-SECTION---------------" << endl;
  gSX->Fit("pol2");

  TF1 *fSX=(TF1*)gSX->GetFunction("pol2");
  Double_t SXp0=fSX->GetParameter(0);
  Double_t SXp1=fSX->GetParameter(1);
  Double_t SXp2=fSX->GetParameter(2);
  fSX->Draw("same");
  c1->SaveAs("XSecDiag.pdf");
  //cout << SXp0 << endl;
  //cout << SXp1 << endl;
  //cout << SXp2 << endl;
  fv.SX0 = SXp0;
  fv.SX1 = SXp1;
  fv.SX2 = SXp2;
  //The linear fit for the Entry vs CoM. This should be linear unless you mucked it up pretty badly
  Int_t m = t->Draw("sqs:Entry$>>SEGraph","Entry$ % 100 == 0","");
  TGraph *gSE = new TGraph(m,t->GetV1(),t->GetV2());
  //gSE->Draw("ap");
  //cout << "\n---------------FIT VALUES FOR ENTRY NUMBERS---------------" << endl;
  gSE->Fit("pol1");

  TF1 *fSE=(TF1*)gSE->GetFunction("pol1");
  Double_t SEp0=fSE->GetParameter(0);
  Double_t SEp1=fSE->GetParameter(1);
  fSE->Draw("same");
  c1->SaveAs("EntryDiag.pdf");

  //cout << SEp0 << endl;
  //cout << SEp1 << endl;
  fv.SE0 = SEp0;
  fv.SE1 = SEp1;
  f->Close();
  delete c1;
  delete gSE;
  delete gSX;
  return fv;
}

FtVal GetFitVals(char *gname, char *kname)
{
  //Does both the fitvals at once!
  FtVal fvg;
  FtVal fvOut;
  //FtVal from GuineaPig
  fvg = GetGFitVals(gname);
  //FtVal from KKMC (or BHWIDE)
  fvOut = GetKFitVals(kname);
  //Copy over the relevant GuineaPig Values
  fvOut.PXMEAN = fvg.PXMEAN;
  fvOut.PXSTD = fvg.PXSTD;
  fvOut.PYMEAN = fvg.PYMEAN;
  fvOut.PYSTD = fvg.PYSTD;
  fvOut.PZMEAN = fvg.PZMEAN;
  fvOut.PZSTD = fvg.PZSTD;
  //Get the nominal (mode) energy value
  TFile *f = new TFile(gname);
  TTree *t = (TTree*)f->Get("t");
  Double_t E1,E2;
  t->SetBranchAddress("E1",&E1);
  t->SetBranchAddress("E2",&E2);

  t->Draw(Form("E1>>hE1(250,%f,%f)",0.5*fvOut.ELO,0.5*fvOut.EHI),"","goff");
  TH1 *hE1 = (TH1*) gDirectory->Get("hE1");
  fvOut.ENOM1 = round(1000.0*hE1->GetBinCenter(hE1->GetMaximumBin()))/1000.0;

  t->Draw(Form("E2>>hE2(250,%f,%f)",0.5*fvOut.ELO,0.5*fvOut.EHI),"","goff");
  TH1 *hE2 = (TH1*) gDirectory->Get("hE2");
  fvOut.ENOM2 = round(1000.0*hE2->GetBinCenter(hE2->GetMaximumBin()))/1000.0;
  f->Close();
  return fvOut;
}

// Function for generating a GuineaPig Lumifile with the desired center of mass energy
void GPGenerate(char *gpname, Double_t CoM)
{
  Double_t E1, E2, x1, y1, x2, y2, PosZ, PosY, PosX;
  //Generate the tree
  TFile *fg = new TFile(gpname,"RECREATE");
  TTree *tg = new TTree("t","GuineaPig Tree");
  tg->Branch("E1",&E1,"E1/D");
  tg->Branch("E2",&E2,"E2/D");
  tg->Branch("x1",&x1,"x1/D");
  tg->Branch("y1",&y1,"y1/D");
  tg->Branch("x2",&x2,"x2/D");
  tg->Branch("y2",&y2,"y2/D");
  tg->Branch("PosX",&PosX,"PosX/D");
  tg->Branch("PosY",&PosY,"PosY/D");
  tg->Branch("PosZ",&PosZ,"PosZ/D");

  //Now we generate 10k events with sqrt(s) == CoM , E1 = CoM/2 , E2 = CoM/2 , 95% of momentum along z, 5% along x,y
  TRandom3 *rand1 = new TRandom3();
  TRandom3 *rand2 = new TRandom3();
  TRandom3 *rand3 = new TRandom3();
  TTimeStamp *tstamp = new TTimeStamp();
  tstamp->Set();
  rand1->SetSeed(tstamp->GetNanoSec());
  tstamp->Set();
  rand2->SetSeed(tstamp->GetNanoSec());
  tstamp->Set();
  rand3->SetSeed(tstamp->GetNanoSec());

  for (Int_t i = 0; i < 100000; i++)
    {
      //We randomly generate x1, y1, x2, y2 as being gaussian with 0.1% of the energy
      //And positions as rms of 100 nm
      E1 = CoM/2.0;
      E2 = CoM/2.0;

      //tstamp->Set();
      //rand1->SetSeed(tstamp->GetNanoSec());
      //tstamp->Set();
      //rand2->SetSeed(tstamp->GetNanoSec());

      //x1 = RanGaus(0.0,0.001,rand1->Uniform());
      //y1 = RanGaus(0.0,0.001,rand2->Uniform());
      x1 = 0.0;
      y1 = 0.0;
      //x1 = rand1->Gaus(0.0,0.001);
      //y1 = rand2->Gaus(0.0,0.001);

      //tstamp->Set();
      //rand1->SetSeed(tstamp->GetNanoSec());
      //tstamp->Set();
      //rand2->SetSeed(tstamp->GetNanoSec());

      //x2 = RanGaus(0.0,0.001,rand1->Uniform());
      //y2 = RanGaus(0.0,0.001,rand2->Uniform());
      //x2 = rand1->Gaus(0.0,0.001);
      //y2 = rand2->Gaus(0.0,0.001);
      x2 = 0.0;
      y2 = 0.0;

      //tstamp->Set();
      //rand1->SetSeed(tstamp->GetNanoSec());
      //tstamp->Set();
      //rand2->SetSeed(tstamp->GetNanoSec());
      //tstamp->Set();
      //rand3->SetSeed(tstamp->GetNanoSec());

      //PosX = RanGaus(0.0,100.0,rand1->Uniform());
      //PosY = RanGaus(0.0,100.0,rand2->Uniform());
      //PosZ = RanGaus(0.0,100.0,rand3->Uniform());
      PosX = rand1->Gaus(0.0,100.0);
      PosY = rand2->Gaus(0.0,100.0);
      PosZ = rand3->Gaus(0.0,100.0);
      tg->Fill();
    }
  fg->Write();
  fg->Close();
}

void GP2X_v2(Int_t ngen,Double_t BES1,Double_t BES2, char *fname, char *gpname, char *kkname)
{
  //Set precision of print statements
  cout.precision(17);
  //Setup the GP input file
  TFile *f = new TFile(gpname);
  TTree *t1 = (TTree*)f->Get("t");
  //Variables we use to read in GP values
  Double_t E1, E2, x1, y1, x2, y2;
  Double_t PosZ, PosY, PosX;
  Double_t gpsqs,px1,px2,py1,py2,pz1,pz2,gamboost,bsqs;
  Double_t bpx, bpy, bpz, ben;
  Double_t nE1,nE2;
  Int_t ev, ri, ent_lo, ent_hi;
  //PHYSICS CONSTANTS
  Double_t PI = 3.1415926;
  Double_t MMU = 0.1056583755;
  Double_t ME = 0.00051099895;
  Double_t cang = 0.014;
  //Fill the FtVal structure for use later
  FtVal FITVALS = GetFitVals(gpname,kkname);
  FITVALS.PrintValues();
  //FOR THE CROSS-SECTION FROM xsec:sqs
  //Double_t XSEC_MAX = FITVALS.SX2*(FITVALS.ELO-7.0*FITVALS.DELENE)*(FITVALS.ELO-7.0*FITVALS.DELENE) + FITVALS.SX1*(FITVALS.ELO-7.0*FITVALS.DELENE) + FITVALS.SX0;
  Double_t XSEC_MAX = FITVALS.SX2*(FITVALS.ELO*0.98)*(FITVALS.ELO*0.98) + FITVALS.SX1*(FITVALS.ELO*0.98) + FITVALS.SX0;
  cout << "MAXIMUM CROSS-SECTION FOR SIMULATION: " << XSEC_MAX << endl;
  //THIS IS THE COUNT OF HOW MANY EVENTS YOU'VE MADE
  //IT IS IN CAPS BECAUSE IT IS IMPORTANT
  Int_t COUNT = 0;
  //Random number generators and timestamps to reseed them
  TRandom3 *rand1 = new TRandom3();
  TRandom3 *rand2 = new TRandom3();
  TRandom3 *rand3 = new TRandom3();
  TRandom3 *rand4 = new TRandom3();
  TRandom3 *rand5 = new TRandom3();
  TTimeStamp *tstamp = new TTimeStamp();
  TTimeStamp *codestamp = new TTimeStamp();
  //TVector3 *GPBoost = new TVector3();
  TVector3 GPBoost;
  TVector3 *CangBoost = new TVector3();
  tstamp->Set();
  rand1->SetSeed(tstamp->GetNanoSec());
  tstamp->Set();
  rand2->SetSeed(tstamp->GetNanoSec());
  tstamp->Set();
  rand3->SetSeed(tstamp->GetNanoSec());
  tstamp->Set();
  rand4->SetSeed(tstamp->GetNanoSec());
  tstamp->Set();
  rand5->SetSeed(tstamp->GetNanoSec());
  //Simulation values , initialize gpsqs
  gpsqs = 999999.9;
  t1->SetBranchAddress("E1",&E1);
  t1->SetBranchAddress("E2",&E2);
  t1->SetBranchAddress("x1",&x1);
  t1->SetBranchAddress("x2",&x2);
  t1->SetBranchAddress("y1",&y1);
  t1->SetBranchAddress("y2",&y2);
  t1->SetBranchAddress("PosX",&PosX);
  t1->SetBranchAddress("PosY",&PosY);
  t1->SetBranchAddress("PosZ",&PosZ);
  
  //read all entries and fill the histograms
  Long64_t nentries = t1->GetEntries();
  Int_t nent = t1->GetEntries();
  t1->Draw(">>+elistgp",Form("2.0*sqrt(E1*E2) > %f && 2.0*sqrt(E1*E2) <= %f && Entry$ < 100001",FITVALS.ELO,FITVALS.EHI),"goff");
  TEventList *elistgp = (TEventList*)gDirectory->Get("elistgp");
  Int_t nlistgp = elistgp->GetN();

  //Set the KKMC input file and its tree
  //TFile *f2 = new TFile("../../KKMC/KKMCResults/RunPolRunP/KKMC_0.3_0.8.root");
  TFile *f2 = new TFile(kkname);
  TTree *t2 = (TTree*)f2->Get("t");
  TLorentzVector *Mup4 = new TLorentzVector();
  TLorentzVector *AMup4 = new TLorentzVector();
  TLorentzVector *GamTotp4 = new TLorentzVector();
  Float_t sqs, xsec, randomroll, maxxsec;
  Int_t ev2, ri2;
  t2->SetBranchAddress("Mup4.",&Mup4);
  t2->SetBranchAddress("AMup4.",&AMup4);
  t2->SetBranchAddress("GamTotp4.",&GamTotp4);
  t2->SetBranchAddress("sqs",&sqs);
  t2->SetBranchAddress("xsec",&xsec);

  maxxsec = t2->GetMaximum("xsec");
  cout << "Maximum cross-section: " << t2->GetMaximum("xsec") << endl;
  cout << "Minimum cross-section: " << t2->GetMinimum("xsec") << endl;

  //Variables for reading all entries and filling the histograms
  Long64_t nentries2 = t2->GetEntries();
  Int_t nent2 = t2->GetEntries();
  bool DiscardFlag = true;
  bool SqsFlag = true;
  // Variables for computing the output
  Double_t sqr, sqrm, truecm, oldgpsqs, inilooptime, finlooptime, tmpBES;
  Double_t K1,K1g,FK2,FpK2,K2,FK3,FpK3,K3, OPosX, OPosY, OPosZ;
  TLorentzVector *DiMup4 = new TLorentzVector();
  TLorentzVector *BoostMup4 = new TLorentzVector();
  TLorentzVector *BoostAMup4 = new TLorentzVector();
  TLorentzVector *BoostGamTotp4 = new TLorentzVector();
  TLorentzVector *BoostDiMup4 = new TLorentzVector();
  TLorentzVector *BoostElep4 = new TLorentzVector();
  TLorentzVector *BoostPosp4 = new TLorentzVector();
  TLorentzVector *SmearMup4 = new TLorentzVector();
  TLorentzVector *SmearAMup4 = new TLorentzVector();
  TLorentzVector *SmearDiMup4 = new TLorentzVector();
  TLorentzVector *SmearGamTotp4 = new TLorentzVector();
  //Track objects
  PartTrk MuTrk;
  PartTrk AMuTrk;
  PartTrk GamTrk;

  //Output file and tree
  TFile *fout = new TFile(fname,"RECREATE");
  TTree *t3 = new TTree("t","t");
  t3->Branch("Elep4","TLorentzVector",&BoostElep4);
  t3->Branch("Posp4","TLorentzVector",&BoostPosp4);
  t3->Branch("Mup4","TLorentzVector",&BoostMup4);
  t3->Branch("AMup4","TLorentzVector",&BoostAMup4);
  t3->Branch("GamTotp4","TLorentzVector",&BoostGamTotp4);
  t3->Branch("DiMup4","TLorentzVector",&BoostDiMup4);
  t3->Branch("SMup4","TLorentzVector",&SmearMup4);
  t3->Branch("SAMup4","TLorentzVector",&SmearAMup4);
  t3->Branch("SDiMup4","TLorentzVector",&SmearDiMup4);
  t3->Branch("SGamTotp4","TLorentzVector",&SmearGamTotp4);
  t3->Branch("TrueCM",&truecm, "TrueCM/D");
  t3->Branch("PosX",&OPosX, "PosX/D");
  t3->Branch("PosY",&OPosY, "PosY/D");
  t3->Branch("PosZ",&OPosZ, "PosZ/D");
  //Begin the loop to generate events now that everything is initialized etc.
  while (COUNT < ngen){
    gpsqs = 9999.9;
    //Randomly get a GP event that is in the energy range accessible by the KKMC events...
    //In a perfect setup there would be KKMC events over the entire GP range
    ri = rand1->Integer(nlistgp);
    //ri = RanIntR(nlistgp,rand1,tstamp);
    t1->GetEntry(elistgp->GetEntry(ri));
    //Comput the CoM of GP event
    px1 = x1*sqrt(E1*E1 - ME*ME);
    px2 = -1.0*x2*sqrt(E2*E2 - ME*ME);
    py1 = y1*sqrt(E1*E1 - ME*ME);
    py2 = -1.0*y2*sqrt(E2*E2 - ME*ME);
    pz1 = sqrt((E1*E1 - ME*ME)*(1.0 - x1*x1 - y1*y1));
    pz2 = -1.0*sqrt((E2*E2 - ME*ME)*(1.0 - x2*x2 - y2*y2));
    oldgpsqs = sqrt(2.0*(E1*E2 - px1*px2 - py1*py2 - pz1*pz2));
    //Roll for beam energy spread and then rescale the initial energies
    while (gpsqs < FITVALS.ELO || gpsqs > FITVALS.EHI)
      {
	//tstamp->Set();
	//rand4->SetSeed(tstamp->GetNanoSec());
	//tstamp->Set();
	//rand5->SetSeed(tstamp->GetNanoSec());
	nE1 = E1 + rand4->Gaus(0.0,BES1);
	nE2 = E2 + rand5->Gaus(0.0,BES2);
	//nE1 = E1 + RanGaus(0.0,BES1,rand4->Uniform());
	//nE2 = E2 + RanGaus(0.0,BES2,rand5->Uniform());
	px1 = x1*sqrt(nE1*nE1 - ME*ME);
	px2 = -1.0*x2*sqrt(nE2*nE2 - ME*ME);
	py1 = y1*sqrt(nE1*nE1 - ME*ME);
	py2 = -1.0*y2*sqrt(nE2*nE2 - ME*ME);
	pz1 = sqrt((nE1*nE1 - ME*ME)*(1.0 - x1*x1 - y1*y1));
	pz2 = -1.0*sqrt((nE2*nE2 - ME*ME)*(1.0 - x2*x2 - y2*y2));
	gpsqs = sqrt(2.0*(nE1*nE2 - px1*px2 - py1*py2 - pz1*pz2));
      }
    //Store these values in the Ele and Pos 4vectors
    BoostElep4->SetXYZT(px1,py1,pz1,nE1);
    BoostPosp4->SetXYZT(px2,py2,pz2,nE2);
    ben = nE1 + nE2;
    bpx = (px1 + px2);
    bpy = (py1 + py2);
    bpz = (pz1 + pz2);
    codestamp->Set();
    //Get the net momenta of the GP system
    //These need to be beta values and a 3-vector
    //gamboost = sqrt(1 + ((px1+px2)*(px1+px2) + (py1+py2)*(py1+py2) + (pz1+pz2)*(pz1+pz2))/(ME*ME) );
    GPBoost.SetXYZ((px1+px2)/(nE1+nE2),(py1+py2)/(nE1+nE2),(pz1+pz2)/(nE1+nE2));
    //Not currently fully supported, but get the crossing angle boost
    //CangBoost->SetXYZ(gpsqs*tan(cang/2.0)/(E1+E2),0.0,0.0);
    //Getting a random tree entry with CoM in acceptable range
    ent_lo = FITVALS.SE1*(gpsqs - FITVALS.DELENE) + FITVALS.SE0;
    ent_hi = FITVALS.SE1*(gpsqs + FITVALS.DELENE) + FITVALS.SE0;
    sqs = 9999.0;
    DiscardFlag = false;
    SqsFlag = true;
    //sometimes this loop gets stuck because it rarely randomly rolls outside its usable range! 
    //so we are going to fix that using a time out
    codestamp->Set();
    inilooptime = codestamp->AsDouble();
    codestamp->Set();
    finlooptime = codestamp->AsDouble() - inilooptime;
    //cout << "TEST1" << endl;
    while (DiscardFlag == false && finlooptime < 1.0)
      {
	codestamp->Set();
	finlooptime = codestamp->AsDouble() - inilooptime;
	//tstamp->Set();
	//rand2->SetSeed(tstamp->GetNanoSec());
	//cout << tstamp->GetSec() << endl;
	ri2 = ent_lo + rand2->Integer(ent_hi - ent_lo);
	t2->GetEntry(ri2);
	//check if we are in an acceptable range to keep
	if (gpsqs - sqs < 1.0*FITVALS.DELENE && gpsqs - sqs >= -1.0*FITVALS.DELENE)
	  {
	//tstamp->Set();
	    //cout << tstamp->GetSec() << endl;
	    //rand3->SetSeed(tstamp->GetNanoSec());
	    randomroll = (rand3->Uniform(1.0));
	    //roll to see if we discard for x-section reasons
	    //cout << randomroll << endl;
	    //cout << gpsqs << endl;
	    //cout << (FITVALS.SX2*gpsqs*gpsqs + FITVALS.SX1*gpsqs + FITVALS.SX0)/XSEC_MAX << endl;
	    if (randomroll <= (FITVALS.SX2*gpsqs*gpsqs + FITVALS.SX1*gpsqs + FITVALS.SX0)/XSEC_MAX)
	      {
		DiscardFlag = true;
	      }
	    else
	      {
		break;
	      }
	      }
      }
    //If we time out let the user know and break the while loop
    if (finlooptime > 1.0)
      {
	cout << "Event timed out: " << gpsqs << endl;
        DiscardFlag = false;
      }
    //Proceed from here if you have a usable event to calculate final state values
    if (DiscardFlag == true)
      {
	if (COUNT >= ngen)
	  {
	    break;
	  }
	//Possibly use the FinalFourVec struct to handle the final state four vectors in the future?
	//Most of these computations are degenerate so it would be nice to have it look nice
	//FinalFourVec ffv;
	//ffv.Mu = (TLorentzVector*)Mup4->Clone("Mu");
	//ffv.AMu = (TLorentzVector*)AMup4->Clone("AMu");
	//ffv.Gam = (TLorentzVector*)GamTotp4->Clone("Gam");
	//ffv.DiMu = ffv.GetDiMu();
	//Reseed and randomly blur the position values so they are more unique
	//cout << "TEST2" << endl;
	//tstamp->Set();
	//cout << tstamp->GetSec() << endl;
	//rand4->SetSeed(tstamp->GetNanoSec());
	OPosX = PosX + rand4->Gaus(0.0,FITVALS.PXSTD/10.0);
	//OPosX = PosX + RanGaus(0.0,FITVALS.PXSTD/10.0,rand4->Uniform());
	//tstamp->Set();
	//rand4->SetSeed(tstamp->GetNanoSec());
	OPosY = PosY + rand4->Gaus(0.0,FITVALS.PYSTD/10.0);
	//OPosY = PosY + RanGaus(0.0,FITVALS.PYSTD/10.0,rand4->Uniform());
	//tstamp->Set();
	//rand4->SetSeed(tstamp->GetNanoSec());
	OPosZ = PosZ + rand4->Gaus(0.0,FITVALS.PZSTD/10.0);
	//OPosZ = PosZ + RanGaus(0.0,FITVALS.PZSTD/10.0,rand4->Uniform());
	//Rescale muon momenta according to the rescaling relationship ... K1g is in the massless limit
	//K1g = gpsqs/sqs;
	//K1 = K1g;
	BoostMup4->SetXYZM(Mup4->Px(),Mup4->Py(),Mup4->Pz(),MMU);
	BoostAMup4->SetXYZM(AMup4->Px(),AMup4->Py(),AMup4->Pz(),MMU);
	BoostGamTotp4->SetXYZT(GamTotp4->Px(),GamTotp4->Py(),GamTotp4->Pz(),GamTotp4->E());
	K1 = sqrt(((gpsqs*gpsqs) - 2.0*MMU*MMU + GamTotp4->M() * GamTotp4->M())/(Mup4->P()*Mup4->P() + AMup4->P()*AMup4->P() + GamTotp4->P() * GamTotp4->P())) / (sqrt(((sqs*sqs) - 2.0*MMU*MMU + GamTotp4->M() * GamTotp4->M())/(Mup4->P()*Mup4->P() + AMup4->P()*AMup4->P() + GamTotp4->P() * GamTotp4->P())));
	BoostMup4->SetXYZM(K1*BoostMup4->Px(),K1*BoostMup4->Py(),K1*BoostMup4->Pz(),MMU);
	BoostAMup4->SetXYZM(K1*BoostAMup4->Px(),K1*BoostAMup4->Py(),K1*BoostAMup4->Pz(),MMU);
	BoostGamTotp4->SetXYZT(K1*BoostGamTotp4->Px(),K1*BoostGamTotp4->Py(),K1*BoostGamTotp4->Pz(),K1*BoostGamTotp4->E()); 
	//In the massive case we use K1
	vBoost(BoostMup4,bpx,bpy,bpz,ben);
	vBoost(BoostAMup4,bpx,bpy,bpz,ben);
	vBoost(BoostGamTotp4,bpx,bpy,bpz,ben);
	BoostDiMup4->SetXYZM(BoostMup4->Px()+BoostAMup4->Px(),BoostMup4->Py()+BoostAMup4->Py(),BoostMup4->Pz()+BoostAMup4->Pz(),sqrt((BoostMup4->E() + BoostAMup4->E()) * (BoostMup4->E() + BoostAMup4->E()) - (BoostMup4->Px()+BoostAMup4->Px())*(BoostMup4->Px()+BoostAMup4->Px()) - (BoostMup4->Py()+BoostAMup4->Py())*(BoostMup4->Py()+BoostAMup4->Py()) - (BoostMup4->Pz()+BoostAMup4->Pz())*(BoostMup4->Pz()+BoostAMup4->Pz())));
	//BoostDiMup4->SetXYZT(BoostMup4->Px()+BoostAMup4->Px(),BoostMup4->Py()+BoostAMup4->Py(),BoostMup4->Pz()+BoostAMup4->Pz(),BoostMup4->E() + BoostAMup4->E());
	//bsqs = sqrt((BoostDiMup4->E() + BoostGamTotp4->E())*(BoostDiMup4->E() + BoostGamTotp4->E()) - (BoostDiMup4->Px() + BoostGamTotp4->Px())*(BoostDiMup4->Px() + BoostGamTotp4->Px()) - (BoostDiMup4->Py() + BoostGamTotp4->Py())*(BoostDiMup4->Py() + BoostGamTotp4->Py()) - (BoostDiMup4->Pz() + BoostGamTotp4->Pz())*(BoostDiMup4->Pz() + BoostGamTotp4->Pz()));
	//BoostMup4->Boost((px1+px2)/(nE1+nE2),(py1+py2)/(nE1+nE2),(pz1+pz2)/(nE1+nE2));
	//BoostAMup4->Boost((px1+px2)/(nE1+nE2),(py1+py2)/(nE1+nE2),(pz1+pz2)/(nE1+nE2));
	//BoostGamTotp4->Boost((px1+px2)/(nE1+nE2),(py1+py2)/(nE1+nE2),(pz1+pz2)/(nE1+nE2));
	//K1 = sqrt(((gpsqs*gpsqs) - 2.0*MMU*MMU + GamTotp4->M() * GamTotp4->M())/(Mup4->P()*Mup4->P() + AMup4->P()*AMup4->P() + GamTotp4->P() * GamTotp4->P())) / (sqrt(((sqs*sqs) - 2.0*MMU*MMU + GamTotp4->M() * GamTotp4->M())/(Mup4->P()*Mup4->P() + AMup4->P()*AMup4->P() + GamTotp4->P() * GamTotp4->P())));
	//K1 = bsqs/sqs;
	//K1 = sqrt(((bsqs*bsqs) - 2.0*MMU*MMU + BoostGamTotp4->M() * BoostGamTotp4->M())/(BoostMup4->P()*BoostMup4->P() + BoostAMup4->P()*BoostAMup4->P() + BoostGamTotp4->P() * BoostGamTotp4->P())) / (sqrt(((sqs*sqs) - 2.0*MMU*MMU + BoostGamTotp4->M() * BoostGamTotp4->M())/(BoostMup4->P()*BoostMup4->P() + BoostAMup4->P()*BoostAMup4->P() + BoostGamTotp4->P() * BoostGamTotp4->P())));
	//BoostElep4->SetXYZM(K1g*px1,K1g*py1,K1g*pz1,ME);
	//BoostPosp4->SetXYZM(K1g*px2,K1g*py2,K1g*pz2,ME);
	//BoostMup4->SetXYZM(K1*BoostMup4->Px(),K1*BoostMup4->Py(),K1*BoostMup4->Pz(),MMU);
	//BoostMup4->SetXYZT(K1*Mup4->Px(),K1*Mup4->Py(),K1*Mup4->Pz(),K1*Mup4->E());
	//BoostAMup4->SetXYZM(K1*BoostAMup4->Px(),K1*BoostAMup4->Py(),K1*BoostAMup4->Pz(),MMU);
	//BoostMup4->Boost(*CangBoost);
	//BoostAMup4->SetXYZM(K1*AMup4->Px(),K1*AMup4->Py(),K1*AMup4->Pz(),MMU);
	//BoostAMup4->Boost(*CangBoost);
	//BoostGamTotp4->SetXYZT(K1*BoostGamTotp4->Px(),K1*BoostGamTotp4->Py(),K1*BoostGamTotp4->Pz(),K1*BoostGamTotp4->E()); 
	//The total photon system is, in general, not massless!
	//BoostGamTotp4->Boost(*CangBoost);

	//Set the boosted,rescaled dimuon 4vector
	//BoostDiMup4->SetXYZM(BoostMup4->Px()+BoostAMup4->Px(),BoostMup4->Py()+BoostAMup4->Py(),BoostMup4->Pz()+BoostAMup4->Pz(),sqrt((BoostMup4->E() + BoostAMup4->E()) * (BoostMup4->E() + BoostAMup4->E()) - (BoostMup4->Px()+BoostAMup4->Px())*(BoostMup4->Px()+BoostAMup4->Px()) - (BoostMup4->Py()+BoostAMup4->Py())*(BoostMup4->Py()+BoostAMup4->Py()) - (BoostMup4->Pz()+BoostAMup4->Pz())*(BoostMup4->Pz()+BoostAMup4->Pz())));
	//cout << "TEST3" << endl;
	//
	//Now that we have the boosted values, let us also get the tracker smeared values for the muons
	//
	MuTrk.PartTrack = (TLorentzVector*)BoostMup4->Clone("BoostMu");
	SmearMup4 = MuTrk.SmearTrack();
	AMuTrk.PartTrack = (TLorentzVector*)BoostAMup4->Clone("BoostAMu");
	SmearAMup4 = AMuTrk.SmearTrack();
	//We have to catch negative energy states for photons because sometimes KKMC likes to 
	//include negative energy photons as a part of its "this is the beam energy spread of the initial system!!!"
	//which I obviously do not like
	if (BoostGamTotp4->E() > 0.0 && BoostGamTotp4->M() >= 0.0)
	  {
	    GamTrk.PartTrack = (TLorentzVector*)BoostGamTotp4->Clone("BoostGamTot");
	    SmearGamTotp4 = GamTrk.SmearECAL();
	  }
	SmearDiMup4->SetXYZM(SmearMup4->Px()+SmearAMup4->Px(),SmearMup4->Py()+SmearAMup4->Py(),SmearMup4->Pz()+SmearAMup4->Pz(),sqrt((SmearMup4->E() + SmearAMup4->E()) * (SmearMup4->E() + SmearAMup4->E()) - (SmearMup4->Px()+SmearAMup4->Px())*(SmearMup4->Px()+SmearAMup4->Px()) - (SmearMup4->Py()+SmearAMup4->Py())*(SmearMup4->Py()+SmearAMup4->Py()) - (SmearMup4->Pz()+SmearAMup4->Pz())*(SmearMup4->Pz()+SmearAMup4->Pz())));
	//cout << "TEST4" << endl;
	//The CoM from the momenta and energy values present ...
	truecm = gpsqs;
	//fill that tree and reseed things!
	t3->Fill();
	COUNT += 1;
	//tstamp->Set();
	//rand1->SetSeed(tstamp->GetNanoSec());
	//tstamp->Set();
	//rand2->SetSeed(tstamp->GetNanoSec());
	//tstamp->Set();
	//rand3->SetSeed(tstamp->GetNanoSec());
	//tstamp->Set();
	//rand4->SetSeed(tstamp->GetNanoSec());
      }
  }
  fout->Write();
  fout->Close();
}
