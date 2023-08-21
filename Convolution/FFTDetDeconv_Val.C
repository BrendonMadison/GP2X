#include <complex>

Double_t GetEntropy(TH1D *hEnt)
{
  Double_t entropy = 0.0;
  Double_t binnorm = 0.0;
  Double_t integral = hEnt->Integral();
  for (Int_t i = 1; i <= hEnt->GetNbinsX(); i++)
    {
      binnorm = hEnt->GetBinContent(i)/integral;
      entropy += -1.0 * binnorm * log(abs(binnorm));
    }
  return entropy;
}

Double_t ComputeSysErr(TH1D *hmodel, TH1D *hdata, Int_t NDOF, Double_t CHI2)
{

  Double_t s_sys = 0.0;
  Double_t cchi2 = 0.0;
  //cout << "---" << endl;
  while (abs(1.0*NDOF - cchi2) > 0.001)
    {
      s_sys += 0.01;
      cchi2 = 0.0;
      for (Int_t i = 1; i <= hmodel->GetNbinsX(); i++)
        {
	  cchi2 += pow(hmodel->GetBinContent(i) - hdata->GetBinContent(i),2.0)/(pow(hmodel->GetBinError(i),2.0) + pow(hdata->GetBinError(i),2.0) + TMath::Sign(1.0,CHI2-1.0*NDOF)*pow(s_sys,2.0));
        }
      //cout << " " << cchi2 << " " << NDOF << " " << s_sys << endl;
    }
  return s_sys*TMath::Sign(1.0,CHI2-1.0*NDOF);
}

//nbins should always be the length of xbin!
TH1D* GetEquipartition(TH1D *hin, Int_t nbins, Double_t XMAX, Double_t XMIN)
{
  //your input histogram should have lots of bins (comparable to total number of events)
  Double_t xbin[nbins];
  Double_t sum = 0.0;
  for (Int_t j = 1; j < hin->GetNbinsX(); j++) 
    {
      if (hin->GetBinCenter(j) > XMIN && hin->GetBinCenter(j) < XMAX) sum += hin->GetBinContent(j);
    }
  Double_t perbin = sum/nbins * 1.0;
  cout << hin->Integral() << endl;
  cout << nbins << endl;
  cout << perbin << endl;
  Double_t binsintegral = 0.0;
  Double_t checkint = 0.0;
  Int_t prevedge = 0;
  Int_t bincnt = 1;
  Int_t printcnt = 0;
  xbin[0] = XMIN;
  for (Int_t i = 1; i < hin->GetNbinsX()-1; i++)
    {
      if (hin->GetBinLowEdge(i) > XMIN && hin->GetBinLowEdge(i) < XMAX)
	{
	  binsintegral = hin->Integral(prevedge,i);
	  if (binsintegral > perbin)
	    {
	      cout << binsintegral << endl;
	      xbin[bincnt] = hin->GetBinLowEdge(i);
	      cout << xbin[bincnt] << endl;
	      prevedge = i;
	      bincnt++;
	    }
	}
    }

  //ensure the last bin is <<our>> max bin 
  xbin[bincnt] = XMAX;
  bincnt++;

  cout << bincnt << endl;

  TH1D *hnew = new TH1D("hnew","hnew",bincnt-1,xbin);

  for (Int_t i = 1; i < hin->GetNbinsX()-1; i++) hnew->Fill(hin->GetBinCenter(i),hin->GetBinContent(i));

  for (Int_t i = 1; i < hnew->GetNbinsX(); i++) cout << hnew->GetBinCenter(i) << " " << hnew->GetBinContent(i) << endl;

  return hnew;
}

TH1D* GetPeakHist(Double_t PEAK, Double_t SIGMA, TH1D *hDATA)
{
  //So you don't have to send as many arguments                                                                                    
  Int_t N = hDATA->GetNbinsX();
  Double_t XMIN = hDATA->GetBinLowEdge(1);
  Double_t XMAX = hDATA->GetBinLowEdge(N+1);
  
  //Create fit for the right side of the peak to discount the radiative effects
  TF1 *fgaus = new TF1("fgaus","gaus",PEAK,PEAK + 4.0*SIGMA);
  fgaus->SetParameters(1.0,PEAK,SIGMA);
  fgaus->SetParLimits(0,0.1,10000.0);
  fgaus->SetParLimits(1,PEAK*0.99,PEAK*1.01);
  fgaus->SetParLimits(2,SIGMA*0.99,SIGMA*1.01);
  fgaus->SetNpx(N);
  
  hDATA->Fit(fgaus,"L M R");
  TH1D *hgaus1 = (TH1D*)fgaus->GetHistogram();

  TF1 *fgaus2 = new TF1("fgaus2","gaus",XMIN,XMAX);
  fgaus2->SetParameters(fgaus->GetParameter(0),fgaus->GetParameter(1),fgaus->GetParameter(2));
  fgaus2->SetNpx(N);
  TH1D *hgaus2 = (TH1D*)fgaus2->GetHistogram();

  //delete fgaus2;
  //delete fgaus;
  //delete hgaus1;
  
  return hgaus2;
}

TH1D* GetMCFromDet(TH1D *hMC1, TH1D *hDET1, TH1D *hMC2, TH1D *hDET2)
{
  //So you don't have to send as many arguments
  Int_t N = hMC1->GetNbinsX();
  Double_t XMIN = hMC1->GetBinLowEdge(1); 
  Double_t XMAX = hMC1->GetBinLowEdge(N+1);
  
  //Temporary histograms for reading in complex, transformed, data
  TH1 *him1 = 0;
  TH1 *hre1 = 0;
  TH1 *him2 = 0;
  TH1 *hre2 = 0;
  TH1 *him3 = 0;
  TH1 *hre3 = 0;
  TH1 *him4 = 0;
  TH1 *hre4 = 0;
  TVirtualFFT::SetTransform(0);

  //for creating the convolution template from the first (simulation) file
  him1 = hMC1->FFT(him1,"IM");
  Double_t *im1 = new Double_t[N];
  for (Int_t i = 0; i < N; i++) im1[i] = him1->GetBinContent(i+1);

  hre1 = hMC1->FFT(hre1,"RE");
  Double_t *re1 = new Double_t[N];
  for (Int_t i = 0; i < N; i++) re1[i] = hre1->GetBinContent(i+1);  

  him2 = hDET1->FFT(him2,"IM");
  Double_t *im2 = new Double_t[N];
  for (Int_t i = 0; i < N; i++) im2[i] = him2->GetBinContent(i+1);  

  hre2 = hDET1->FFT(hre2,"RE");
  Double_t *re2 = new Double_t[N];
  for (Int_t i = 0; i < N; i++) re2[i] = hre2->GetBinContent(i+1);  

  //for subtracting out detector effects from the second (data) file
  him3 = hDET2->FFT(him3,"IM");
  Double_t *im3 = new Double_t[N];
  for (Int_t i = 0; i < N; i++) im3[i] = him3->GetBinContent(i+1);

  hre3 = hDET2->FFT(hre3,"RE");
  Double_t *re3 = new Double_t[N];
  for (Int_t i = 0; i < N; i++) re3[i] = hre3->GetBinContent(i+1);  

  him4 = hMC2->FFT(him4,"IM");
  Double_t *im4 = new Double_t[N];
  for (Int_t i = 0; i < N; i++) im4[i] = him4->GetBinContent(i+1);

  hre4 = hMC2->FFT(hre4,"RE");
  Double_t *re4 = new Double_t[N];
  for (Int_t i = 0; i < N; i++) re4[i] = hre4->GetBinContent(i+1);  

  //Compute complex difference
  //These first values are the "detector response"
  Double_t *reNew = new Double_t[N];
  Double_t *imNew = new Double_t[N];
  Double_t *reMC = new Double_t[N];
  Double_t *imMC = new Double_t[N];
  std::complex<Double_t> tmp1;
  std::complex<Double_t> tmp2;
  std::complex<Double_t> tmp3;
  std::complex<Double_t> tmp4;
  //If you wanted to get the detector response use DetRes
  //Otherwise the default is MCVal
  std::complex<Double_t> DetRes;
  std::complex<Double_t> MCVal;
  for (Int_t i = 0; i < N; i++)
    {
      tmp1 = re1[i] + 1i*im1[i]; //MC1
      tmp2 = re2[i] + 1i*im2[i]; //DET1
      tmp3 = re3[i] + 1i*im3[i]; //DET2
      tmp4 = re4[i] + 1i*im4[i]; //MC2

      DetRes = tmp4 / tmp3;
      MCVal = DetRes * tmp2;

      reNew[i] = DetRes.real();
      imNew[i] = DetRes.imag();

      reMC[i] = MCVal.real();
      imMC[i] = MCVal.imag();
    }
 
  //Now let's make a backward transform:
  //Use the following if you believe you are only real valued
  //TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  //Use the following if you believe you are complex valued.
  //Make sure to use magnitude later on!
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &N, "C2CBACKWARD");
  fft_back->SetPointsComplex(reMC,imMC);
  fft_back->Transform();
  TH1 *hb = 0;
  //Let's look at the output
  hb = TH1::TransformHisto(fft_back,hb,"Mag");
  TH1D *hbd = new TH1D("hbd","hbd",N,XMIN,XMAX);

  for (Int_t i = 0; i < N; i++) hbd->SetBinContent(i+1,hb->GetBinContent(i+1));

  return hbd;
  
} 

TH1D* GetGausFromDet(Double_t PEAK, Double_t SIGMA, TH1D *hMC1, TH1D *hDET1, TH1D *hMC2, TH1D *hDET2, char *ver, Int_t nbin, Double_t barr[])
{
  //So you don't have to send as many arguments                                                                                    
  Int_t N = hMC1->GetNbinsX();
  Double_t XMIN = hMC1->GetBinLowEdge(1);
  Double_t XMAX = hMC1->GetBinLowEdge(N+1);
  //First we get the MC deconvolved from detector using our simulation dataset
  TH1D *hConvMC = GetMCFromDet(hMC1,hDET1,hMC2,hDET2);

  TH1D *hIntMC1 = (TH1D*)hMC1->GetCumulative();
  TH1D *hIntDET1 = (TH1D*)hDET1->GetCumulative();
  TH1D *hIntMC2 = (TH1D*)hMC2->GetCumulative();
  TH1D *hIntDET2 = (TH1D*)hDET2->GetCumulative();

  TH1D *hIntConvMC = GetMCFromDet(hIntMC1,hIntDET1,hIntMC2,hIntDET2);

  //Now we need to get the new PEAK, SIGMA values since they changed from DET to MC level
  TF1 *fg = new TF1("fg","gaus",PEAK,PEAK + 2.0*SIGMA);
  fg->SetParameters(1.0,PEAK,SIGMA);
  fg->SetNpx();
  hConvMC->Fit(fg,"L M R");

  Double_t NEWPEAK = fg->GetParameter(1);
  Double_t NEWSIGMA = fg->GetParameter(2);

  //cout << NEWPEAK << endl;
  //cout << NEWSIGMA << endl;

  //Then we get the gaussian distribution corresponding to the deconvolved MC (using fits)
  //TH1D *hPeak = GetPeakHist(NEWPEAK,NEWSIGMA,hConvMC);

  TF1 *fg2 = new TF1("fg2","gaus",XMIN,XMAX);
  //fg2->SetParameters(fg->GetParameter(0),fg->GetParameter(1),fg->GetParameter(2));
  fg2->SetParameters(fg->GetParameter(0),PEAK,fg->GetParameter(2));
  fg2->SetNpx(N);
  TH1D *hGaus = (TH1D*)fg2->GetHistogram();

  //We also have the model/ideal gaussian
  //Honestly, how correct sigma is doesn't seem to matter that much as the outcome is more reliant on the data
  TF1 *fg3 = new TF1("fg3","gaus",XMIN,XMAX);
  fg3->SetParameters(fg->GetParameter(0),PEAK,SIGMA);
  fg3->SetNpx(N);
  TH1D *hIdeal = (TH1D*)fg3->GetHistogram();
  
  //Then we use this to deconvolve the radiative effects from the remaining MC level data
  //So the result should be approximately the gaussian BES distribution
  //Technically it will be a fourier series of it so you will likely see harmonic noise
  TH1D *hConvPeak = GetMCFromDet(hIdeal,hMC1,hGaus,hConvMC);

  hConvMC->SetName("hmc");
  //hConvMC->SaveAs(Form("HistNorebinMC_%s.root",ver));
  hConvMC = (TH1D*)hConvMC->Rebin(nbin,"hmc",barr);
  hConvMC->SaveAs(Form("HistMC_%s.root",ver));

  hIntConvMC->SetName(Form("HistIntMC_%s",ver));
  hIntConvMC = (TH1D*)hIntConvMC->Rebin(nbin,Form("HistIntMC_%s",ver),barr);
  //hIntConvMC->SaveAs(Form("HistIntMC_%s.root",ver));

  //Restore the Deconvolved histogram that used the CDF
  TH1D *hFromInt = (TH1D*)hIntConvMC->Clone();
  hFromInt->SetBinContent(1,hConvMC->GetBinContent(1));
  for(Int_t i = 2; i <= hFromInt->GetNbinsX(); i++) hFromInt->SetBinContent(i,hIntConvMC->GetBinContent(i) - hIntConvMC->GetBinContent(i-1));

  hFromInt->SetName(Form("HistFromInt_%s",ver));
  hFromInt = (TH1D*)hFromInt->Rebin(nbin,Form("HistFromInt_%s",ver),barr);
  //hFromInt->SaveAs(Form("HistFromInt_%s.root",ver));

  hIntMC2->SetName(Form("HistIntMC2_%s",ver));
  hIntMC2 = (TH1D*)hIntMC2->Rebin(nbin,Form("HistIntMC2_%s",ver),barr);
  //hIntMC2->SaveAs(Form("HistIntMC2_%s.root",ver));

  return hConvPeak;
}

TH1D* GetDetFromGaus(Double_t PEAK, Double_t SIGMA, TH1D *hMC1, TH1D *hDET1, TH1D *hMC2, TH1D *hDET2, TH1D *hGAUS)
{
  //So you don't have to send as many arguments                                                                                    
  Int_t N = hMC1->GetNbinsX();
  Double_t XMIN = hMC1->GetBinLowEdge(1);
  Double_t XMAX = hMC1->GetBinLowEdge(N+1);

  //We also have the model/ideal gaussian
  //Honestly, how correct sigma is doesn't seem to matter that much as the outcome is more reliant on the data
  TF1 *fg3 = new TF1("fg3","gaus",XMIN,XMAX);
  fg3->SetParameters(hGAUS->GetMaximum()/SIGMA,PEAK,SIGMA);
  fg3->SetNpx(N);
  TH1D *hIdeal = (TH1D*)fg3->GetHistogram();

  hIdeal->Scale(hGAUS->GetMaximum()/hIdeal->GetMaximum());

  //First we get the MC convolved from gaussian to MC level using our model and simulation dataset
  //TH1D *hConvMC = GetMCFromDet(hMC1,hIdeal,hMC2,hGAUS);

  //First we get the MC convolved from gaussian to MC level using our model and simulation dataset
  //TH1D *hConvDET = GetMCFromDet(hDET1,hMC1,hDET2,hConvMC);
  TH1D *hConvDET = GetMCFromDet(hDET1,hIdeal,hDET2,hGAUS);
  //TH1D *hConvDET = GetMCFromDet(hIdeal,hDET1,hGAUS,hDET2);

  //TCanvas *ctest = new TCanvas();
  //hIdeal->Draw("hist");
  //hGAUS->Draw("same");
  //ctest->SaveAs("HistCheck.pdf");

  return hConvDET;
}

void FFTDetDeconv_Val(Int_t dbins, Double_t loedge, Double_t hiedge, Int_t N, Double_t XMIN,Double_t XMAX, Double_t PEAK, Double_t SIGMA, char *fname1, char *fname2, char *version, char *MCSTR, char *DETSTR)
{

  //Int_t n = 10000;
  //Double_t XMIN = 250.0-10.0;
  //Double_t XMAX = 250.0+30.0;

  //We will adjust the binning for presentation purposes
  //Double_t loedge = PEAK * 0.968;
  //Double_t loedge = PEAK * (1.0 - 15.0 * SIGMA/PEAK);
  //Double_t hiedge = PEAK * 1.008;
  //Double_t hiedge = PEAK * (1.0 + 5.0 * SIGMA/PEAK);
  //Double_t loedge = 241.0;
  //Double_t hiedge = 252.0;
  //Double_t loedge = 242.0;
  //Double_t hiedge = 251.0;
  //number of bins for drawing
  //Int_t dbins = 180;
  //bins array for rebinning later
  Double_t binarr[dbins+1];
  for(Int_t i = 0; i <= dbins; i++) binarr[i] = loedge + i*(hiedge-loedge)/dbins;

  //TFile *f = new TFile("MuCutLarge.root");
  TFile *f = new TFile(fname1);
  TTree *t = (TTree*)f->Get("t");

  //TFile *f2 = new TFile("WHIZ_MuCut_7.root");
  TFile *f2 = new TFile(fname2);
  //TTree *t2 = (TTree*)f2->Get("t");
  TTree *t2 = (TTree*)f2->Get("t");

  TCanvas *cMain = new TCanvas();

  cout << "l284" << endl;

  //t->Draw(Form("DiMup4->E() + abs(DiMup4->P())>>hd1(%i,%f,%f)",N,XMIN,XMAX),Form("DiMup4->E() + abs(DiMup4->P()) > %f && DiMup4->E() + abs(DiMup4->P()) < %f",XMIN,XMAX),"");
  t->Draw(Form("%s>>hd1(%i,%f,%f)",MCSTR,N,XMIN,XMAX),Form("%s > %f && %s < %f",MCSTR,XMIN,MCSTR,XMAX),"");
  TH1D *hd1 = (TH1D*)gDirectory->Get("hd1");

  //t->Draw(Form("SDiMup4->E() + abs(SDiMup4->P())>>hs1(%i,%f,%f)",N,XMIN,XMAX),Form("SDiMup4->E() + abs(SDiMup4->P()) > %f && SDiMup4->E() + abs(SDiMup4->P()) < %f",XMIN,XMAX),"same");
  t->Draw(Form("%s>>hs1(%i,%f,%f)",DETSTR,N,XMIN,XMAX),Form("%s > %f && %s < %f",DETSTR,XMIN,DETSTR,XMAX),"same");
  TH1D *hs1 = (TH1D*)gDirectory->Get("hs1");

  //t2->Draw(Form("DiMup4->E() + abs(DiMup4->P())>>hd2(%i,%f,%f)",N,XMIN,XMAX),Form("DiMup4->E() + abs(DiMup4->P()) > %f && DiMup4->E() + abs(DiMup4->P()) < %f",XMIN,XMAX),"same");
  t2->Draw(Form("%s>>hd2(%i,%f,%f)",MCSTR,N,XMIN,XMAX),Form("%s > %f && %s < %f",MCSTR,XMIN,MCSTR,XMAX),"same");
  TH1D *hd2 = (TH1D*)gDirectory->Get("hd2");

  //t2->Draw(Form("SDiMup4->E() + abs(SDiMup4->P())>>hs2(%i,%f,%f)",N,XMIN,XMAX),Form("SDiMup4->E() + abs(SDiMup4->P()) > %f && SDiMup4->E() + abs(SDiMup4->P()) < %f",XMIN,XMAX),"same");
  t2->Draw(Form("%s>>hs2(%i,%f,%f)",DETSTR,N,XMIN,XMAX),Form("%s > %f && %s < %f",DETSTR,XMIN,DETSTR,XMAX),"same");
  TH1D *hs2 = (TH1D*)gDirectory->Get("hs2");

  //hd1->Smooth(256);
  //hd2->Smooth(256);
  //hs1->Smooth(256);
  //hs2->Smooth(256);

  //Save data plots for reference
  cMain->SaveAs(Form("InputPlots_%s.pdf",version));

  //Get hplot for the purpose of plotting
  //Get hnorm for the purpose of normalization
  t2->Draw(Form("%s>>hplot(%i,%f,%f)",DETSTR,200,loedge,hiedge),Form("%s > %f && %s < %f",DETSTR,loedge,DETSTR,hiedge),"goff");
  TH1D *hplot = (TH1D*)gDirectory->Get("hplot");
  hplot->SetLineColor(2);

  //t2->Draw(Form("SDiMup4->E() + abs(SDiMup4->P())>>hnorm(%i,%f,%f)",200,0.0,PEAK*2.0),Form("SDiMup4->E() + abs(SDiMup4->P()) > %f && SDiMup4->E() + abs(SDiMup4->P()) < %f",0.0,PEAK*2.0),"goff");
  t2->Draw(Form("%s>>hnorm(%i,%f,%f)",DETSTR,200,0.0,PEAK*2.0),Form("%s > %f && %s < %f",DETSTR,0.0,DETSTR,PEAK*2.0),"goff");
  TH1D *hnorm = (TH1D*)gDirectory->Get("hnorm");

  //Get the deconvolved gaussian distribution
  TH1D *h1 = GetGausFromDet(PEAK, SIGMA, hd1, hs1, hd2 , hs2,version,dbins,binarr);
  //Then we will rebin it for presentation purposes
  TH1D *h5 = new TH1D("h5","h5",200,loedge,hiedge);
  for (Int_t i = 0; i < N; i++) h5->Fill(h1->GetBinCenter(i+1),h1->GetBinContent(i+1));

  cout << "l284" << endl;

  //Rescale to have the correct integral
  //Use h1 as h5 has a smaller data window and doesn't include all events
  h5->Scale(hnorm->Integral() /h1->Integral());

  h5->SetLineColor(3);
  h5->SetTitle("Deconvolved Gaussian Compared Detector Level");
  h5->Draw("hist");
  hplot->Draw("same");

  cMain->SaveAs(Form("GaussDetDeconv_%s.pdf",version));

  //Now we will 'go backwards' by fitting this to a gaussian and then re-convolving to get the
  //detector level result that has the fitted parameters.
  //This will allow us to get a better idea on the residual and chi2

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPalette(109);
  gStyle->SetGridStyle(7);
  gStyle->SetGridColor(17);
  //gStyle->SetFitFormat("5.8g");
  TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.33);
  pad1->SetBottomMargin(0.2);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.2);
  pad2->SetBottomMargin(0.2);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  h5->Draw("hist");
  //h5->GetXaxis()->SetRangeUser(245.0,255.0);
  h5->GetXaxis()->SetTitle("#sqrt{s} Deconvolved from #sqrt{s_{p}} [GeV]");
  h5->GetYaxis()->SetTitle("Bin Counts for ~100 fb^{-1} Dataset");
  h5->SetTitle("#sqrt{s} Deconvolution from Detector Data Test");

  //Now we need to get the new PEAK, SIGMA values since they changed from DET to MC level
  TF1 *fFit = new TF1("fFit","gaus",PEAK - 2.0 * SIGMA ,PEAK + 3.0*SIGMA);
  fFit->SetParameters(1.0,PEAK,SIGMA);
  fFit->SetNpx(N);

  //Now we make a copy that has the full window
  TF1 *fFitFull = new TF1("fFitFull","gaus",XMIN ,XMAX);
  fFitFull->SetParameters(fFit->GetParameter(0),fFit->GetParameter(1),fFit->GetParameter(2));
  fFitFull->SetNpx(N);
  TH1D *hFit = (TH1D*)fFitFull->GetHistogram();

  //TH1D *hDetFit = GetDetFromGaus(PEAK, SIGMA, hd1, hs1, hd2 , hs2, hFit);  

  TH1D *hDetFit = GetMCFromDet(hs1,hFit,hs2,h1);

  hDetFit->Scale(hs2->Integral() / hDetFit->Integral());
  //hDetFit->SaveAs(Form("MCfromDet_%s.root",version));
  //hDetFit->Scale(hs2->Integral(hs2->GetXaxis()->FindBin(loedge),hs2->GetXaxis()->FindBin(hiedge)) / hDetFit->Integral(hDetFit->GetXaxis()->FindBin(loedge),hDetFit->GetXaxis()->FindBin(hiedge)));
  //hDetFit->Rebin(10);
  //hs2->Rebin(10);
  hDetFit->SetLineColorAlpha(30,0.5);
  hs2->SetLineColorAlpha(46,0.5);

  //h5->Scale(hs2->Integral()/h5->Integral());
  //fFit->SetParameter(0,h5->GetMaximum());
  h5 = (TH1D*)h5->Rebin(50,"h5",binarr);
  h5->Fit(fFit,"L M R");
  h5->Draw("hist");

  fFit->SetRange(loedge,hiedge);
  fFit->Draw("same");

  hDetFit->Scale(h5->Integral()/hDetFit->Integral());
  hs2->Scale(h5->Integral()/hs2->Integral());

  hDetFit = (TH1D*)hDetFit->Rebin(50,"hDetFit",binarr);
  hs2 = (TH1D*)hs2->Rebin(50,"hs2",binarr);

  hDetFit->Draw("same");
  hs2->Draw("same");
  hDetFit->GetXaxis()->SetRangeUser(loedge,hiedge);
  hDetFit->GetXaxis()->SetTitle("#sqrt{s} Deconvolved from #sqrt{s_{p}} [GeV]");
  hDetFit->GetYaxis()->SetTitle("Bin Counts for ~100 fb^{-1} Dataset");
  hDetFit->SetTitle("#sqrt{s} Deconvolution from Detector Data Test");

  //Get entropy of the fit
  Double_t entFit = GetEntropy(hDetFit);
  cout << "Entropy of the deconvolved (Gaus) fit: " << entFit << endl;

  //We need to get the number of bins in our range...
  Int_t nbins = 0;
  for (Int_t i=1; i < hDetFit->GetNbinsX(); i++)
    {
      if (hDetFit->GetBinCenter(i) >= loedge && hDetFit->GetBinCenter(i) <= hiedge)
	{
	  nbins += 1;
	}
    }

  //hRes->GetYaxis()->SetLabelFont(63);
  //hRes->GetYaxis()->SetLabelSize(16);
  Double_t pull = 0.0;
  Double_t unc = 0.0;
  Double_t Chi2 = 0.0;
  Int_t cnt = 1;

  //cout << "Successful?" << endl;
  //heq->Scale(1.0/heq->GetBinContent(10));
  //TH1D *hseq = (TH1D*)heq->Clone();
  //for (Int_t i=0; i <= hseq->GetNbinsX(); i++) hseq->SetBinContent(i,0.0);
  //for (Int_t i=1; i < hs2->GetNbinsX(); i++) hseq->Fill(hs2->GetBinCenter(i),hs2->GetBinContent(i));
  //hseq->Scale(heq->Integral()/hseq->Integral());

  // for (Int_t i=1;i<=heq->GetNbinsX();i++) 
  //   {
  //     if (heq->GetBinCenter(i) >= loedge && heq->GetBinCenter(i) <= hiedge)
  // 	{
  // 	  unc = sqrt(pow(heq->GetBinError(i),2.0) + pow(hseq->GetBinError(i),2.0));
  // 	  pull = (heq->GetBinContent(i) - hseq->GetBinContent(i));
  // 	  Chi2 += (pull*pull)/(unc*unc);
  // 	  hRes->SetBinContent(cnt,pull/unc);
  // 	  cnt += 1;
  // 	}
  //   }

  //fix the bin errors as rescaling and transforms have made these super screwy

  //we need to normalize but only over the plot range!

  Double_t sumdf = 0.0;
  Double_t sums2 = 0.0;

  for (Int_t i = 1; i <=hDetFit->GetNbinsX();i++)
    {
      if (hDetFit->GetBinCenter(i) >= loedge && hDetFit->GetBinCenter(i) <= hiedge)
	{
	  sumdf += hDetFit->GetBinContent(i);
	  sums2 += hs2->GetBinContent(i);
	}
    }

  cout << "Sum of Detector Fit: " << sumdf << endl;
  cout << "Sum of data (hs2): " << sums2 << endl;
  hDetFit->Scale(sums2/sumdf);
  for (Int_t i=1; i<=hDetFit->GetNbinsX();i++) hDetFit->SetBinError(i,sqrt(hDetFit->GetBinContent(i)));
  for (Int_t i=1; i<=hs2->GetNbinsX();i++) hs2->SetBinError(i,sqrt(hs2->GetBinContent(i)));

  for (Int_t i=2;i<=hDetFit->GetNbinsX();i++) 
    {
      if (hDetFit->GetBinCenter(i) >= loedge && hDetFit->GetBinCenter(i) <= hiedge)
  	{
  	  unc = sqrt(pow(hDetFit->GetBinError(i),2.0) + pow(hs2->GetBinError(i),2.0));
  	  pull = (hDetFit->GetBinContent(i) - hs2->GetBinContent(i));
  	  Chi2 += (pull*pull)/(unc*unc);
  	  //hRes->SetBinContent(cnt,pull/unc);
  	  cnt += 1;
	  cout << i << " " << hDetFit->GetBinCenter(i) << " " << hs2->GetBinCenter(i) << endl;
  	}
    }

  //save detector level fit for averaging later
  hDetFit->SetName(Form("HistDetFit_%s",version));
  hDetFit->SaveAs(Form("HistDetFit_%s.root",version));
  hs2->SaveAs(Form("HistData_%s.root",version));

  hd2 = (TH1D*)hd2->Rebin(50,"hd2",binarr);
  hd2->SaveAs(Form("HistMCData_%s.root",version));

  //TH1D *heq = (TH1D*)GetEquipartition(hDetFit, 40, hiedge, loedge);
  //TH1D *hseq = (TH1D*)heq->Clone();
  //for (Int_t i=0; i <= hseq->GetNbinsX(); i++) hseq->SetBinContent(i,0.0);
  //for (Int_t i=1; i < hs2->GetNbinsX(); i++) hseq->Fill(hs2->GetBinCenter(i),hs2->GetBinContent(i));
  //heq->Scale(1.0/heq->GetBinContent(10));
  //cout << "check" << endl;
  //for (Int_t i=1; i < heq->GetNbinsX(); i++) cout << heq->GetBinCenter(i) << " " << heq->GetBinContent(i) << endl;
  //hseq->Scale(heq->GetBinContent(5)/hseq->GetBinContent(5));
  //for (Int_t i=1; i < hseq->GetNbinsX(); i++) cout << hseq->GetBinCenter(i) << " " << hseq->GetBinContent(i) << endl;
  //heq->Draw("hist e2 PLC");
  //hseq->Draw("hist e2 same PLC");

  Double_t equichi2 = 0.0;
  TH1D *hRes = new TH1D("hRes","hRes",40,loedge,hiedge);
  //hRes->GetXaxis()->SetLabelFont(63);
  //hRes->GetXaxis()->SetLabelSize(16);
  // hRes->GetXaxis()->SetTitle("Pulls");
  // cnt = 1;

  // for (Int_t i=1;i<heq->GetNbinsX();i++) 
  //   {
  // 	  unc = sqrt(pow(heq->GetBinError(i),2.0) + pow(hseq->GetBinError(i),2.0));
  // 	  pull = (heq->GetBinContent(i) - hseq->GetBinContent(i));
  // 	  equichi2 += (pull*pull)/(unc*unc);
  // 	  hRes->SetBinContent(cnt,pull);
  // 	  cnt += 1;
  // 	  cout << i << " " << heq->GetBinCenter(i) << " " << hseq->GetBinCenter(i) << endl;
  //   }

  // cout << "Equipartition uncertainty: " << equichi2 << endl;

  // pad2->cd();

  // hRes->GetXaxis()->SetRangeUser(loedge,hiedge);
  // hRes->GetXaxis()->SetTitle("");
  // hRes->GetYaxis()->SetTitle("");

  // hRes->SetTitle("Pulls Distribution");
  // //hRes->SetTitleSize(0.4);
  // //hRes->GetXaxis()->SetLabelSize(0.4);
  // hRes->GetYaxis()->SetTitle("Pulls");
  // hRes->GetXaxis()->SetTitle("#sqrt{s_{p}} [GeV]");
  // //hRes->GetYaxis()->SetLabelSize(0.4);
  // hRes->Draw("same");

  pad1->cd();

  //Double_t syserr = ComputeSysErr(hDetFit, hs2, nbins, Chi2);
  Double_t syserr = (Chi2/(nbins-2.0) - 1.0) * fFit->GetParError(1);

  TLegend *legend = new TLegend(0.1,0.48,0.48,0.9);
 legend->SetHeader(Form("%s",version),"C"); // option "C" allows to center the header
  legend->AddEntry(h5,"Deconvolved Gaus","l");
  legend->AddEntry(fFit,"Fit To Gaus","l");
  legend->AddEntry(hs2,"Det. Data","l");
  legend->AddEntry(hDetFit,"Fit Convolved to Det.","l");
  legend->AddEntry((TObject*)0,Form("#mu = %f #pm %f",fFit->GetParameter(1),fFit->GetParError(1)),"");
  legend->AddEntry((TObject*)0,Form("#sigma = %f #pm %f",fFit->GetParameter(2),fFit->GetParError(2)),"");
  legend->AddEntry((TObject*)0,Form("#chi^{2} / NDoF = %f/%i",Chi2,nbins-2),"");
  legend->AddEntry((TObject*)0,Form("Systematic Unc. = %f",syserr),"");
  legend->Draw();

  cout << "Chi2/NDoF Statistic: " << Chi2 << "/" << nbins << endl;

  cout << "Systematic Uncertainty: " << syserr  << endl;

  cMain->SaveAs(Form("FitPlot_%s.pdf",version));

  //t2->Draw(Form("SDiMup4->E() + abs(SDiMup4->P())>>hsPlot(%i,%f,%f)",200,loedge,hiedge),Form("SDiMup4->E() + abs(SDiMup4->P()) > %f && SDiMup4->E() + abs(SDiMup4->P()) < %f",loedge,hiedge),"same");
  //TH1D *hsPlot = (TH1D*)gDirectory->Get("hsPlot");
  //TH1D *h6 = new TH1D("h6","h6",200,loedge,hiedge);
  //for (Int_t i = 0; i < N; i++) h6->Fill(h6->GetBinCenter(i+1),hDetFit->GetBinContent(i+1));

  //h6->Scale(hsPlot->Integral() / h6->Integral());
  //h6->SetLineColorAlpha(30,0.5);
  //hsPlot->SetLineColorAlpha(46,0.5);

  //h6->Draw("hist");
  //hsPlot->Draw("same");
  //Notes
  //1.) The amount of data used in the "deconvolution template" is directly correlated to how well it performs. With small data set you get lots of high frequency noise and harmonic noise.
  //
  //2.) With that said, it is indeed feasible to use this to deconvolve the detector effects out of the detector level data to get data that is similar to the monte carlo data WITHOUT USING THE MONTE CARLO DATA. Meaning its feasible to do with the actual experiment.
  //
  //3.) We should be able to do a similar thing to deconvolve the sqrt(s_p) into a gaussian of stdev determined by the original E_+ and E_- values. This then gets around the radiative effects and the "sqrt(s_p) isn't perfect estimator" effects.
  //
  //4.) So next steps? 
  //a.) Clean up this code
  //b.) Make code parametric so it can run off arguments
  //c.) Make new function to do deconvolve radiative effects etc. after removing detector effects
  //d.) Fit the final deconvolution to a gaussian and get the estimated precision
  //e.) Possibly do a convolution (so invert the convolution) of the fitted gaussian so you get the fitted detector level results but you'd only have 2 parameters, stdev and mean.
  //f.) Redo normalization to not occur until the end (when its gaussian or simpler distribution)

}
