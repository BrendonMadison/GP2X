//Function for filling FtVal with values from kname (KKMC or BHWIDE file)

//IDEA: STICH DECONVOLUTION REGIONS TOGETHER TO MINIMIZE PEAK REGION HARMONIC NOISE

using namespace RooFit;

void ValSPMC(char* infile, Double_t MEAN, Double_t VAR)
{
  cout << "Running ValSP with  \n" << infile << endl;
  cout << MEAN << endl;
  cout << VAR << endl;
  //RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("maxSteps",100);
  //RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("epsRel",0.0000001);
  //RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("epsAbs",0.0000001);

  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","FUMILI");
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSL","Levenberg-Marquardt");
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.0000000001);
  //ROOT::Math::MinimizerOptions::SetDefaultPrecision(0.0000000001);
  //ROOT::Math::MinimizerOptions().SetMaxIterations(10000000);
  //ROOT::Math::MinimizerOptions().SetMaxFunctionCalls(10000000);
  ROOT::Math::MinimizerOptions().Print();
  
  //Plot the dimuon and electron boost (beta) values according to P/E
  //Create canvas and turn the stat box OFF
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  gStyle->SetOptStat(0);

  string xtitle="Center-of-Mass Energy, #sqrt{s_{p}} [GeV]";
  //string xtitle="Center";
  //Double_t ymaxh = 40000.0;
  Double_t yminh = 0.0;
  Double_t ymaxp = 5.5;
  //Double_t XMIN = MEAN - 10.0*VAR;
  //Double_t XMAX = MEAN + 3.0*VAR;
  Double_t XMIN = 242.0;
  Double_t XMAX = 251.0;
  Double_t bwidth = (XMAX-XMIN)/100.0;
  Double_t Escale_ini = MEAN;
  Double_t sigma_ini = 100.0*VAR/MEAN;
  Double_t msig = 0.0022 /2.0;
  Double_t psig = 0.001875 /2.0;
  Double_t dEhalf = 0.5*bwidth;
  Size_t msize = 1.5;

  gStyle->SetTitleFont(62, ""); // set the pad title font
  gStyle->SetMarkerSize(msize);
  gStyle->SetLegendTextSize(0.045);

  RooRealVar x("x","x",XMIN,XMAX,"GeV");             // Fit the beam energy in specified range (Last argument is the unit)

  
  TFile *f = new TFile(infile);
  //If you are using the normal data settings
  TTree *t = (TTree*)f->Get("t"); //If using GP2X use this
  //TTree *t = (TTree*)f->Get("DiMuonTree"); //If using ILCOSFT use this
  //Int_t k = t->Draw(Form("250.0 - sqrt((250.0)**2 - 2.0*250.0*DiMup4->E() + 91.1876**2)>>hBetaSP(100,%f,%f)",XMIN,XMAX),Form("250.0 - sqrt((250.0)**2 - 2.0*250.0*DiMup4->E() + 91.1876**2) > %f && 250.0 - sqrt((250.0)**2 - 2.0*250.0*DiMup4->E() + 91.1876**2) < %f",XMIN,XMAX),"");
  Int_t k = t->Draw(Form("DiMup4->E() + abs(DiMup4->P())>>hBetaSP(90,%f,%f)",XMIN,XMAX),Form("DiMup4->E() + abs(DiMup4->P()) > %f && DiMup4->E() + abs(DiMup4->P()) < %f",XMIN,XMAX),""); //USE THIS IF GP2X
  //Int_t k = t->Draw(Form("mcpSqrtsPlus>>hBetaSP(100,%f,%f)",XMIN,XMAX),Form("abs(mcpMuonCosTheta[0]) < 0.99619 && abs(mcpMuonCosTheta[1]) < 0.99619 && mcpSqrtsPlus > %f && mcpSqrtsPlus < %f",XMIN,XMAX),""); //USE THIS IF ILCSOFT
  TH1D *h1 = (TH1D*)gDirectory->Get("hBetaSP");
  //TH1D *h1 = (TH1D*)gDirectory->Get("AveMC");
  //TH1D *h1 = (TH1D*)gDirectory->Get("hs1");
  //TH1D *hd2 = (TH1D*)gDirectory->Get("hs1");

  //h1->Scale(hd2->Integral()/h1->Integral());
  //for(Int_t i = 1; i < h1->GetNbinsX(); i++) h1->SetBinError(i,sqrt(h1->GetBinContent(i)));

  Double_t ymaxh = 1.5 * h1->GetMaximum();

  h1->Print();
  RooDataHist db("db", "db", x, Import(*h1));
  x.setRange("FitRange",XMIN,XMAX);

  RooRealVar Escale("E_{0}","Escale",Escale_ini,0.999*Escale_ini,1.001*Escale_ini,"GeV");           // Energy scale of full energy
  //RooRealVar Escale1("E_{1}","Escale1",Escale_ini,0.99*Escale_ini,1.01*Escale_ini,"GeV");           // Energy scale of full energy
  //RooRealVar Escale2("E2","Escale2",Escale_ini*0.99,0.98*Escale_ini,1.01*Escale_ini,"GeV");
  //RooRealVar Escale("E_{0}","Escale",Escale_ini,1.0*Escale_ini,1.00*Escale_ini,"GeV");           // Energy scale of full energy
  //RooRealVar Escale("E_{0}","Escale",Escale_ini,1.0*Escale_ini,1.0*Escale_ini,"GeV");           // Energy scale of full energy
  RooRealVar sigma("#sigma","sigma",sigma_ini,0.8*sigma_ini,1.25*sigma_ini,"%");                                       // Fractional resolution in per cent of Gaussian used in convolution
  RooRealVar s_peak("s_{peak}","s_{peak}",sigma_ini,0.95*sigma_ini,1.05*sigma_ini,"%");                                       // Fractional resolution in per cent of Gaussian used in convolution
  //RooRealVar m_peak("m_{peak}","m_{peak}",msig,0.99*msig,1.01*msig,"%");                                       // Fractional resolution in per cent of Gaussian used in convolution
  //RooRealVar p_peak("p_{peak}","p_{peak}",psig,0.99*psig,1.01*psig,"%");                                       // Fractional resolution in per cent of Gaussian used in convolution
  //RooRealVar sigma("#sigma","sigma",0.164,0.164,0.164,"%");                                       // Fractional resolution in per cent of Gaussian used in convolution
  //RooRealVar sigma("#sigma","sigma",sigma_ini,sigma_ini,sigma_ini,"%");
  // RooRealVar alpha("#alpha","alpha",15.0,5.0,25.0);
  //RooRealVar alpha("#alpha","alpha",24.3,0.1,28.0);// Alpha parameter of Beta distribution
  RooRealVar alpha("#alpha","alpha",1.0,0.5,35.0);
  //RooRealVar alpha("#alpha","alpha",28.0,28.0,28.0);// Alpha parameter of Beta distribution
  
  //RooRealVar alpha2("#alpha2","alpha2",15.0,5.0,25.0);                                                // Alpha parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",12.1,12.1,12.1);                                                // Alpha parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",11.66,11.66,11.66);                                                // Alpha parameter of Beta distribution
  //RooRealVar beta("#beta","beta",15.0,10.0,45.0);                                                // Alpha parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",0.20,0.15,0.4);                                                   // Beta parameter of Beta distribution
  //RooRealVar beta("#beta","beta",0.2695,0.2695,0.2695);                                                   // Beta parameter of Beta distribution
  RooRealVar beta("#beta","beta",0.2695,0.15,0.7);                                                   // Beta parameter of Beta distribution
  //RooRealVar beta("#beta","beta",0.33334,0.33334,0.33334);                                                   // Beta parameter of Beta distribution
  //RooRealVar beta2("#beta2","beta2",0.20,0.15,0.4);                                                   // Beta parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",0.11710645899206922,0.11710645899206922,0.11710645899206922);                                                // Alpha parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",0.1130971624410773335262428677792249,0.11309716244792249,0.11309716244792249);                                                // Alpha parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",20,5,50);                                                // Alpha parameter of Beta distribution
  //RooRealVar beta("#beta","beta",0.2,0.2,0.10773335262428677);                                                   // Beta parameter of Beta distribution
  RooRealVar fpeak("f_{peak}","fpeak",0.11,0.01,0.35);                                               // Fraction of the Gaussian peak in the model
  //RooRealVar fpeak("f_{peak}","fpeak",0.2573,0.2573,0.2573);                                               // Fraction of the Gaussian peak in the model
  //RooRealVar fpeak2("f_{peak2}","fpeak2",0.35,0.1,0.95);                                               // Fraction of the Gaussian peak in the model

  //RooRealVar vXMAX("vXMAX","vXMAX",XMAX,XMAX,XMAX);

  //RooGenericPdf model("model","(vXMAX + 0.001 - x)^(alpha - 1.0) * (vXMAX + 1.001 - x)^(-1.0*alpha - beta)",RooArgList(x,vXMAX,alpha,beta));

  RooMyConvolvedBetaPdf tail("Convolved Beta","compiled class",x,Escale,alpha,beta,sigma);
  //RooMyConvolvedBetaPdf tail2("Convolved Beta","compiled class",x,Escale2,alpha2,beta2,sigma);
  RooMyNonStandardGaussianPdf peak("peak","peak pdf",x,Escale,s_peak);    // Note sigma here is the fractional resolution in percent
  //RooMyNonStandardGaussianPdf mpeak("mpeak","peak pdf",x,Escale,m_peak);    // Note sigma here is the fractional resolution in percent
  //RooMyNonStandardGaussianPdf ppeak("ppeak","peak pdf",x,Escale,p_peak);    // Note sigma here is the fractional resolution in percent
  //RooProdPdf peak("peak","peak pdf",mpeak,ppeak,0.0);
  RooAddPdf model("Model","model",RooArgList(peak,tail),fpeak);          // fpeak is the peak fraction
  //RooAddPdf model("Model","model",RooArgList(pmodel,tail2),fpeak2);          // fpeak is the peak fraction

  RooArgSet observables(model);

  // Plot data
  RooPlot* frame=x.frame(Title("Generator Level CoM Estimate #sqrt{s_{p}}"));
  db.plotOn(frame);

  RooChi2Var chi2("chi2","chi2",model,db,Range("FitRange"),
		  DataError(RooAbsData::Expected));   // RooChi2Var needs the fit range specified in calculating its statistic
  // Otherwise it may default to input binned dataset range that includes all of the TH1 bins
  // See https://sft.its.cern.ch/jira/browse/ROOT-10038
  double chi2f{};
  double binw{};
  double rchisq{};
  double nevts{};
  int nparam{};
  double bins{};
  int ndof{};

  cout << "Initial Chi-squared " << chi2.getVal() << endl;
  //     RooFitResult *result = model.chi2FitTo(db,Range("FitRange"),Hesse(false),Minos(false),Save(),IntegrateBins(1.0e-12));

  //model.chi2FitTo(db,Range("FitRange"));

  model.fitTo(db,Range("FitRange"));

  cout << "Fitted Chi-squared " << setprecision(10) << chi2.getVal() << endl;
  chi2f = chi2.getVal();
  cout << "Escale  " << setprecision(10) << Escale.getValV()  << " +- " << Escale.getError() << " Fractional error " << Escale.getError()/Escale.getValV() << endl;
  //cout << "Escale  " << setprecision(10) << Escale1.getValV()  << " +- " << Escale1.getError() << " Fractional error " << Escale1.getError()/Escale1.getValV() << endl;
  cout << "sigma " << setprecision(10) << sigma.getValV() << " +- " << sigma.getError() << endl;
  cout << "alpha " << setprecision(10) << alpha.getValV() << " +- " << alpha.getError() << endl;
  cout << "beta " << setprecision(10) << beta.getValV() << " +- " << beta.getError() << endl;
  cout << "fpeak " << setprecision(10) << fpeak.getValV() << " +- " << fpeak.getError() << endl;
  cout << "s_peak " << setprecision(10) << s_peak.getValV() << " +- " << s_peak.getError() << endl;
  //cout << "m_sig " << setprecision(10) << m_peak.getValV() << " +- " << m_peak.getError() << endl;
  //cout << "p_sig " << setprecision(10) << p_peak.getValV() << " +- " << p_peak.getError() << endl;

  //sigma.setConstant(kTRUE);
  alpha.setConstant(kTRUE);
  beta.setConstant(kTRUE);
  fpeak.setConstant(kTRUE);
  //Escale1.setConstant(kTRUE);
  //s_peak.setConstant(kTRUE);

  model.fitTo(db,Range("FitRange"));

  cout << "Fitted Chi-squared " << setprecision(10) << chi2.getVal() << endl;
  chi2f = chi2.getVal();
  cout << "Escale  " << setprecision(10) << Escale.getValV()  << " +- " << Escale.getError() << " Fractional error " << Escale.getError()/Escale.getValV() << endl;
  //cout << "Escale1  " << setprecision(10) << Escale1.getValV()  << " +- " << Escale1.getError() << " Fractional error " << Escale1.getError()/Escale1.getValV() << endl;
  cout << "sigma " << setprecision(10) << sigma.getValV() << " +- " << sigma.getError() << endl;
  cout << "alpha " << setprecision(10) << alpha.getValV() << " +- " << alpha.getError() << endl;
  cout << "beta " << setprecision(10) << beta.getValV() << " +- " << beta.getError() << endl;
  cout << "fpeak " << setprecision(10) << fpeak.getValV() << " +- " << fpeak.getError() << endl;
  //cout << "m_sig " << setprecision(10) << m_peak.getValV() << " +- " << m_peak.getError() << endl;
  //cout << "p_sig " << setprecision(10) << p_peak.getValV() << " +- " << p_peak.getError() << endl;
  cout << "s_peak " << setprecision(10) << s_peak.getValV() << " +- " << s_peak.getError() << endl;

  sigma.setConstant(kTRUE);
  //alpha.setConstant(kTRUE);
  //beta.setConstant(kTRUE);
  //fpeak.setConstant(kTRUE);
  //Escale1.setConstant(kTRUE);
  s_peak.setConstant(kTRUE);

  model.fitTo(db,Range("FitRange"));

  cout << "Fitted Chi-squared " << setprecision(10) << chi2.getVal() << endl;
  chi2f = chi2.getVal();
  cout << "Escale  " << setprecision(10) << Escale.getValV()  << " +- " << Escale.getError() << " Fractional error " << Escale.getError()/Escale.getValV() << endl;
  //cout << "Escale1  " << setprecision(10) << Escale1.getValV()  << " +- " << Escale1.getError() << " Fractional error " << Escale1.getError()/Escale1.getValV() << endl;
  cout << "sigma " << setprecision(10) << sigma.getValV() << " +- " << sigma.getError() << endl;
  cout << "alpha " << setprecision(10) << alpha.getValV() << " +- " << alpha.getError() << endl;
  cout << "beta " << setprecision(10) << beta.getValV() << " +- " << beta.getError() << endl;
  cout << "fpeak " << setprecision(10) << fpeak.getValV() << " +- " << fpeak.getError() << endl;
  //cout << "m_sig " << setprecision(10) << m_peak.getValV() << " +- " << m_peak.getError() << endl;
  //cout << "p_sig " << setprecision(10) << p_peak.getValV() << " +- " << p_peak.getError() << endl;
  cout << "s_peak " << setprecision(10) << s_peak.getValV() << " +- " << s_peak.getError() << endl;

  //compute the chi2 the correct way...
  TH1D *hmodel = (TH1D*)model.createHistogram("hmodel",x,RooFit::Binning(25,XMIN,XMAX));
  TH1D *hdb = (TH1D*)db.createHistogram("hdb",x,RooFit::Binning(25,XMIN,XMAX));
  //hmodel->Scale(500000.0/hmodel->Integral());
  //hdb->Scale(500000.0/hdb->Integral());
  //for (Int_t i = 1; i <= hmodel->GetNbinsX(); i++) hmodel->SetBinError(i,sqrt(hmodel->GetBinContent(i)));
  //for (Int_t i = 1; i <= hdb->GetNbinsX(); i++) hdb->SetBinError(i,sqrt(hdb->GetBinContent(i)) + 0.01*hdb->GetBinContent(i));
  //hmodel->Scale(h1->GetBinContent(10)/hmodel->GetBinContent(5));
  cout << "model hist check " << hmodel->GetBinContent(5) << endl;
  cout << "data hist check " << h1->GetBinContent(5) << endl;
  cout << "data err check " << h1->GetBinError(5) << endl;
  Double_t hchi2 = 0.0;
  TH1D *hres2 = (TH1D*)hdb->Clone();
  for (Int_t i = 0; i <= hmodel->GetNbinsX(); i++) 
    {
      hres2->SetBinContent(i,0);
      if (hmodel->GetBinCenter(i) > XMIN && hmodel->GetBinCenter(i) < XMAX) 
	{
	  hchi2 += pow(hmodel->GetBinContent(i) - hdb->GetBinContent(i),2.0)/(pow(hdb->GetBinError(i),2.0) + hmodel->GetBinContent(i));
	  hres2->SetBinContent(i,(hmodel->GetBinContent(i) - hdb->GetBinContent(i))/sqrt(pow(hdb->GetBinError(i),2.0) + hmodel->GetBinContent(i)));
	}
    }

  cout << "Histogram computed chi2 " << hchi2 << endl;

  model.plotOn(frame,Range("FitRange"),NormRange("FitRange"),Components(tail),LineStyle(kSolid),LineColor(kMagenta));
  model.plotOn(frame,Range("FitRange"),NormRange("FitRange"),Components(peak),LineStyle(kSolid),LineColor(kRed));
  model.plotOn(frame,Range("FitRange"),NormRange("FitRange"));
  model.paramOn(frame,Label(""),Parameters(RooArgSet(Escale,sigma,alpha,beta,fpeak,s_peak)),Layout(0.17,0.485,0.89),Format("NEU",AutoPrecision(2)));

  RooArgSet *flparams = model.getParameters(observables);
  nparam = (flparams->selectByAttrib("Constant",kFALSE))->getSize();
  cout << "nparam = " << nparam << endl;
  rchisq = frame->chiSquare(nparam);
  cout << "Reduced chi-squared?? " << rchisq << endl;
  cout << "Inferred dof " << chi2.getVal()/rchisq << endl;
  nevts=frame->getFitRangeNEvt();
  cout << "nevts = " << nevts << endl;
  binw=frame->getFitRangeBinW();
  cout << "binw = " << binw << endl;
  bins = (XMAX-XMIN)/binw;
  cout << "Number of bins in fit " << bins << endl;
  ndof = int(bins+0.001)-nparam;
  cout << "Chi-squared = " << chi2.getVal() << " ndof = " << ndof << " p-value = " << TMath::Prob(chi2.getVal(),ndof) << endl;
  cout << "Calculated reduced chi-squared of " << chi2.getVal()/double(ndof) << endl;

  // I suspect the fit and the RooChi2Var statistic are using different error definitions.
  // Perhaps one is based on the model value and the other on the data.
  // NO. From the error messages it looks like the function has been renormalized to the reduced intervals.
  // I also noticed an "Extended" option somewhere.

  // https://root.cern/doc/master/rf109__chi2residpull_8C.html

  // S h o w   r e s i d u a l   a n d   p u l l   d i s t s
  // -------------------------------------------------------
  // BEWARE - The comparison is between the data and the latest fit model plotted (could even be a component model ...)

  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist *hresid = frame->residHist();

  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = frame->pullHist();

  // Create a new frame to draw the residual distribution and add the distribution to the frame
  RooPlot *frame2 = x.frame(Title("Residual Distribution"));
  frame2->addPlotable(hresid, "P");

  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot *frame3 = x.frame(Title("Pull Distribution"));
  frame3->addPlotable(hpull, "P");

  TCanvas *c = new TCanvas("FitPlots", "FitPlots", 1600, 1200);
  // See https://root-forum.cern.ch/t/how-to-draw-pad-with-residuals-below-an-histogram/6820

  // try pad stuff
  TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);  // xlow, ylow, xup, yup
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
  pad1->SetBottomMargin(0.0225);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  const double LEFTMARGIN=0.15;
  pad1->SetLeftMargin(LEFTMARGIN);
  pad1->SetTicks(1,1);
  pad1->SetGrid();
  frame->SetMinimum(0.0);
  frame->SetMaximum(ymaxh);
  frame->SetMinimum(yminh);
  frame->SetTitle("Generator Level CoM Estimate #sqrt{s_{p}}");
  frame->GetXaxis()->SetTitle(xtitle.c_str());      // Use command line argument for this

  int ibinw = (binw*1000.0)+0.001;
  cout << "Bin width " << ibinw << " [MeV] " << endl;
  string ytitle = "Events per "+to_string(ibinw)+" MeV bin";
  cout << "ytitle: " << ytitle << endl;
  frame->GetYaxis()->SetTitle(ytitle.c_str());

  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleOffset(1.2);
  gStyle->SetTitleFont(42, ""); // set the pad title font
  gStyle->SetTitleFontSize(0.0725);

  double intlumi=100.0;
  double pxmin=0.495;
  double pdx=0.19;

  // Add chi-squared and degrees of freedom to plot
  int ndofdb=db.numEntries() -nparam;
  TPaveLabel *t1 = new TPaveLabel(pxmin,0.62,pxmin+pdx,0.75, Form("#chi^{2}/ndf = %6.1f / %3d", chi2.getVal(),ndofdb),"brNDC");
  t1->SetTextFont(53);
  t1->SetTextSize(24);
  frame->addObject(t1);

  frame->Draw();

  // Also need a legend ...
  const float lxmin=0.175;
  float lymin=0.26;
  const float ldx=0.32;
  const float ldy=0.22;
  TLegend* leg = new TLegend(lxmin,lymin,lxmin+ldx,lymin+ldy);
  leg->SetTextFont(42);
  //leg->SetHeader("#sqrt{s} = 250 GeV, #rho(E,z)=0.0","C");

  leg->SetHeader(Form("#sqrt{s} = %.1f GeV, e^{+} e^{-} #rightarrow #mu #mu (#gamma)",250.0),"C");
  TLegendEntry *entry=leg->AddEntry("data_obs","GP2X Data","ep");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(2);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1);
  entry->SetTextFont(42);

  TLegendEntry *entry2=leg->AddEntry("Fit","Fit","fl");
  entry2->SetLineColor(kBlue);
  entry2->SetLineStyle(1);
  entry2->SetLineWidth(2);
  entry2->SetMarkerColor(kBlue);
  entry2->SetMarkerStyle(20);
  entry2->SetMarkerSize(1);
  entry2->SetTextFont(42);

  leg->Draw();

  pad2->cd();
  pad2->SetLeftMargin(LEFTMARGIN);
  pad2->SetTicks(1,1);
  pad2->SetGrid();
  frame3->GetXaxis()->SetTitle(xtitle.c_str());
  frame3->GetYaxis()->SetTitle("(Data-Model)/Error");
  frame3->SetTitle("");
  frame3->SetMaximum( ymaxp);
  frame3->SetMinimum(-ymaxp);
  frame3->GetXaxis()->SetTitleOffset(1.2);
  frame3->GetYaxis()->SetTitleOffset(0.6);
  frame3->SetTitleSize(0.1);
  frame3->GetXaxis()->SetLabelSize(0.10);
  frame3->GetYaxis()->SetLabelSize(0.10);
  frame3->GetYaxis()->SetTitleSize(0.09);

  // Draw the pull distribution here now that we have overwritten the pull values with the values calculated as chi above
  frame3->Draw();
  //hres2->Draw("hist e2");
  c->SaveAs("MCspBeta.pdf");
  f->Close();
}

void ValSPDET(char* infile, Double_t MEAN, Double_t VAR)
{
  cout << "Running ValSP with  \n" << infile << endl;
  cout << MEAN << endl;
  cout << VAR << endl;
  //RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("maxSteps",100);
  //RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("epsRel",0.0000001);
  //RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setRealValue("epsAbs",0.0000001);

  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","FUMILI");
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSL","Levenberg-Marquardt");
  //ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.0000000001);
  //ROOT::Math::MinimizerOptions::SetDefaultPrecision(0.0000000001);
  //ROOT::Math::MinimizerOptions().SetMaxIterations(10000000);
  //ROOT::Math::MinimizerOptions().SetMaxFunctionCalls(10000000);
  ROOT::Math::MinimizerOptions().Print();
  
  //Plot the dimuon and electron boost (beta) values according to P/E
  //Create canvas and turn the stat box OFF
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  gStyle->SetOptStat(0);

  string xtitle="Center-of-Mass Energy, #sqrt{s_{p}} [GeV]";
  //string xtitle="Center";
  //Double_t ymaxh = 40000.0;
  Double_t yminh = 0.0;
  Double_t ymaxp = 5.5;
  //Double_t XMIN = MEAN - 10.0*VAR;
  //Double_t XMAX = MEAN + 3.0*VAR;
  Double_t XMIN = 242.0;
  Double_t XMAX = 251.0;
  Double_t bwidth = (XMAX-XMIN)/100.0;
  Double_t Escale_ini = MEAN;
  Double_t sigma_ini = 100.0*VAR/MEAN;
  Double_t msig = 0.0022 /2.0;
  Double_t psig = 0.001875 /2.0;
  Double_t dEhalf = 0.5*bwidth;
  Size_t msize = 1.5;

  gStyle->SetTitleFont(62, ""); // set the pad title font
  gStyle->SetMarkerSize(msize);
  gStyle->SetLegendTextSize(0.045);

  RooRealVar x("x","x",XMIN,XMAX,"GeV");             // Fit the beam energy in specified range (Last argument is the unit)

  
  TFile *f = new TFile(infile);
  TH1D *h1 = (TH1D*)f->Get("hs1");

  //h1->Scale(hd2->Integral()/h1->Integral());
  //for(Int_t i = 1; i < h1->GetNbinsX(); i++) h1->SetBinError(i,sqrt(h1->GetBinContent(i)));

  Double_t ymaxh = 1.5 * h1->GetMaximum();

  h1->Print();
  RooDataHist db("db", "db", x, Import(*h1));
  x.setRange("FitRange",XMIN,XMAX);

  RooRealVar Escale("E_{0}","Escale",Escale_ini,0.999*Escale_ini,1.001*Escale_ini,"GeV");           // Energy scale of full energy
  //RooRealVar Escale1("E_{1}","Escale1",Escale_ini,0.99*Escale_ini,1.01*Escale_ini,"GeV");           // Energy scale of full energy
  //RooRealVar Escale2("E2","Escale2",Escale_ini*0.99,0.98*Escale_ini,1.01*Escale_ini,"GeV");
  //RooRealVar Escale("E_{0}","Escale",Escale_ini,1.0*Escale_ini,1.00*Escale_ini,"GeV");           // Energy scale of full energy
  //RooRealVar Escale("E_{0}","Escale",Escale_ini,1.0*Escale_ini,1.0*Escale_ini,"GeV");           // Energy scale of full energy
  RooRealVar sigma("#sigma","sigma",sigma_ini,0.8*sigma_ini,1.25*sigma_ini,"%");                                       // Fractional resolution in per cent of Gaussian used in convolution
  RooRealVar s_peak("s_{peak}","s_{peak}",sigma_ini,0.95*sigma_ini,1.05*sigma_ini,"%");                                       // Fractional resolution in per cent of Gaussian used in convolution
  //RooRealVar m_peak("m_{peak}","m_{peak}",msig,0.99*msig,1.01*msig,"%");                                       // Fractional resolution in per cent of Gaussian used in convolution
  //RooRealVar p_peak("p_{peak}","p_{peak}",psig,0.99*psig,1.01*psig,"%");                                       // Fractional resolution in per cent of Gaussian used in convolution
  //RooRealVar sigma("#sigma","sigma",0.164,0.164,0.164,"%");                                       // Fractional resolution in per cent of Gaussian used in convolution
  //RooRealVar sigma("#sigma","sigma",sigma_ini,sigma_ini,sigma_ini,"%");
  // RooRealVar alpha("#alpha","alpha",15.0,5.0,25.0);
  //RooRealVar alpha("#alpha","alpha",24.3,0.1,28.0);// Alpha parameter of Beta distribution
  RooRealVar alpha("#alpha","alpha",1.0,0.5,35.0);
  //RooRealVar alpha("#alpha","alpha",28.0,28.0,28.0);// Alpha parameter of Beta distribution
  
  //RooRealVar alpha2("#alpha2","alpha2",15.0,5.0,25.0);                                                // Alpha parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",12.1,12.1,12.1);                                                // Alpha parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",11.66,11.66,11.66);                                                // Alpha parameter of Beta distribution
  //RooRealVar beta("#beta","beta",15.0,10.0,45.0);                                                // Alpha parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",0.20,0.15,0.4);                                                   // Beta parameter of Beta distribution
  //RooRealVar beta("#beta","beta",0.2695,0.2695,0.2695);                                                   // Beta parameter of Beta distribution
  RooRealVar beta("#beta","beta",0.2695,0.15,0.7);                                                   // Beta parameter of Beta distribution
  //RooRealVar beta("#beta","beta",0.33334,0.33334,0.33334);                                                   // Beta parameter of Beta distribution
  //RooRealVar beta2("#beta2","beta2",0.20,0.15,0.4);                                                   // Beta parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",0.11710645899206922,0.11710645899206922,0.11710645899206922);                                                // Alpha parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",0.1130971624410773335262428677792249,0.11309716244792249,0.11309716244792249);                                                // Alpha parameter of Beta distribution
  //RooRealVar alpha("#alpha","alpha",20,5,50);                                                // Alpha parameter of Beta distribution
  //RooRealVar beta("#beta","beta",0.2,0.2,0.10773335262428677);                                                   // Beta parameter of Beta distribution
  RooRealVar fpeak("f_{peak}","fpeak",0.11,0.01,0.35);                                               // Fraction of the Gaussian peak in the model
  //RooRealVar fpeak("f_{peak}","fpeak",0.2573,0.2573,0.2573);                                               // Fraction of the Gaussian peak in the model
  //RooRealVar fpeak2("f_{peak2}","fpeak2",0.35,0.1,0.95);                                               // Fraction of the Gaussian peak in the model

  //RooRealVar vXMAX("vXMAX","vXMAX",XMAX,XMAX,XMAX);

  //RooGenericPdf model("model","(vXMAX + 0.001 - x)^(alpha - 1.0) * (vXMAX + 1.001 - x)^(-1.0*alpha - beta)",RooArgList(x,vXMAX,alpha,beta));

  RooMyConvolvedBetaPdf tail("Convolved Beta","compiled class",x,Escale,alpha,beta,sigma);
  //RooMyConvolvedBetaPdf tail2("Convolved Beta","compiled class",x,Escale2,alpha2,beta2,sigma);
  RooMyNonStandardGaussianPdf peak("peak","peak pdf",x,Escale,s_peak);    // Note sigma here is the fractional resolution in percent
  //RooMyNonStandardGaussianPdf mpeak("mpeak","peak pdf",x,Escale,m_peak);    // Note sigma here is the fractional resolution in percent
  //RooMyNonStandardGaussianPdf ppeak("ppeak","peak pdf",x,Escale,p_peak);    // Note sigma here is the fractional resolution in percent
  //RooProdPdf peak("peak","peak pdf",mpeak,ppeak,0.0);
  RooAddPdf model("Model","model",RooArgList(peak,tail),fpeak);          // fpeak is the peak fraction
  //RooAddPdf model("Model","model",RooArgList(pmodel,tail2),fpeak2);          // fpeak is the peak fraction

  RooArgSet observables(model);

  // Plot data
  RooPlot* frame=x.frame(Title("Deconvolved Detector Level CoM Estimate #sqrt{s_{p}}"));
  db.plotOn(frame);

  RooChi2Var chi2("chi2","chi2",model,db,Range("FitRange"),
		  DataError(RooAbsData::Expected));   // RooChi2Var needs the fit range specified in calculating its statistic
  // Otherwise it may default to input binned dataset range that includes all of the TH1 bins
  // See https://sft.its.cern.ch/jira/browse/ROOT-10038
  double chi2f{};
  double binw{};
  double rchisq{};
  double nevts{};
  int nparam{};
  double bins{};
  int ndof{};

  cout << "Initial Chi-squared " << chi2.getVal() << endl;
  //     RooFitResult *result = model.chi2FitTo(db,Range("FitRange"),Hesse(false),Minos(false),Save(),IntegrateBins(1.0e-12));

  //model.chi2FitTo(db,Range("FitRange"));

  model.fitTo(db,Range("FitRange"));

  cout << "Fitted Chi-squared " << setprecision(10) << chi2.getVal() << endl;
  chi2f = chi2.getVal();
  cout << "Escale  " << setprecision(10) << Escale.getValV()  << " +- " << Escale.getError() << " Fractional error " << Escale.getError()/Escale.getValV() << endl;
  //cout << "Escale  " << setprecision(10) << Escale1.getValV()  << " +- " << Escale1.getError() << " Fractional error " << Escale1.getError()/Escale1.getValV() << endl;
  cout << "sigma " << setprecision(10) << sigma.getValV() << " +- " << sigma.getError() << endl;
  cout << "alpha " << setprecision(10) << alpha.getValV() << " +- " << alpha.getError() << endl;
  cout << "beta " << setprecision(10) << beta.getValV() << " +- " << beta.getError() << endl;
  cout << "fpeak " << setprecision(10) << fpeak.getValV() << " +- " << fpeak.getError() << endl;
  cout << "s_peak " << setprecision(10) << s_peak.getValV() << " +- " << s_peak.getError() << endl;
  //cout << "m_sig " << setprecision(10) << m_peak.getValV() << " +- " << m_peak.getError() << endl;
  //cout << "p_sig " << setprecision(10) << p_peak.getValV() << " +- " << p_peak.getError() << endl;

  //sigma.setConstant(kTRUE);
  alpha.setConstant(kTRUE);
  beta.setConstant(kTRUE);
  fpeak.setConstant(kTRUE);
  //Escale1.setConstant(kTRUE);
  //s_peak.setConstant(kTRUE);

  model.fitTo(db,Range("FitRange"));

  cout << "Fitted Chi-squared " << setprecision(10) << chi2.getVal() << endl;
  chi2f = chi2.getVal();
  cout << "Escale  " << setprecision(10) << Escale.getValV()  << " +- " << Escale.getError() << " Fractional error " << Escale.getError()/Escale.getValV() << endl;
  //cout << "Escale1  " << setprecision(10) << Escale1.getValV()  << " +- " << Escale1.getError() << " Fractional error " << Escale1.getError()/Escale1.getValV() << endl;
  cout << "sigma " << setprecision(10) << sigma.getValV() << " +- " << sigma.getError() << endl;
  cout << "alpha " << setprecision(10) << alpha.getValV() << " +- " << alpha.getError() << endl;
  cout << "beta " << setprecision(10) << beta.getValV() << " +- " << beta.getError() << endl;
  cout << "fpeak " << setprecision(10) << fpeak.getValV() << " +- " << fpeak.getError() << endl;
  //cout << "m_sig " << setprecision(10) << m_peak.getValV() << " +- " << m_peak.getError() << endl;
  //cout << "p_sig " << setprecision(10) << p_peak.getValV() << " +- " << p_peak.getError() << endl;
  cout << "s_peak " << setprecision(10) << s_peak.getValV() << " +- " << s_peak.getError() << endl;

  sigma.setConstant(kTRUE);
  //alpha.setConstant(kTRUE);
  //beta.setConstant(kTRUE);
  //fpeak.setConstant(kTRUE);
  //Escale1.setConstant(kTRUE);
  s_peak.setConstant(kTRUE);

  model.fitTo(db,Range("FitRange"));

  cout << "Fitted Chi-squared " << setprecision(10) << chi2.getVal() << endl;
  chi2f = chi2.getVal();
  cout << "Escale  " << setprecision(10) << Escale.getValV()  << " +- " << Escale.getError() << " Fractional error " << Escale.getError()/Escale.getValV() << endl;
  //cout << "Escale1  " << setprecision(10) << Escale1.getValV()  << " +- " << Escale1.getError() << " Fractional error " << Escale1.getError()/Escale1.getValV() << endl;
  cout << "sigma " << setprecision(10) << sigma.getValV() << " +- " << sigma.getError() << endl;
  cout << "alpha " << setprecision(10) << alpha.getValV() << " +- " << alpha.getError() << endl;
  cout << "beta " << setprecision(10) << beta.getValV() << " +- " << beta.getError() << endl;
  cout << "fpeak " << setprecision(10) << fpeak.getValV() << " +- " << fpeak.getError() << endl;
  //cout << "m_sig " << setprecision(10) << m_peak.getValV() << " +- " << m_peak.getError() << endl;
  //cout << "p_sig " << setprecision(10) << p_peak.getValV() << " +- " << p_peak.getError() << endl;
  cout << "s_peak " << setprecision(10) << s_peak.getValV() << " +- " << s_peak.getError() << endl;

  //compute the chi2 the correct way...
  TH1D *hmodel = (TH1D*)model.createHistogram("hmodel",x,RooFit::Binning(25,XMIN,XMAX));
  TH1D *hdb = (TH1D*)db.createHistogram("hdb",x,RooFit::Binning(25,XMIN,XMAX));
  //hmodel->Scale(500000.0/hmodel->Integral());
  //hdb->Scale(500000.0/hdb->Integral());
  //for (Int_t i = 1; i <= hmodel->GetNbinsX(); i++) hmodel->SetBinError(i,sqrt(hmodel->GetBinContent(i)));
  //for (Int_t i = 1; i <= hdb->GetNbinsX(); i++) hdb->SetBinError(i,sqrt(hdb->GetBinContent(i)) + 0.01*hdb->GetBinContent(i));
  //hmodel->Scale(h1->GetBinContent(10)/hmodel->GetBinContent(5));
  cout << "model hist check " << hmodel->GetBinContent(5) << endl;
  cout << "data hist check " << h1->GetBinContent(5) << endl;
  cout << "data err check " << h1->GetBinError(5) << endl;
  Double_t hchi2 = 0.0;
  TH1D *hres2 = (TH1D*)hdb->Clone();
  for (Int_t i = 0; i <= hmodel->GetNbinsX(); i++) 
    {
      hres2->SetBinContent(i,0);
      if (hmodel->GetBinCenter(i) > XMIN && hmodel->GetBinCenter(i) < XMAX) 
	{
	  hchi2 += pow(hmodel->GetBinContent(i) - hdb->GetBinContent(i),2.0)/(pow(hdb->GetBinError(i),2.0) + hmodel->GetBinContent(i));
	  hres2->SetBinContent(i,(hmodel->GetBinContent(i) - hdb->GetBinContent(i))/sqrt(pow(hdb->GetBinError(i),2.0) + hmodel->GetBinContent(i)));
	}
    }

  cout << "Histogram computed chi2 " << hchi2 << endl;

  model.plotOn(frame,Range("FitRange"),NormRange("FitRange"),Components(tail),LineStyle(kSolid),LineColor(kMagenta));
  model.plotOn(frame,Range("FitRange"),NormRange("FitRange"),Components(peak),LineStyle(kSolid),LineColor(kRed));
  model.plotOn(frame,Range("FitRange"),NormRange("FitRange"));
  model.paramOn(frame,Label(""),Parameters(RooArgSet(Escale,sigma,alpha,beta,fpeak,s_peak)),Layout(0.17,0.485,0.89),Format("NEU",AutoPrecision(2)));

  RooArgSet *flparams = model.getParameters(observables);
  nparam = (flparams->selectByAttrib("Constant",kFALSE))->getSize();
  cout << "nparam = " << nparam << endl;
  rchisq = frame->chiSquare(nparam);
  cout << "Reduced chi-squared?? " << rchisq << endl;
  cout << "Inferred dof " << chi2.getVal()/rchisq << endl;
  nevts=frame->getFitRangeNEvt();
  cout << "nevts = " << nevts << endl;
  binw=frame->getFitRangeBinW();
  cout << "binw = " << binw << endl;
  bins = (XMAX-XMIN)/binw;
  cout << "Number of bins in fit " << bins << endl;
  ndof = int(bins+0.001)-nparam;
  cout << "Chi-squared = " << chi2.getVal() << " ndof = " << ndof << " p-value = " << TMath::Prob(chi2.getVal(),ndof) << endl;
  cout << "Calculated reduced chi-squared of " << chi2.getVal()/double(ndof) << endl;

  // I suspect the fit and the RooChi2Var statistic are using different error definitions.
  // Perhaps one is based on the model value and the other on the data.
  // NO. From the error messages it looks like the function has been renormalized to the reduced intervals.
  // I also noticed an "Extended" option somewhere.

  // https://root.cern/doc/master/rf109__chi2residpull_8C.html

  // S h o w   r e s i d u a l   a n d   p u l l   d i s t s
  // -------------------------------------------------------
  // BEWARE - The comparison is between the data and the latest fit model plotted (could even be a component model ...)

  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist *hresid = frame->residHist();

  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist *hpull = frame->pullHist();

  // Create a new frame to draw the residual distribution and add the distribution to the frame
  RooPlot *frame2 = x.frame(Title("Residual Distribution"));
  frame2->addPlotable(hresid, "P");

  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot *frame3 = x.frame(Title("Pull Distribution"));
  frame3->addPlotable(hpull, "P");

  TCanvas *c = new TCanvas("FitPlots", "FitPlots", 1600, 1200);
  // See https://root-forum.cern.ch/t/how-to-draw-pad-with-residuals-below-an-histogram/6820

  // try pad stuff
  TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);  // xlow, ylow, xup, yup
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
  pad1->SetBottomMargin(0.0225);
  pad1->SetBorderMode(0);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->SetBorderMode(0);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  const double LEFTMARGIN=0.15;
  pad1->SetLeftMargin(LEFTMARGIN);
  pad1->SetTicks(1,1);
  pad1->SetGrid();
  frame->SetMinimum(0.0);
  frame->SetMaximum(ymaxh);
  frame->SetMinimum(yminh);
  frame->SetTitle("Deconvolved Detector Level CoM Estimate #sqrt{s_{p}}");
  frame->GetXaxis()->SetTitle(xtitle.c_str());      // Use command line argument for this

  int ibinw = (binw*1000.0)+0.001;
  cout << "Bin width " << ibinw << " [MeV] " << endl;
  string ytitle = "Events per "+to_string(ibinw)+" MeV bin";
  cout << "ytitle: " << ytitle << endl;
  frame->GetYaxis()->SetTitle(ytitle.c_str());

  frame->GetXaxis()->SetTitleOffset(1.2);
  frame->GetYaxis()->SetTitleOffset(1.2);
  gStyle->SetTitleFont(42, ""); // set the pad title font
  gStyle->SetTitleFontSize(0.0725);

  double intlumi=100.0;
  double pxmin=0.495;
  double pdx=0.19;

  // Add chi-squared and degrees of freedom to plot
  int ndofdb=db.numEntries() -nparam;
  TPaveLabel *t1 = new TPaveLabel(pxmin,0.62,pxmin+pdx,0.75, Form("#chi^{2}/ndf = %6.1f / %3d", chi2.getVal(),ndofdb),"brNDC");
  t1->SetTextFont(53);
  t1->SetTextSize(24);
  frame->addObject(t1);

  frame->Draw();

  // Also need a legend ...
  const float lxmin=0.175;
  float lymin=0.26;
  const float ldx=0.32;
  const float ldy=0.22;
  TLegend* leg = new TLegend(lxmin,lymin,lxmin+ldx,lymin+ldy);
  leg->SetTextFont(42);
  //leg->SetHeader("#sqrt{s} = 250 GeV, #rho(E,z)=0.0","C");

  leg->SetHeader(Form("#sqrt{s} = %.1f GeV, e^{+} e^{-} #rightarrow #mu #mu (#gamma)",250.0),"C");
  TLegendEntry *entry=leg->AddEntry("data_obs","GP2X Data","ep");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(2);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(20);
  entry->SetMarkerSize(1);
  entry->SetTextFont(42);

  TLegendEntry *entry2=leg->AddEntry("Fit","Fit","fl");
  entry2->SetLineColor(kBlue);
  entry2->SetLineStyle(1);
  entry2->SetLineWidth(2);
  entry2->SetMarkerColor(kBlue);
  entry2->SetMarkerStyle(20);
  entry2->SetMarkerSize(1);
  entry2->SetTextFont(42);

  leg->Draw();

  pad2->cd();
  pad2->SetLeftMargin(LEFTMARGIN);
  pad2->SetTicks(1,1);
  pad2->SetGrid();
  frame3->GetXaxis()->SetTitle(xtitle.c_str());
  frame3->GetYaxis()->SetTitle("(Data-Model)/Error");
  frame3->SetTitle("");
  frame3->SetMaximum( ymaxp);
  frame3->SetMinimum(-ymaxp);
  frame3->GetXaxis()->SetTitleOffset(1.2);
  frame3->GetYaxis()->SetTitleOffset(0.6);
  frame3->SetTitleSize(0.1);
  frame3->GetXaxis()->SetLabelSize(0.10);
  frame3->GetYaxis()->SetLabelSize(0.10);
  frame3->GetYaxis()->SetTitleSize(0.09);

  // Draw the pull distribution here now that we have overwritten the pull values with the values calculated as chi above
  frame3->Draw();
  //hres2->Draw("hist e2");
  c->SaveAs("DETspBeta.pdf");
  f->Close();
}

void DeconVal(char *infile, char *infileDECON)
{
  //gROOT->ProcessLine(" .x RooMyConvolvedBetaPdf.cxx ");
  //gROOT->ProcessLine(" .x RooMyNonStandardGaussianPdf.cxx ");
  //gROOT->ProcessLine(" .x ValMac.C ");
  //gROOT->ProcessLine(" .x GP2X_v3.C ");
  //ValidationPlots(infile);
  //ValidationPlots(infile);
  //Pass true to use the deconvolved data file settings
  ValSPDET(infileDECON,250.0,0.379);
  //Pass false to use the normal data settings
  ValSPMC(infile,250.0,0.379);
}

