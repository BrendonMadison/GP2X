void DECONSVG_VALM(char *INFILE, char *DINFILE)
{
  gROOT->ProcessLine(" .x /panfs/pfs.local/work/wilson/b047m507/GP2X_dev/RooMyConvolvedBetaPdf.cxx ");
  gROOT->ProcessLine(" .x /panfs/pfs.local/work/wilson/b047m507/GP2X_dev/RooMyNonStandardGaussianPdf.cxx ");
  gROOT->ProcessLine(Form(" .x /panfs/pfs.local/work/wilson/b047m507/GP2X_dev/DeconValM.C(\"%s\",\"%s\")",INFILE,DINFILE));
}
