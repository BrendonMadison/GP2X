#this uses the pyroot as in the py2.7 installation native to the b047m507 user
from ROOT import TFile, TH1, TH1D
from math import sqrt
import pickle as pkl
import sys

#default settings are (for dimuons at 100fb-1)
#hs1 200 242 251 550000
#the histogram settings much match those that come out of the deconovlution program!
def RootifySVG(INAME,RNAME,HNAME,NBINS,LEND,HEND,NDATA):
    inf = open(INAME,"rb")
    data = pkl.load(inf)
    h = TH1D(HNAME,HNAME,NBINS,LEND,HEND)

    for i in range(len(data)):
        h.SetBinContent(i+1,data[i])
    h.Scale(NDATA/h.Integral())
    for i in range(len(data)):
        h.SetBinError(i,sqrt(h.GetBinContent(i)))
    h.Rebin(2)

    #output to file
    f = TFile.Open(RNAME,"RECREATE")
    h.Write()
    f.Close()


RootifySVG(sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]),float(sys.argv[5]),float(sys.argv[6]),float(sys.argv[7]))
#RootifySVG("OutSVGDatMC.pkl","DatMC.root","hd1")
#RootifySVG("OutSVGRefDet.pkl","RefDet.root","hs2")
#RootifySVG("OutSVGRefMC.pkl","RefMC.root","hd2")
