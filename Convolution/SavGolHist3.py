from ROOT import TH1, TH1D, TFile, TTree, gDirectory
import pickle as pkl
import sys

def Picklify(FNAME,HNAME,PNAME):
    #f = TFile.Open("HistMC_SVGTest_sp.root")
    #hs1 = f.Get("HistMC_SVGTest_sp")
    f = TFile.Open(FNAME)
    hs1 = f.Get(HNAME)
    print("Done loading in histograms!")

    hs1arr = hs1.GetArray()
    #output data list
    datas1 = []
    #for i in range(hs1.GetNbinsX()):
    for i in range(int(hs1.GetNbinsX())):
        print(hs1.GetArray()[i+1])
        datas1.append(float(hs1.GetArray()[i+1]))
    print(datas1)
    print("Done converting histograms to lists!")
    #output pickle file
    #outf = open('OutharrSVGTest.pkl','wb')
    outf = open(PNAME,'wb')
    pkl.dump(datas1,outf)
    outf.close()

fname = sys.argv[1]
hname = sys.argv[2]
pname = sys.argv[3]

Picklify(fname,hname,pname)

print("Wrote" + str(fname) + " to pickle file: " + str(pname))
