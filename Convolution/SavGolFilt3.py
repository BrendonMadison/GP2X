#this must use python3 as installed at 
#module load python
import pickle as pkl
import sys
import numpy as np
from scipy.signal import savgol_filter as svg

#default settings are 71 3 22
def SavGolify(INAME,ONAME,WINSIZE,POLYORDER,FROMPEAK):
    #infile with the detector data to be savgol filtered
    inf = open(INAME,"rb")
    data = pkl.load(inf)
    #print(np.argmax(data))
    inf.close()    
    #original test using range of [0.01,500] GeV for sqrt(s_p) with 100000 bins
    #x = np.linspace(0.01,500,len(data))
    print(np.argmax(data))
    #pass1 = svg(data,5,1,mode='nearest')
    #pass1 = svg(data,11,1,mode='nearest')
    svgdata = data[0:np.argmax(data)-FROMPEAK]
    #svgdata = data[0:np.argmax(pass1)-FROMPEAK]
    #svgdata = svg(data,WINSIZE,POLYORDER,mode='nearest')
    print(len(svgdata))
    print(svgdata[-1])
    svgdata = svg(svgdata,WINSIZE,POLYORDER,mode='nearest')
    svgdata = svg(svgdata,WINSIZE,POLYORDER,mode='nearest')
    print(len(svgdata))
    print(svgdata[-1])
    #svgdata2 = svg(data[np.argmax(data)-44:np.argmax(data)-25],11,3)
    #svgdata = svg(svgdata,61,2)
    #svgdata = svg(svgdata,71,4)
    #svgdata = np.append(svgdata,svgdata2)
    svgdata = np.append(svgdata[:-5],(svgdata[-5]+data[np.argmax(data)-FROMPEAK-4])*0.5)
    svgdata = np.append(svgdata,data[np.argmax(data)-FROMPEAK-4:len(data)])
    #pkdata = svg(data[np.argmax(data)-FROMPEAK-4:len(data)],21,4,mode='nearest')
    #pkdata = svg(pkdata,21,4,mode='nearest')
    #svgdata = np.append(svgdata[:-5],(svgdata[-5]+pkdata[0])*0.5)
    #svgdata = np.append(svgdata,pkdata[2:len(pkdata))
    svglist = svgdata.tolist()
    #print(svglist)
    print(np.argmax(svglist))
    print(len(svglist))
    #output pickle file to read back into root
    outf = open(ONAME,"wb")
    pkl.dump(svglist,outf,protocol=2)
    outf.close()

pname = sys.argv[1]
svgname = sys.argv[2]
wsize = int(sys.argv[3])
porder = int(sys.argv[4])
frmpk = int(sys.argv[5])

SavGolify(pname,svgname,wsize,porder,frmpk)
#SavGolify("OutharrDatMC.pkl","OutSVGDatMC.pkl")
#SavGolify("OutharrRefDet.pkl","OutSVGRefDet.pkl")
#SavGolify("OutharrRefMC.pkl","OutSVGRefMC.pkl")



