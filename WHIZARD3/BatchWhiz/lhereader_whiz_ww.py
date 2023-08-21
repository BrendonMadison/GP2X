import ROOT as r
from array import array
import sys
from locale import atof
from math import sqrt

#print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

ARGSQS = atof(sys.argv[1])
ARGIFILE = str(sys.argv[2])
ARGOFILE = str(sys.argv[3])
ARGXSEC = atof(sys.argv[4])
ARGXSECERR = atof(sys.argv[5])

def convert_lhe(sqs, xsecin, xsecerrin, fname_in, fname_out="lhe_example.root"):

    MME = 0.00051099895000
    MMU = 0.1056583755
    #get the cross-section
    #The following should be in picobarns
    #WHIZARD uses femtobarns and KKMC uses picobarns
    xs = xsecin/1000.0
    xe = xsecerrin/1000.0
    print(xs)
    print(xe)
    print(type(xs))
    print(type(xe))
    events = []
    event = ""
    in_event = False
    with open(fname_in, "r") as fhin:
        iline = 0
        for line in fhin:

            if in_event and line.startswith("<"):
                in_event = False
                events.append(event)
                event = ""

            if in_event:
                event += line

            if line.startswith("<event>"):
                in_event = True

        if event:
            events.append(event)

    f1 = r.TFile(fname_out, "recreate")
    t1 = r.TTree("t","t")

    insqs = array( 'f', [ 0.0 ] )
    inxs = array( 'f', [ 0.0 ] )    
    inxe = array( 'f', [ 0.0 ] )
    bevent = array( 'i', [ 0 ] )
    PartNum = array( 'i', [ 0 ] )
    GamNum = array( 'i', [ 0 ] )
    
    Elep4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    Posp4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    Mup4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    AMup4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    Nup4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    ANup4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    Wp4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    AWp4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    Gam1p4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    Gam2p4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    Gam3p4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    Gam4p4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    Gam5p4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    Gam6p4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    GamTotp4 = r.TLorentzVector(0.0,0.0,0.0,0.0)
    
    t1.Branch("PartNum",PartNum, "PartNum/I")
    t1.Branch("GamNum",GamNum, "GamNum/I")
    t1.Branch("sqs",insqs, "sqs/F")
    t1.Branch("xsec",inxs, "xsec/F")
    t1.Branch("xsecerr",inxe, "xsecerr/F")
    t1.Branch("event",bevent, "event/I")
    t1.Branch("Elep4.","TLorentzVector",Elep4)
    t1.Branch("Posp4.","TLorentzVector",Posp4)
    t1.Branch("Mup4.","TLorentzVector",Mup4)
    t1.Branch("AMup4.","TLorentzVector",AMup4)
    t1.Branch("Nup4.","TLorentzVector",Nup4)
    t1.Branch("ANup4.","TLorentzVector",ANup4)
    t1.Branch("Wp4.","TLorentzVector",Wp4)
    t1.Branch("AWp4.","TLorentzVector",AWp4)
    t1.Branch("Gam1p4.","TLorentzVector",Gam1p4)
    t1.Branch("Gam2p4.","TLorentzVector",Gam2p4)
    t1.Branch("Gam3p4.","TLorentzVector",Gam3p4)
    t1.Branch("Gam4p4.","TLorentzVector",Gam4p4)
    t1.Branch("Gam5p4.","TLorentzVector",Gam5p4)
    t1.Branch("Gam6p4.","TLorentzVector",Gam6p4)
    t1.Branch("GamTotp4.","TLorentzVector",GamTotp4)
    
    for ievt,evt in enumerate(events):
        particle_lines = evt.splitlines()[1:]
        pcnt = int(0)
        gcnt = int(0)
        ecnt = int(0)
        acnt = int(0)
        insqs[0] = sqs
        inxs[0] = xs
        inxe[0] = xe
        sumgamp4 = [0.0,0.0,0.0,0.0]
        Elep4.SetXYZM(0.0,0.0,0.0,0.0)
        Posp4.SetXYZM(0.0,0.0,0.0,0.0)
        Mup4.SetXYZM(0.0,0.0,0.0,0.0)
        AMup4.SetXYZM(0.0,0.0,0.0,0.0)
        Nup4.SetXYZM(0.0,0.0,0.0,0.0)
        ANup4.SetXYZM(0.0,0.0,0.0,0.0)
        Wp4.SetXYZM(0.0,0.0,0.0,0.0)
        AWp4.SetXYZM(0.0,0.0,0.0,0.0)
        Gam1p4.SetXYZM(0.0,0.0,0.0,0.0)
        Gam2p4.SetXYZM(0.0,0.0,0.0,0.0)
        Gam3p4.SetXYZM(0.0,0.0,0.0,0.0)
        Gam4p4.SetXYZM(0.0,0.0,0.0,0.0)
        Gam5p4.SetXYZM(0.0,0.0,0.0,0.0)
        Gam6p4.SetXYZM(0.0,0.0,0.0,0.0)
        GamTotp4.SetXYZM(0.0,0.0,0.0,0.0)
        for particle_line in particle_lines:
            parts = particle_line.split()
            evt_pdgid = int(parts[0])
            #evt_status = int(parts[1])
            #evt_parent1 = int(parts[2])
            #evt_parent2 = int(parts[3])
            #evt_color1 = int(parts[4])
            #evt_color2 = int(parts[5])
            evt_px, evt_py, evt_pz, evt_e = map(float,parts[6:10])
            #evt_mass = float(parts[10])
            #evt_spin = float(parts[11])
            bevent = ievt
            if evt_pdgid == 13:
                Mup4.SetXYZM(evt_px,evt_py,evt_pz,MMU)
            if evt_pdgid == -13:
                AMup4.SetXYZM(evt_px,evt_py,evt_pz,MMU)
            if evt_pdgid == 14:
                Nup4.SetXYZM(evt_px,evt_py,evt_pz,0.0)
            if evt_pdgid == -14:
                ANup4.SetXYZM(evt_px,evt_py,evt_pz,0.0)
            if evt_pdgid == 11 and ecnt == 1:
                Mup4.SetXYZM(evt_px,evt_py,evt_pz,MME)
            elif evt_pdgid == 11 and ecnt == 0:
                Elep4.SetXYZM(evt_px,evt_py,evt_pz,MME)
                ecnt = 1
            if evt_pdgid == -11 and acnt == 1:
                AMup4.SetXYZM(evt_px,evt_py,evt_pz,MME)
            elif evt_pdgid == -11 and acnt == 0:
                Posp4.SetXYZM(evt_px,evt_py,evt_pz,MME)
                acnt = 1
            if evt_pdgid == 22:
                GamTotp4.SetXYZM(evt_px + GamTotp4.Px(),evt_py + GamTotp4.Py(),evt_pz + GamTotp4.Pz(),sqrt(abs((GamTotp4.E()+evt_e)**2 - (evt_px + GamTotp4.Px())**2 - (evt_py + GamTotp4.Py())**2 - (evt_pz + GamTotp4.Pz())**2)))
                gcnt += 1
            pcnt += 1
        PartNum[0] = pcnt
        GamNum[0] = gcnt
        Wp4.SetXYZT(Mup4.Px() + ANup4.Px(),Mup4.Py() + ANup4.Py(),Mup4.Pz() + ANup4.Pz(),Mup4.E()+ANup4.E())
        AWp4.SetXYZT(AMup4.Px() + Nup4.Px(),AMup4.Py() + Nup4.Py(),AMup4.Pz() + Nup4.Pz(),AMup4.E()+Nup4.E())
        #GamM = sqrt((sqs+0.0001)**2 - 2.0*(MME**2) - 2.0*(Mup4.E()*AMup4.E() - Mup4.Px()*AMup4.Px() - Mup4.Py()*AMup4.Py() - Mup4.Pz()*AMup4.Pz()) - 2.0*(Mup4.E()*(sqs-Mup4.E()-AMup4.E()) + Mup4.Px()*(Mup4.Px() + AMup4.Px()) +  Mup4.Py()*(Mup4.Py() + AMup4.Py()) +  Mup4.Pz()*(Mup4.Pz() + AMup4.Pz())) - 2.0*(AMup4.E()*(sqs-Mup4.E()-AMup4.E()) + AMup4.Px()*(Mup4.Px() + AMup4.Px()) +  AMup4.Py()*(Mup4.Py() + AMup4.Py()) +  AMup4.Pz()*(Mup4.Pz() + AMup4.Pz())))
        #GamTotp4.SetXYZM(-1.0*(Mup4.Px()+AMup4.Px()),-1.0*(Mup4.Py()+AMup4.Py()),-1.0*(Mup4.Pz()+AMup4.Pz()),GamM)
        t1.Fill()
    t1.Print()
    t1.Write()
    f1.Close()

if __name__ == "__main__":

    #convert_lhe(250.0,"KKMCResults/RunTEST/LHE_OUT.LHE.500.0")
    convert_lhe(ARGSQS,ARGXSEC,ARGXSECERR,ARGIFILE,ARGOFILE)
