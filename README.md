# GuineaPig To X (GP2X) (readme under construction for v2)
With X being a difermion event generator . So far confirmed works with KKMC, BHWIDE, WHIZARD v3

## Relevant research published here: https://arxiv.org/abs/2308.09676

Takes initial beam dynamics from GuineaPig and uses them to boost the final state difermions according to their center of mass energies. Main use is for precision generation of dimuons or Bhabha scattered electrons and positrons.

## Requires: 

ROOT v6 , std library of C to run the .C files  
A computing cluster with SLURM to run the .sh file if you want to run parallel jobs  

## Includes:

MakeLumiRoot.C
-- Takes a GuineaPig luminosity file, typically called lumi.ee.out , and converts it into a ROOT file that can be used with GP2X.C

GP2X.C
-- Takes a ROOT file of format given by MakeLumiRoot.C and a difermion ROOT file with the following tree branches:  
----"Mup4." of type TLorentzVector  
----"AMup4." of type TLorentzVector  
----"GamTotp4." of type TLorentzVector  
----"sqs" of type Float_t  
----"xsec" of type Float_t  
--Outputs a ROOT file that has randomly generated final state particles that are boosted into the GuineaPig beam particle's momenta  
--Includes option for beam energy spread (BES) as well as spreading in the beam vertex.  
--Struct ParticleTrack that stores boosted particle four vector and has functions for smearing if it is measured in the ECAL or particle tracker

runGP2X.sh
-- Runs GP2X in parallel setting so you can generate millions of output events  

Also includes an input data files of 100klumi.ee.out (a lumi file). There is no difermion file as they are simply too large to attach to github. Using a difermion file and a decent computer you should see about 100k events every 10 minutes. The code has been optimized about 20 times with various tricks that you can see in GP2X.C . It used to take days to get 100k events...


### To Do:
1.) Include tracker resolution effects on difermion tracks (FINISHED)
2.) Include ECAL resolution effects on photons (ALMOST FINISHED)
3.) Include displaced vertex effects on tracker measurements (if any measurable)  
4.) Add support for hadronic difermion final states  
5.) Implement (optional) diagnostic fits of the final output  
