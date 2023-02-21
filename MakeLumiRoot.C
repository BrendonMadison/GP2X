#include <iostream>
#include <string>
#include <fstream>
using namespace std;

//Function for cutting the lumi file to a certain number of lines
void CutLumiFile(ifstream& inLumi, ofstream& outLumi, Int_t NumLines)
{
  string tmp;
  for (Int_t i = 0; i < NumLines; i++)
    {
      getline(inLumi,tmp);
      outLumi << tmp << endl;
    }
}

//For of MakeLumiRoot for when the user simply provides the GuineaPig Lumi file and the Output file name
void MakeLumiRoot(char *gname,char *oname)
{
  auto f = TFile::Open(Form("%s.root",oname),"RECREATE");
  auto t = new TTree("t","");
  t->ReadFile(gname,"E1/F:E2/F:PosZ/F:PosY/F:PosX/F:D:x1/F:y1/F:x2/F:y2/F:AA:BB:CC:AB:BA:AC:BC",' ');
  f->Write();
  f->Close();
}

//For of MakeLumiRoot for when the user wants to cut the Lumi file to a certain length (can save time of random sampling later on)
 void MakeLumiRoot(char *gname,char *oname, Int_t nlines)
{
  //Open the file stream objects
  ifstream inputLumi;
  inputLumi.open(gname);
  ofstream outputLumi;
  outputLumi.open(Form("%s_%i.lumi",gname,nlines));
  //Cut the lumi file down
  CutLumiFile(inputLumi,outputLumi,nlines);
  //Close those files!
  outputLumi.close();
  inputLumi.close();
  //Open in root and turn it into a tree!
  //Huh... doesn't that usually go the other way around?
  //Who would've guessed that you could chop something down and create a tree.
  auto f = TFile::Open(Form("%s.root",oname),"RECREATE");
  auto t = new TTree("t","");
  t->ReadFile(Form("%s_%i.lumi",gname,nlines),"E1/F:E2/F:PosZ/F:PosY/F:PosX/F:D:x1/F:y1/F:x2/F:y2/F:AA:BB:CC:AB:BA:AC:BC",' ');
  f->Write();
  f->Close();
}
