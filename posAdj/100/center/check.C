#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <sstream>

void check()
{

  cout << "Starting Calculation" << endl;

  TRandom3* randomGen = new TRandom3(12191982);

  TFile* file = new TFile("complete.root");

  TTree* tree = (TTree *)file->Get("Signal");

  int nEvents = tree->GetEntries();
  
  double addresses[1225] = {};
  for (int k=0; k<1225; k++) {
    std::stringstream ss2;
    ss2 << k;
    string str = "Crystal_"+ss2.str();
    const char* charstr = str.c_str();
    tree->SetBranchAddress(charstr, &addresses[k]);
  }

  cout << "The number of events registered is: " << tree->GetEntries() << endl;
  
  TH1D* totalEnergy = new TH1D("totalEnergy", "total energy", 50, 0, 100);
  
  double total = 0.;
  
  for (int i =0; i < nEvents; i++)
    { total = 0;
      tree->GetEntry(i);
      for (int j=0; j< 1225; j++)
	{total+= addresses[j];}
      totalEnergy->Fill(total);
    }

  totalEnergy->Draw();

  

}
