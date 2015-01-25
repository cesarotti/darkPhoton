#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>

/*
 * THINGS TO CHECK BEFORE YOU RUN
 * 1. Energy
 * 2. Position
 * 3. Binning
 * 4. Title
 */

/*
  This program is used for checking how different noise and number of
  crystals affects the energy resolution. Doesn't check for several photons.
  "True" energy is now calculated by the total energy recovered from the
  Geant simulation.
*/



const double WEIGHT = 1./1000.;

// Changes x and y coordinates in crystal units to the crystal ID
int CoordtoID(int x, int y)
{
  return(y+17)*35+(x+17);
}
 

// Gets x position from crystalID
int IDtoX(int crystalID)
{
  return crystalID%35-17;
}


//Gets y position from crystalID
int IDtoY(int crystalID)
{
  return crystalID/35-17;
}




void spreadPlots() 
{ 

  cout << "Starting plots..." << endl;

  TFile* file = new TFile("complete.root");

  TTree* tree = (TTree *)file->Get("Signal");



  int nEvents = tree->GetEntries();

   double addresses[1225] = {};
  for (int k=0; k<1225; k++){
    std::stringstream ss2;
    ss2 << k; 
    string str = "Crystal_"+ss2.str();
    const char* charstr = str.c_str(); 
    tree->SetBranchAddress(charstr, &addresses[k]);
  }


  TH2D* spread = new TH2D("spread", "Energy_Spread", 12, 6, 18, 20, -10, 10);
  TH1D* energyColl = new TH1D("energy", "Energy_Recovered", 26, 15, 40);
  double totEnergy = 0.0;
  //iterate through all events
  for (int i = 0; i < nEvents; i++)
    {
      tree->GetEntry(i);
      totEnergy= 0.;
      for(int w = 0; w < 1225; w++)
	{
	  
	  spread->Fill(IDtoX(w), IDtoY(w), addresses[w]*WEIGHT);
	  totEnergy+= addresses[w];
	  
	}
      energyColl->Fill(totEnergy);
    }

  for (int x=0; x<25; x++){
    for (int y=-9; y<25; y++){
      if (spread->GetBinContent(x, y)<.2)
	{spread->SetBinContent(x, y, 0.);}
    }
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 1000, 500);
  canvas->Divide(2,1);
  canvas->cd(1);
  spread->Draw("BOX");
  canvas->cd(2);
  spread->Draw("CONTZ");
}

  

