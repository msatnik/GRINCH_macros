#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <math.h>

void grtrackmatch(){//MAIN

  //Load Rootfile
  TFile *f1 = TFile::Open("../output/sbs8.root"); // Load rootfile

  //Load Histograms
  TH2D* h_BBgr_clus_xmean_projx= (TH2D*)f1->Get("h_BBgr_clus_xmean_projx");
  TH2D* h_BBgr_clus_xmean_projx_trackmatch= (TH2D*)f1->Get("h_BBgr_clus_xmean_projx_trackmatch");
  TH2D* h_BBgr_allclus_xmean_projx = (TH2D*)f1->Get("h_BBgr_allclus_xmean_projx");
 


}
