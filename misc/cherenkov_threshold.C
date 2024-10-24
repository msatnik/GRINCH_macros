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
#include<math.h>
#include <stack>
#include "TLorentzVector.h"
#include "TCut.h"
#include "TLatex.h"
#include "TLine.h"
#include <TArrayD.h>


// globals
double M_pi_minus = 139.57039;// MeV/C^2
double M_e = 0.51099895069;// MeV/C^2

double n_C4F8 = 1.00132;// 405 nm
double n_C3F8 = 1.00111;// 405 nm
double n_CO2=   1.00045;// 589.29 nm
double n_air =     1.000293;// 589.29 nm

  //// functions
  double calculate_momentum_threshold(double mass, double n);


void cherenkov_threshold(){// main
  
  double p_pi_C4F8 = calculate_momentum_threshold(M_pi_minus, n_C4F8);
  double p_pi_C3F8 = calculate_momentum_threshold(M_pi_minus, n_C3F8);
  double p_pi_CO2 = calculate_momentum_threshold(M_pi_minus, n_CO2);
  double p_pi_air = calculate_momentum_threshold(M_pi_minus, n_air);

  double p_e_C4F8 = calculate_momentum_threshold(M_e, n_C4F8);
  double p_e_C3F8 = calculate_momentum_threshold(M_e, n_C3F8);
  double p_e_CO2 = calculate_momentum_threshold(M_e, n_CO2);
  double p_e_air = calculate_momentum_threshold(M_e, n_air);
  cout<<"------------------------------------------------------------------------------------------"<<endl;
  cout<<"Pion Threshold in C4F8 = "<< p_pi_C4F8<<" MeV/ C^{2}"<<endl;
  cout<<"Pion Threshold in C3F8 = "<< p_pi_C3F8<<" MeV/ C^{2}"<<endl;
  cout<<"Pion Threshold in CO2 = "<< p_pi_CO2<<" MeV/ C^{2}"<<endl;
  cout<<"Pion Threshold in Air = "<< p_pi_air<<" MeV/ C^{2}"<<endl;
  cout<<"------------------------------------------------------------------------------------------"<<endl;
  cout<<"Electron Threshold in C4F8 = "<< p_e_C4F8<<" MeV/ C^{2}"<<endl;
  cout<<"Electron Threshold in C3F8 = "<< p_e_C3F8<<" MeV/ C^{2}"<<endl;
  cout<<"Electron Threshold in CO2 = "<< p_e_CO2<<" MeV/ C^{2}"<<endl;
  cout<<"Electron Threshold in Air = "<< p_e_air<<" MeV/ C^{2}"<<endl;
  cout<<"------------------------------------------------------------------------------------------"<<endl;
}// end main

  double calculate_momentum_threshold(double M, double n)
  {
    double p = M / sqrt(n*n -1);
    return p;
  }
