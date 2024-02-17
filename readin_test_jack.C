#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include <TStopwatch.h>

//Delcare Global Parameters
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event

void readin_test_jack(Int_t entries_input = -1,Int_t kine = 8){// MAIN

  TChain *C = new TChain("T");

  Double_t HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th = 35.0; // Angle that the center of HCal is at 

  string configfilename = Form("config/test_config%d.cfg",kine);

  cout<<"reading in from a config file"<<endl;

  ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      C->Add(currentline);
      cout << "Loaded file at: " << currentline << endl;
    }    
  }
  // TCut globalcut = "";
  // while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
  //   if( !currentline.BeginsWith("#") ){
  //     globalcut += currentline;
  //   }    
  //   cout<< "Global Cut: "<<globalcut<<endl;
  // }
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
	cout << "Loading HCal distance: " << HCal_d << endl;
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof();	
	cout << "Loading HCal angle: " << HCal_th << endl;
      }
    }
    delete tokens;
  }

  // //Apply Global Cut
  // cout<<endl<<"Populating list with global cut. May take a few minutes...."<<endl;
  // TEventList *elist = new TEventList("elist","Elastic Event List");
  // C->Draw(">>elist",globalcut);
  // cout << endl << "Event list populated with cut: "<<globalcut << endl;

 

  // cout<<endl << "Opened up TChain with nentries: " << C->GetEntries() << ". After globalcut: " << Nevents << "." << endl << endl;
  // cout<<"Entries: "<<Nevents<<endl;


  // DECLARE PARAMETERS
  // HCAL
  Double_t HCAL_e;
  // BBCAL
  Double_t BBps_e, BBsh_e;
 // Track params 
  Double_t BBtr_x[maxTracks], BBtr_y[maxTracks],  BBtr_p[maxTracks], BBtr_th[maxTracks], BBtr_ph[maxTracks];
  Double_t BBtr_n;
  // Physics
  Double_t kineW2;

  // Declare root tree variabls and set values to memeory locations in root file 
  // Turn off all branches
  C->SetBranchStatus( "*", 0 ); 
  // Turn on only branches we need. 
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.x", 1 );
  C->SetBranchStatus( "bb.tr.y", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "bb.tr.th", 1 );
  C->SetBranchStatus( "bb.tr.ph", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );


 // Map branches to the variables 
  C->SetBranchAddress( "sbs.hcal.e", &HCAL_e );
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
  C->SetBranchAddress( "bb.tr.x", &BBtr_x );
  C->SetBranchAddress( "bb.tr.y", &BBtr_y );
  C->SetBranchAddress( "bb.tr.p", &BBtr_p );
  C->SetBranchAddress( "bb.tr.th", &BBtr_th );
  C->SetBranchAddress( "bb.tr.ph", &BBtr_ph );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "e.kine.W2", &kineW2 );
  
  // Declare outfile
  TFile *fout = new TFile( Form("output/sbs_test_jack%d.root",kine), "RECREATE" );


  // Declare Histograms
  TH1D* h_W2 =  new TH1D("h_W2",";W2;" ,200, 0, 2);
  TH1D* h_BBps_e =  new TH1D("h_BBps_e",";BBps_e;" ,200, 0, 2);
  TH1D* h_HCAL_e = new TH1D("h_HCAL_e",";HCAL e;",200,0,3);
  TH1D* h_BBsh_e =  new TH1D("h_BBsh_e",";BBsh_e;" ,200, 0, 3);
  TH2D* h_BBsh_ps_e = new TH2D("h_BBsh_ps_e", "; sh_e  ; ps_e ", 100,0,3,100,0,3 );
  TH1D* h_BB_e_p =  new TH1D("h_BB_e_p",";total shower energy/ p;" ,200, 0, 2);

  //Bool_t globalcut = "BBtr_n==1&&BBps_e>0.2";

   // Set long int to keep track of total entries
  cout<<"Loading branches. Hold tight! "<<endl;
  Long64_t Nevents = C->GetEntries();
  UInt_t run_number = 0;
  cout<<"Entries: "<<Nevents<<endl;

  //check if input for the number of events is valid
  Int_t max = 0;
  if (entries_input == -1){
    max = Nevents;
  }
  else{
    max = entries_input;
  }
  if (max > Nevents){ max = Nevents;}
  cout<<"max = "<<max<<endl; 

  // Loop over events
  for(Long64_t nevent = 0; nevent<max; nevent++){
   
    C->GetEntry(nevent); 
    if (nevent%50000==0) cout << " Entry = " << nevent << endl;

    // Apply global cut 
    if ( BBtr_n!=1) continue; // doing a global cut here for now instead of the elist so things load faster
    

    // General Histos
    h_W2 ->Fill(kineW2);
    h_HCAL_e ->Fill(HCAL_e);
    h_BBps_e ->Fill(BBps_e);
    h_BBsh_e ->Fill(BBsh_e);
    h_BBsh_ps_e ->Fill(BBsh_e,BBps_e);
    Double_t e_over_p  = (BBps_e + BBsh_e)/BBtr_p[0];
    h_BB_e_p -> Fill(e_over_p);
   
  }//end loop over entries

  cout<< "Looped over " << max<< " entries." <<endl;
  fout -> Write();

}
//end main
