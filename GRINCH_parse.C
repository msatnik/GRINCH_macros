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


//global params
const Int_t maxTracks = 1000; // Reasonable limit on tracks to be stored per event

void GRINCH_e_p_eff(Int_t entries_input = -1, Int_t kine = 8){
  gStyle->SetNumberContours(255); 
  TChain *C = new TChain("T");

  // these aren't helpful for me at the moment but are a good template for reading in vars 
  Double_t tFitMin = 30; // Minimum number of entries per channel to calibrate ADC/TDC time
  Double_t t_trig = 510; // Mean tdc trig value (HCAL - BB) 
  Double_t E_e = 1.92; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d = 14.5; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th = 35.0; // Angle that the center of HCal is at  
  Double_t W2_mean = 0.93; // Mean of W at current kinematic
  Double_t W2_sig = 0.039; // Width of W at current kinematic
  Double_t GR_cut = 30; //GRINCH ToT sum cutoff 
  Double_t PS_cut = 0.2; //PS energy cutoff



  cout<<"kine: "<<kine<<endl;

  string configfilename = Form("config/sbs%d.cfg",kine);

ifstream configfile(configfilename);
  TString currentline;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      C->Add(currentline);
      cout << "Loaded file at: " << currentline << endl;
    }    
  }
  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }    
    cout<< "Global Cut: "<<globalcut<<endl;
  }
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "tFitMin" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	tFitMin = sval.Atof();
	cout << "Loading timing fit min entries: " << tFitMin << endl;
      }
      if( skey == "t_trig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	t_trig = sval.Atof();
	cout << "Loading mean timing difference BB/HCal trigger: " << t_trig << endl;
      }
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
	cout << "Loading beam energy: " << E_e << endl;
      }
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
      if( skey == "W2_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_mean = sval.Atof();
	cout << "Loading W2 mean cut: " << W2_mean << endl;
      }
      if( skey == "W2_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        W2_sig = sval.Atof();
	cout << "Loading W2 sigma cut: " << W2_sig << endl;
      }
      if( skey == "GR_cut" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        GR_cut = sval.Atof();
	cout << "Loading Grinch ToT cut: " << GR_cut << endl;
      }
      if( skey == "PS_cut" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        PS_cut = sval.Atof();
	cout << "Loading Preshower e cut: " << PS_cut << endl;
      }
    }
    delete tokens;
  }


// Declare parameters
  // HCAL params
  Double_t HCAL_e;
  // BBCAL params 
  Double_t BBps_e, BBsh_e;
  // Track params 
  Double_t BBtr_x[maxTracks], BBtr_y[maxTracks],  BBtr_p[maxTracks], BBtr_th[maxTracks], BBtr_ph[maxTracks];
  Double_t BBtr_n;
  // Hodoscope
  Double_t HODOtmean[1000];
  // Physics
  Double_t kineW2;


  Double_t BBgr_allclus_tmean[1000], BBgr_allclus_adc[1000], BBgr_allclus_size[1000], BBgr_allclus_trms[1000], BBgr_allclus_tot_mean[1000], BBgr_allclus_trackindex[1000], BBgr_allclus_xmean[1000], BBgr_allclus_ymean[1000], BBgr_allclus_dx[1000], BBgr_allclus_dy[1000]; ; 
  Double_t BBgr_clus_tmean, BBgr_clus_adc, BBgr_clus_size, BBgr_clus_trms, BBgr_clus_tot_mean, BBgr_clus_trackindex, BBgr_clus_xmean, BBgr_clus_ymean; 
  Double_t BBgr_hit_amp[1000], BBgr_hit_clustindex[1000], BBgr_hit_col[1000], BBgr_hit_row[1000], BBgr_hit_pmtnum[1000], BBgr_hit_trackindex[1000], BBgr_hit_xhit[1000], BBgr_hit_yhit[1000], BBgr_hit_time[1000];
  Int_t hitsGR; 

  
 Double_t BBgr_bestcluster,  BBgr_ngoodhits, BBgr_ntrackmatch;

  Int_t Ndata_BBgr_tdcelemID;
  Double_t  BBgr_tdc_tdcelemID[1000], BBgr_tdc_tdc[1000],BBgr_tdc_te[1000], BBgr_tdc_mult[1000], BBgr_tdc_tot[1000];


  Double_t BBgr_clus_mirrorindex;
  Double_t BBgr_allclus_mirrorindex[1000];


  // Declare root tree variables and set values to memory locations in root file
  // Turn off all branches
  C->SetBranchStatus( "*", 0 );
  // Turn on specifc branches
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.x", 1 );
  C->SetBranchStatus( "bb.tr.y", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "bb.tr.th", 1 );
  C->SetBranchStatus( "bb.tr.ph", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
 C->SetBranchStatus( "e.kine.W2", 1 );
  C->SetBranchStatus("Ndata.bb.grinch_tdc.allclus.adc",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.adc",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.size",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.t_mean",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.t_rms",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.tot_mean",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.trackindex",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.x_mean",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.y_mean",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.mirrorindex",1);
  C->SetBranchStatus("bb.grinch_tdc.allclus.dx",1);//
  C->SetBranchStatus("bb.grinch_tdc.allclus.dy",1);// C->SetBranchStatus("bb.grinch_tdc.clus.adc",1); // TDC LE sum 
  C->SetBranchStatus("bb.grinch_tdc.clus.size",1); // number of PMTs in the cluster
  C->SetBranchStatus("bb.grinch_tdc.clus.t_mean",1); // average LE time pf the PMTs
  C->SetBranchStatus("bb.grinch_tdc.clus.t_rms",1); // RMS of the average LE time. 
  C->SetBranchStatus("bb.grinch_tdc.clus.tot_mean",1); // mean of the time-over-threshold
  C->SetBranchStatus("bb.grinch_tdc.clus.trackindex",1);// which track the cluster matches (-1 if none)
  C->SetBranchStatus("bb.grinch_tdc.clus.x_mean",1); // the mean x position of the PMTs in the cluster
  C->SetBranchStatus("bb.grinch_tdc.clus.y_mean",1); //  the mean y position of the PMTs in the cluster
  C->SetBranchStatus("bb.grinch_tdc.clus.mirrorindex",1); // which mirror it was matched to  

  C->SetBranchStatus("bb.grinch_tdc.hit.amp",1);
  C->SetBranchStatus("bb.grinch_tdc.hit.clustindex",1);
  C->SetBranchStatus("bb.grinch_tdc.hit.col",1);
  C->SetBranchStatus("bb.grinch_tdc.hit.row",1);
  C->SetBranchStatus("bb.grinch_tdc.hit.pmtnum",1);
  C->SetBranchStatus("bb.grinch_tdc.hit.time",1);
  C->SetBranchStatus("bb.grinch_tdc.hit.trackindex",1);
  C->SetBranchStatus("bb.grinch_tdc.hit.xhit",1);
  C->SetBranchStatus("bb.grinch_tdc.hit.yhit",1);
  C->SetBranchStatus("Ndata.bb.grinch_tdc.hit.pmtnum",1);

  C->SetBranchStatus("bb.grinch_tdc.bestcluster",1);//
  C->SetBranchStatus("bb.grinch_tdc.ngoodhits",1);
  C->SetBranchStatus("bb.grinch_tdc.ntrackmatch",1);

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
  C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );
 C->SetBranchAddress( "e.kine.W2", &kineW2 );
 C->SetBranchAddress("Ndata.bb.grinch_tdc.allclus.adc",&Nclusters);
  C->SetBranchAddress("bb.grinch_tdc.allclus.adc",&BBgr_allclus_adc);
  C->SetBranchAddress("bb.grinch_tdc.allclus.size",&BBgr_allclus_size);
  C->SetBranchAddress("bb.grinch_tdc.allclus.t_mean",&BBgr_allclus_tmean);
  C->SetBranchAddress("bb.grinch_tdc.allclus.t_rms",&BBgr_allclus_trms);
  C->SetBranchAddress("bb.grinch_tdc.allclus.tot_mean",&BBgr_allclus_tot_mean);
  C->SetBranchAddress("bb.grinch_tdc.allclus.trackindex",&BBgr_allclus_trackindex);
  C->SetBranchAddress("bb.grinch_tdc.allclus.x_mean",&BBgr_allclus_xmean);
  C->SetBranchAddress("bb.grinch_tdc.allclus.y_mean",&BBgr_allclus_ymean);
  C->SetBranchAddress("bb.grinch_tdc.allclus.mirrorindex",&BBgr_allclus_mirrorindex);
  C->SetBranchAddress("bb.grinch_tdc.allclus.dx",&BBgr_allclus_dx);//
  C->SetBranchAddress("bb.grinch_tdc.allclus.dy",&BBgr_allclus_dy);//
C->SetBranchAddress("bb.grinch_tdc.clus.adc",&BBgr_clus_adc);
  C->SetBranchAddress("bb.grinch_tdc.clus.size",&BBgr_clus_size);
  C->SetBranchAddress("bb.grinch_tdc.clus.t_mean",&BBgr_clus_tmean);
  C->SetBranchAddress("bb.grinch_tdc.clus.t_rms",&BBgr_clus_trms);
  C->SetBranchAddress("bb.grinch_tdc.clus.tot_mean",&BBgr_clus_tot_mean);
  C->SetBranchAddress("bb.grinch_tdc.clus.trackindex",&BBgr_clus_trackindex);
  C->SetBranchAddress("bb.grinch_tdc.clus.x_mean",&BBgr_clus_xmean);
  C->SetBranchAddress("bb.grinch_tdc.clus.y_mean",&BBgr_clus_ymean);
  C->SetBranchAddress("bb.grinch_tdc.clus.mirrorindex",&BBgr_clus_mirrorindex);

  C->SetBranchAddress("bb.grinch_tdc.hit.amp",&BBgr_hit_amp);
  C->SetBranchAddress("bb.grinch_tdc.hit.clustindex",&BBgr_hit_clustindex);
  C->SetBranchAddress("bb.grinch_tdc.hit.col",&BBgr_hit_col);
  C->SetBranchAddress("bb.grinch_tdc.hit.row",&BBgr_hit_row);
  C->SetBranchAddress("bb.grinch_tdc.hit.pmtnum",&BBgr_hit_pmtnum);
  C->SetBranchAddress("bb.grinch_tdc.hit.time",&BBgr_hit_time);
  C->SetBranchAddress("bb.grinch_tdc.hit.trackindex",&BBgr_hit_trackindex);
  C->SetBranchAddress("bb.grinch_tdc.hit.xhit",&BBgr_hit_xhit);
  C->SetBranchAddress("bb.grinch_tdc.hit.yhit",&BBgr_hit_yhit);
  C->SetBranchAddress("Ndata.bb.grinch_tdc.hit.pmtnum",&hitsGR);

  C->SetBranchAddress("bb.grinch_tdc.bestcluster",&BBgr_bestcluster);//
  C->SetBranchAddress("bb.grinch_tdc.ngoodhits",&BBgr_ngoodhits);
  C->SetBranchAddress("bb.grinch_tdc.ntrackmatch",&BBgr_ntrackmatch);

 // Declare outfile
  //TFile *fout = new TFile( Form("output/sbs%d.root",kine), "RECREATE" ); //need to chage back
  TFile *fout = new TFile("output/sbstest.root", "RECREATE" );
  
  // Histograms
  TH1D* h_W2 =  new TH1D("h_W2",";W2;" ,200, 0, 2);

   TH1D* h_W2_gr_anticut =  new TH1D("h_W2_gr_anticut",";W2 with anticut on grinch;",100,0,2);
  TH1D* h_W2_gr_ps_anticut= new TH1D("h_W2_gr_ps_anticut",";W2 with anticut on grinch and ps;",100,0,2);
  TH1D* h_W2_gr_ps_cut= new TH1D("h_W2_gr_ps_cut",";W2 with cut on grinch and ps;",200,0,2);
  TH1D* h_W2_ps_anticut= new TH1D("h_W2_ps_anticut",";W2 with anticut on ps;",100,0,2);
  TH1D* h_W2_ps_cut= new TH1D("h_W2_ps_cut",";W2 with cut on ps;",200,0,2);
  TH1D* h_W2_gr_cut =  new TH1D("h_W2_gr_cut",";W2 with cut on grinch;",200,0,2);
  TH1D* h_W2_allclus_anticut =  new TH1D("h_W2_allclus_anticut",";W2 with anticut on allclus branch;",100,0,2);
  TH1D* h_W2_allclus_ps_cut =  new TH1D("h_W2_allclus_ps_cut",";W2 with cut on allclus branch and ps;",200,0,2);
  TH1D* h_W2_allclus_ps_anticut =  new TH1D("h_W2_allclus_ps_anticut",";W2 with anticut on allclus branch and ps;",100,0,2);
