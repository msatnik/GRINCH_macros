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


void GRINCH_e_p_eff(Int_t entries_input = -1, Int_t kine = 8){//main
  gStyle->SetNumberContours(255); 
  gStyle->SetCanvasPreferGL(1);
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
  Double_t GR_tw_slope = -0.317; // slope for GRINCH linear timewalk correction 
  Double_t GR_tw_intercept = 5.24; // intecept for GRINCH linear timewalk correction
  Double_t GR_clusx_projx_slope_mirror2=1.354; // slope for clusx (y-axis) vs projx (x-axis) for mirror2
  Double_t GR_clusx_projx_intercept_mirror2 = 0.022; // intercept for clusx (y-axis) vs projx (x-axis) for mirror2
  Double_t GR_clusx_projx_slope_mirror3=1.344; // slope for clusx (y-axis) vs projx (x-axis) for mirror3
  Double_t GR_clusx_projx_intercept_mirror3 = 0.031; // intercept for clusx (y-axis) vs projx (x-axis) for mirror3
  Double_t GR_mirror4_min = 0.6; // cut on projx for mirror3
  Double_t GR_mirror4_max = 0.75;// cut on projx for mirror3
  Double_t GR_mirror3_min = 0.1; // cut on projx for mirror3
  Double_t GR_mirror3_max = 0.45;// cut on projx for mirror3
  Double_t GR_mirror2_min = -0.4; // cut on projx for mirror2
  Double_t GR_mirror2_max = -0.05; // cut on projx for mirror2
  Double_t GR_mirror1_min = -0.57; // cut on projx for mirror3
  Double_t GR_mirror1_max = -0.49;// cut on projx for mirror3

  //// fit parameters for grinch dx vs th. Second try with better mirror cuts based on the projection to the mirror
  const Double_t GR_dx_th_intercept_mirror1 = 0.170;
  const Double_t GR_dx_th_slope_mirror1=1.361;
  const Double_t GR_dx_th_sigma_mirror1=0.013;
  const Double_t GR_dx_th_intercept_mirror2 = 0.005;
  const Double_t GR_dx_th_slope_mirror2 = 1.424;
  const Double_t GR_dx_th_sigma_mirror2=0.012;
  const Double_t GR_dx_th_intercept_mirror3 = 0.0095;
  const Double_t GR_dx_th_slope_mirror3 = 1.339;
  const Double_t GR_dx_th_sigma_mirror3=0.013;
  const Double_t GR_dx_th_intercept_mirror4 = -0.270;
  const Double_t GR_dx_th_slope_mirror4 = 1.466;
  const Double_t GR_dx_th_sigma_mirror4=0.011;
  

  ///// fist try with less good mirror cuts
  // const Double_t GR_dx_th_slope_mirror1= 1.252; 
  // const Double_t GR_dx_th_intercept_mirror1 = 0.143;
  // const Double_t GR_dx_th_slope_mirror2 = 1.432;
  // const Double_t GR_dx_th_intercept_mirror2 = 0.006;
  // const Double_t GR_dx_th_slope_mirror3 = 1.334;
  // const Double_t GR_dx_th_intercept_mirror3 = 0.01;
  // const Double_t GR_dx_th_slope_mirror4 = 1.323;
  // const Double_t GR_dx_th_intercept_mirror4 = -0.244;

  const Double_t SIGMA_dxth = 0.0134;
  const Double_t nSIGMA_dxth = 3;


  // dy ph fit pars with the better mirror cuts 
  const Double_t GR_dy_ph_par0_mirror1 = -0.002;
  const Double_t GR_dy_ph_par1_mirror1 = -1.5;
  const Double_t GR_dy_ph_par2_mirror1 = 10;
  const Double_t GR_dy_ph_par3_mirror1 = -150;
  const Double_t GR_dy_ph_sigma_mirror1 = 0.023;

  const Double_t GR_dy_ph_par0_mirror2 = 0.041;
  const Double_t GR_dy_ph_par1_mirror2 = -1.4;
  const Double_t GR_dy_ph_par2_mirror2 = -5;
  const Double_t GR_dy_ph_par3_mirror2 = -185;
  const Double_t GR_dy_ph_sigma_mirror2 = 0.020;

  const Double_t GR_dy_ph_par0_mirror3 = 0.051;
  const Double_t GR_dy_ph_par1_mirror3 = -1.9;
  const Double_t GR_dy_ph_par2_mirror3 = 0;
  const Double_t GR_dy_ph_par3_mirror3 = -120;
  const Double_t GR_dy_ph_sigma_mirror3 = 0.025;

  const Double_t GR_dy_ph_par0_mirror4 = 0.081;
  const Double_t GR_dy_ph_par1_mirror4 = -1.8;
  const Double_t GR_dy_ph_par2_mirror4 = -25;
  const Double_t GR_dy_ph_par3_mirror4 = 10;
  const Double_t GR_dy_ph_sigma_mirror4 = 0.020;


  ///// first try at dy ph pars 
  // const Double_t GR_dy_ph_par0_mirror1 = -0.004;
  // const Double_t GR_dy_ph_par1_mirror1 = -1.565;
  // const Double_t GR_dy_ph_par2_mirror1 = 10.468;
  // const Double_t GR_dy_ph_par3_mirror1 = -118.935;

  // const Double_t GR_dy_ph_par0_mirror2 = 0.040;
  // const Double_t GR_dy_ph_par1_mirror2 = -1.223;
  // const Double_t GR_dy_ph_par2_mirror2 = -5.231;
  // const Double_t GR_dy_ph_par3_mirror2 = -308.1;

  // const Double_t GR_dy_ph_par0_mirror3 = 0.049;
  // const Double_t GR_dy_ph_par1_mirror3 = -1.708;
  // const Double_t GR_dy_ph_par2_mirror3 = 0.096;
  // const Double_t GR_dy_ph_par3_mirror3 = -283.709;

  // const Double_t GR_dy_ph_par0_mirror4 = 0.084;
  // const Double_t GR_dy_ph_par1_mirror4 = -2.036;
  // const Double_t GR_dy_ph_par2_mirror4 = -31.496;
  // const Double_t GR_dy_ph_par3_mirror4 = 0;

  const Double_t mirror1_mirror2_boundary_location = -0.59;
  const Double_t mirror1_mirror2_boundary_uncert = 0.04;
  const Double_t mirror2_mirror3_boundary_location = 0.01;
  const Double_t mirror2_mirror3_boundary_uncert = 0.04;
  const Double_t mirror3_mirror4_boundary_location = 0.63;
  const Double_t mirror3_mirror4_boundary_uncert = 0.04;

  const Double_t mirror1_max  = mirror1_mirror2_boundary_location - mirror1_mirror2_boundary_uncert;
  const Double_t mirror2_min = mirror1_mirror2_boundary_location + mirror1_mirror2_boundary_uncert;
  const Double_t mirror2_max = mirror2_mirror3_boundary_location - mirror2_mirror3_boundary_uncert;
  const Double_t mirror3_min = mirror2_mirror3_boundary_location + mirror2_mirror3_boundary_uncert;
  const Double_t mirror3_max = mirror3_mirror4_boundary_location -  mirror3_mirror4_boundary_uncert;
  const Double_t mirror4_min = mirror3_mirror4_boundary_location +  mirror3_mirror4_boundary_uncert;

  cout<<"mirror1_max " <<mirror1_max<<endl;
  cout<<"mirror2_min "<< mirror2_min <<endl;
  cout<<"mirror2_max "<< mirror2_max<<endl;
  cout<<"mirror3_min "<< mirror3_min <<endl;
  cout<<"mirror3_max "<< mirror3_max<<endl;
  cout<<"mirror4_min "<< mirror4_min<<endl;

  //SBS9
  //#GR_clusx_projx_slope_mirror2 1.362
  //#GR_clusx_projx_intercept_mirror2 -0.085
  //#GR_clusx_projx_slope_mirror3 1.377
  //#GR_clusx_projx_intercept_mirror3 -0.077

  //SBS8
  // GR_clusx_projx_slope_mirror2 1.354
  // GR_clusx_projx_intercept_mirror2 0.022
  // GR_clusx_projx_slope_mirror3 1.344
  // GR_clusx_projx_intercept_mirror3 0.031



  cout<<"kine: "<<kine<<endl;

  string configfilename = Form("config/sbs%d.cfg",kine);

  // cout<<"testing reading in from a config file"<<endl;

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
      if( skey == "GR_tw_slope" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        GR_tw_slope = sval.Atof();
	cout << "Loading GRINCH time walk slope: " << GR_tw_slope << endl;
      }
      if( skey == "GR_tw_intercept" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        GR_tw_intercept = sval.Atof();
	cout << "Loading GRINCH time walk intercept: " << GR_tw_intercept << endl;
      }
      if( skey == "GR_clusx_projx_slope_mirror2" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	GR_clusx_projx_slope_mirror2= sval.Atof();
	cout << "Loading GRINCH clusterx to projx slope mirror2: " << GR_clusx_projx_slope_mirror2<< endl;
      }
      if( skey == "GR_clusx_projx_intercept_mirror2" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        GR_clusx_projx_intercept_mirror2 = sval.Atof();
	cout << "Loading GRINCH clusterx to projx intercept mirror2: " << GR_clusx_projx_intercept_mirror2<< endl;
      }
      if( skey == "GR_clusx_projx_slope_mirror3" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        GR_clusx_projx_slope_mirror3 = sval.Atof();
	cout << "Loading GRINCH clusterx to projx slope mirror3: " << GR_clusx_projx_slope_mirror3<< endl;
      }
      if( skey == "GR_clusx_projx_intercept_mirror3" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
        GR_clusx_projx_intercept_mirror3 = sval.Atof();
	cout << "Loading GRINCH clusterx to projx intercept mirror3: " << GR_clusx_projx_intercept_mirror3<< endl;
      }
      // if( skey == "GR_mirror1_min" ){
      // 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
      //   GR_mirror1_min = sval.Atof();
      // 	cout << "Loading GRINCH mirror1 min cutoff: " << GR_mirror1_min<< endl;
      // }
      // if( skey == "GR_mirror1_max" ){
      // 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
      // 	GR_mirror1_max = sval.Atof();
      // 	cout << "Loading GRINCH mirror1 max cutoff: " << GR_mirror1_max<< endl;
      // }
      // if( skey == "GR_mirror2_min" ){
      // 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
      //   GR_mirror2_min = sval.Atof();
      // 	cout << "Loading GRINCH mirror2 min cutoff: " << GR_mirror2_min<< endl;
      // }
      // if( skey == "GR_mirror2_max" ){
      // 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
      // 	GR_mirror2_max = sval.Atof();
      // 	cout << "Loading GRINCH mirror2 max cutoff: " << GR_mirror2_max<< endl;
      // }
      // if( skey == "GR_mirror3_min" ){
      // 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
      // 	GR_mirror3_min = sval.Atof();
      // 	cout << "Loading GRINCH mirror3 min cutoff: " << GR_mirror3_min<< endl;
      // }
      // if( skey == "GR_mirror3_max" ){
      // 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
      // 	GR_mirror3_max = sval.Atof();
      // 	cout << "Loading GRINCH mirror3 max cutoff: " << GR_mirror3_max<< endl;
      // }
      // if( skey == "GR_mirror4_min" ){
      // 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
      // 	GR_mirror4_min = sval.Atof();
      // 	cout << "Loading GRINCH mirror4 min cutoff: " << GR_mirror4_min<< endl;
      // }
      // if( skey == "GR_mirror4_max" ){
      // 	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
      // 	GR_mirror4_max = sval.Atof();
      // 	cout << "Loading GRINCH mirror4 max cutoff: " << GR_mirror4_max<< endl;
      // }


    }
    delete tokens;
  }

  //cout<<endl<<"Populating list with global cut. May take a few minutes.... hold tight!"<<endl;
  //TEventList *elist = new TEventList("elist","Elastic Event List");
  //C->Draw(">>elist",globalcut);

  // cout << endl << "Event list populated with cut: "<<globalcut << endl;


  // Declare parameters
  // HCAL params
  Double_t HCAL_e;
  // BBCAL params 
  Double_t BBps_e, BBsh_e;
  // Track params 
  Double_t BBtr_x[maxTracks], BBtr_y[maxTracks],  BBtr_p[maxTracks], BBtr_th[maxTracks], BBtr_ph[maxTracks],BBtr_vz[maxTracks];
  Double_t BBtr_n;
  // Hodoscope
  Double_t HODOtmean[1000];
  // Physics
  Double_t kineW2;
  // GRINCH
  Double_t BBgr_allclus_tmean[1000], BBgr_allclus_adc[1000], BBgr_allclus_size[1000], BBgr_allclus_trms[1000], BBgr_allclus_tot_mean[1000], BBgr_allclus_trackindex[1000], BBgr_allclus_xmean[1000], BBgr_allclus_ymean[1000], BBgr_allclus_dx[1000], BBgr_allclus_dy[1000]; ; 
  Double_t BBgr_clus_tmean, BBgr_clus_adc, BBgr_clus_size, BBgr_clus_trms, BBgr_clus_tot_mean, BBgr_clus_trackindex, BBgr_clus_xmean, BBgr_clus_ymean; 
  Double_t BBgr_hit_amp[1000], BBgr_hit_clustindex[1000], BBgr_hit_col[1000], BBgr_hit_row[1000], BBgr_hit_pmtnum[1000], BBgr_hit_trackindex[1000], BBgr_hit_xhit[1000], BBgr_hit_yhit[1000], BBgr_hit_time[1000];
  Int_t hitsGR; 

  Int_t Ndata_BBgr_tdcelemID;
  Double_t  BBgr_tdc_tdcelemID[1000], BBgr_tdc_tdc[1000],BBgr_tdc_te[1000], BBgr_tdc_mult[1000], BBgr_tdc_tot[1000];


  Double_t BBgr_clus_mirrorindex;
  Double_t BBgr_allclus_mirrorindex[1000];

  Double_t bb_tdctrig_tdc[1000], bb_tdctrig_tdcelemID[1000];
  Int_t Nclusters; 
  Int_t Ndata_bb_tdctrig_tdcelemID;

  Double_t BBgr_bestcluster,  BBgr_ngoodhits, BBgr_ntrackmatch;

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
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.hodotdc.clus.tmean", 1 );
  C->SetBranchStatus("bb.tdctrig.tdc",1);
  C->SetBranchStatus("bb.tdctrig.tdcelemID",1);
  C->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
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
  C->SetBranchStatus("bb.grinch_tdc.allclus.dy",1);//

  // track-matched cluster variables
  C->SetBranchStatus("bb.grinch_tdc.clus.adc",1); // TDC LE sum 
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
  C->SetBranchAddress( "bb.tr.vz", &BBtr_vz );
  C->SetBranchAddress( "bb.tr.ph", &BBtr_ph );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.hodotdc.clus.tmean", &HODOtmean );
  C->SetBranchAddress("bb.tdctrig.tdc",&bb_tdctrig_tdc);
  C->SetBranchAddress("bb.tdctrig.tdcelemID",&bb_tdctrig_tdcelemID);
  C->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID",&Ndata_bb_tdctrig_tdcelemID);
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

  // // the hit closest to the good time cut
  //  C->SetBranchStatus("Ndata.bb.grinch_tdc.tdcelemID",1) ;
  //  C->SetBranchStatus("bb.grinch_tdc.tdcelemID",1) ; // The PMT number 
  //  C->SetBranchStatus("bb.grinch_tdc.tdc",1) ; // leading edge
  //  C->SetBranchStatus("bb.grinch_tdc.tdc_te",1) ; // trailing edge
  //  C->SetBranchStatus("bb.grinch_tdc.tdc_mult",1) ; //tdc multiplicity
  //  C->SetBranchStatus("bb.grinch_tdc.tdc_tot",1) ; // tdc time over threshold (te-le)
  //  C->SetBranchAddress("Ndata.bb.grinch_tdc.tdcelemID",&Ndata_BBgr_tdcelemID) ;
  //  C->SetBranchAddress("bb.grinch_tdc.tdcelemID",&BBgr_tdc_tdcelemID) ; // The PMT number
  //  C->SetBranchAddress("bb.grinch_tdc.tdc",&BBgr_tdc_tdc) ; // leading edge
  //  C->SetBranchAddress("bb.grinch_tdc.tdc_te",&BBgr_tdc_te) ; // trailing edge
  //  C->SetBranchAddress("bb.grinch_tdc.tdc_mult",&BBgr_tdc_mult) ; //tdc multiplicity
  //  C->SetBranchAddress("bb.grinch_tdc.tdc_tot",&BBgr_tdc_tot) ; // tdc time over threshold (te-le)


  Bool_t gr_zeros_cut = BBgr_clus_size>=1&&BBgr_clus_adc>0&&BBgr_clus_trackindex!=-1;
  
  // Declare outfile
  TFile *fout = new TFile( Form("output/sbs%d.root",kine), "RECREATE" ); //need to chage back
  //TFile *fout = new TFile("output/sbstest.root", "RECREATE" );

  // Histograms

  // physics 
  TH1D* h_W2 =  new TH1D("h_W2",";W2;" ,200, 0, 2);
  TH1D* h_W2_gr_anticut =  new TH1D("h_W2_gr_anticut",";W2 with anticut on grinch;",100,0,2);
  TH1D* h_W2_gr_ps_anticut= new TH1D("h_W2_gr_ps_anticut",";W2 with anticut on grinch and ps;",100,0,2);
  TH1D* h_W2_gr_ps_cut= new TH1D("h_W2_gr_ps_cut",";W2 with cut on grinch and ps;",200,0,2);
  TH1D* h_W2_ps_anticut= new TH1D("h_W2_ps_anticut",";W2 with anticut on ps;",100,0,2);
  TH1D* h_W2_ps_cut= new TH1D("h_W2_ps_cut",";W2 with cut on ps;",200,0,2);
  TH1D* h_W2_gr_cut =  new TH1D("h_W2_gr_cut",";W2 with cut on grinch;",200,0,2);
  TH1D* h_W2_elastic =  new TH1D("h_W2_elastic",";W2 with elastic cuts;",200,0,2);
  TH1D* h_W2_elastic_trcut =  new TH1D("h_W2_elastic_trcut",";W2 with elastic cuts and grinch track cut;",200,0,2);
  TH1D* h_W2_allclus_cut =  new TH1D("h_W2_allclus_cut",";W2 with cut on allclus branch;",200,0,2);
  TH1D* h_W2_allclus_anticut =  new TH1D("h_W2_allclus_anticut",";W2 with anticut on allclus branch;",100,0,2);
  TH1D* h_W2_allclus_ps_cut =  new TH1D("h_W2_allclus_ps_cut",";W2 with cut on allclus branch and ps;",200,0,2);
  TH1D* h_W2_allclus_ps_anticut =  new TH1D("h_W2_allclus_ps_anticut",";W2 with anticut on allclus branch and ps;",100,0,2);


  // BBCAL
  TH1D* h_BBps_e =  new TH1D("h_BBps_e",";BBps_e;" ,300, 0, 3);
  TH1D* h_HCAL_e = new TH1D("h_HCAL_e",";HCAL e;",300,0,3);
  TH1D* h_BBsh_e =  new TH1D("h_BBsh_e",";BBsh_e;" ,400, 0, 4);
  TH2D* h_BBsh_ps_e = new TH2D("h_BBsh_ps_e", "; sh_e  ; ps_e ", 300,0,3,300,0,3 );
  TH1D* h_BB_e_p =  new TH1D("h_BB_e_p",";total shower energy/ p;" ,200, 0, 2);
  TH1D* h_BBps_sh_e = new TH1D("h_BBps_sh_e",";sh_e + ps_e ;",200,0,6);
  TH1D* h_BBps_grcut = new TH1D("h_BBps_grcut",";preshower e (cut on grinch clus size) ; ",300,0,3);
  TH1D* h_BBps_grcut_totsum = new TH1D("h_BBps_grcut_totsum",";preshower e (cut on grinch tot sum) ; ",300,0,3);
  TH1D* h_BBps_granticut = new TH1D("h_BBps_granticut",";preshower e (anticut on clus size) ; ",300,0,3);
  TH1D* h_BBps_granticut_totsum = new TH1D("h_BBps_granticut_totsum",";preshower e (anticut on clus totsum) ; ",300,0,3);
  TH1D* h_BBps_nogrtrack = new TH1D("h_BBps_nogrtrack",";preshower e (anticut on gr track match) ; ",300,0,3);

  TH1D* h_BB_e_p_granticut =  new TH1D("h_BB_e_p_granticut",";total shower energy/ p;" ,200, 0, 2);
  TH1D* h_BB_e_p_granticut_totsum =  new TH1D("h_BB_e_p_granticut_totsum",";total shower energy/ p;" ,200, 0, 2);
  TH1D* h_BB_e_p_grcut =  new TH1D("h_BB_e_p_grcut",";total shower energy/ p;" ,200, 0, 2);
  TH1D* h_BB_e_p_grcut_totsum =  new TH1D("h_BB_e_p_grcut_totsum",";total shower energy/ p;" ,200, 0, 2);

  // Track
  TH1D* h_BBtr_x = new TH1D("h_BBtr_x", "bb.tr.x[0]",100,-1,1);
  TH1D* h_BBtr_th =  new TH1D("h_BBtr_th","bb.tr.th[0]", 60,-0.3,0.3);
  TH1D* h_BBtr_y = new TH1D("h_BBtr_y", "bb.tr.y[0]",100,-1,1);
  TH1D* h_BBtr_ph =  new TH1D("h_BBtr_ph","bb.tr.ph[0]", 60,-0.3,0.3);
  TH1D* h_BBtr_p =  new TH1D("h_BBtr_p","bb.tr.p[0]", 500,0,5);



  // RF test
  TH1D *h_bbtrigger =  new TH1D("h_bbtrigger", ";bb trigger;", 250, 335,385);
  TH1D *h_rftime =  new TH1D("h_rftime", ";rftime;", 500, 0,200);
  TH1D *h_rf_bbtrig_diff = new TH1D("h_rf_bbtrig_diff", ";bb trigger - rftime;", 2000,0,500);
  TH1D *h_test = new TH1D("h_test", ";bb trigger - rftime -(gr clus tmean - hodo tmean);", 2000,0,500);
  TH1D *h_rf_bbtrig_gr_diff = new TH1D("h_rf_bbtrig_gr_diff ", "; bb trigger - rftime -(gr clus tmean - hodo tmean);", 2000,0,500); 
  TH1D *h_bbtrigger_diff =  new TH1D("h_bbtrigger_diff", ";bb trigger - (gr clus tmean - hodo tmean);", 250, 335,385);
  TH1D *h_rftime_diff =  new TH1D("h_rftime_diff", ";rftime - (gr clus tmean - hodo tmean);", 500, 0,250);
  TH2D *h_rftime_BBgr_clus_tmean = new TH2D("h_rftime_BBgr_clus_tmean", "; rftime ; gr clus tmean",500,0,250,60,-30,30);
  TH2D *h_rftime_BBgr_clus_tmean_hodosub = new TH2D("h_rftime_BBgr_clus_tmean_hodosub", "; rftime ; gr clus tmean - hodo tmean",500,0,250,60,-30,30);
  TH2D *h_bbtrigger_BBgr_clus_tmean_hodosub = new TH2D("h_bbtrigger_BBgr_clus_tmean_hodosub", "; bb trigger ; gr clus tmean - hodo tmean",250,335,385,60,-30,30);
  TH2D *h_bbtrigger_BBgr_clus_tmean = new TH2D("h_bbtrigger_BBgr_clus_tmean", "; bb trigger; gr clus tmean",250,335,385,60,-30,30);
  TH2D *h_rf_bbtrig_gr_tmean = new TH2D("h_rf_bbtrig_gr_tmean", "; bb trigger - rftime  ; gr clus tmean - hodo tmean",2000,0,500,60,-30,30);

  // GRINCH best cluster variables 
  TH1D* h_BBgr_clus_adc = new TH1D("h_BBgr_clus_adc", ";grinch clus ToT sum;",200,0,200);
  TH1D* h_BBgr_clus_size =  new TH1D("h_BBgr_clus_size", ";bb.grinch_tdc.clus.size;",20,0,20);
  TH1D* h_BBgr_clus_size_pscut =  new TH1D("h_BBgr_clus_size_pscut", ";bb.grinch_tdc.clus.size;",20,0,20);
  TH1D* h_BBgr_clus_size_pscut_notrackmatch =  new TH1D("h_BBgr_clus_size_pscut_notrackmatch", ";bb.grinch_tdc.clus.size;",20,0,20);
  TH1D* h_BBgr_clus_size_psanticut =  new TH1D("h_BBgr_clus_size_psanticut", ";bb.grinch_tdc.clus.size;",20,0,20);
  TH1D* h_BBgr_clus_size_psanticut_notrackmatch =  new TH1D("h_BBgr_clus_size_psanticut_notrackmatch", ";bb.grinch_tdc.clus.size;",20,0,20);
  TH1D* h_BBgr_clus_size_notrackmatch =  new TH1D("h_BBgr_clus_size_notrackmatch", ";bb.grinch_tdc.clus.size;",20,0,20);
  TH1D* h_BBgr_clus_tmean =  new TH1D("h_BBgr_clus_tmean", ";bb.grinch_tdc.clus.t_mean - hodo_tmean;",60,-30,30);
  TH1D* h_BBgr_clus_tmean_nohodo =  new TH1D("h_BBgr_clus_tmean_nohodo", ";bb.grinch_tdc.clus.t_mean;",60,-30,30);
  TH1D* h_BBgr_clus_trms =  new TH1D("h_BBgr_clus_trms", ";bb.grinch_tdc.clus.t_rms;",50,0,5);
  TH1D* h_BBgr_clus_tot_mean =  new TH1D("h_BBgr_clus_tot_mean", ";bb.grinch_tdc.clus.tot_mean;",100,0,50);
  TH2D* h_BBgr_clus_tmean_tot =  new TH2D("h_BBgr_clus_tmean_tot",";clus tot; clus tmean ;",50,0,50,60,-30,30);
  TH2D* h_BBgr_clus_tot_mean_size =  new TH2D("h_BBgr_clus_tot_mean_size", ";cluster size;bb.grinch_tdc.clus.tot_mean; clus;",15,0,15,100,0,50);
  TH1D* h_BBgr_clus_trackindex =  new TH1D("h_BBgr_clus_trackindex", ";bb.grinch_tdc.clus.trackindex;",5,-1,4);
  TH1D* h_BBgr_elas_trackindex =  new TH1D("h_BBgr_elas_trackindex", ";track index w elastic cut;",5,-1,4);
  TH1D* h_BBgr_clus_xmean =  new TH1D("h_BBgr_clus_xmean", ";bb.grinch_tdc.clus.xmean;",200,-2,2);
  TH1D* h_BBgr_clus_ymean =  new TH1D("h_BBgr_clus_ymean", ";bb.grinch_tdc.clus.ymean;",100,-0.25,0.25);
  TH2D* h_BBgr_clus_xmean_projx = new TH2D("h_BBgr_clus_xmean_projx", "; projected x at grinch window from track ;cluster x position",200,-1,1,200,-1,1);
  TH2D* h_BBgr_clus_xmean_projx_trackmatch = new TH2D("h_BBgr_clus_xmean_projx_trackmatch", "; projected x at grinch window from track. trackmatched only ;cluster x position",200,-1,1,200,-1,1);
  TH2D* h_BBgr_clus_mirror2 = new TH2D("h_BBgr_clus_mirror2", "Mirror 2; projected y at grinch window from track ;cluster y position",100,-0.25,0.25,100,-0.25,0.25);
  TH2D* h_BBgr_clus_mirror1 = new TH2D("h_BBgr_clus_mirror1", "Mirror 1; projected y at grinch window from track ;cluster y position",100,-0.25,0.25,100,-0.25,0.25);
  TH2D* h_BBgr_clus_mirror3 = new TH2D("h_BBgr_clus_mirror3", "Mirror 3; projected y at grinch window from track ;cluster y position",100,-0.25,0.25,100,-0.25,0.25);
  TH2D* h_BBgr_clus_mirror4 = new TH2D("h_BBgr_clus_mirror4", "Mirror 4; projected y at grinch window from track ;cluster y position",100,-0.25,0.25,100,-0.25,0.25);
  TH1D* h_BBgr_projx = new TH1D("h_BBgr_projx", "; projected x from track;", 100,-1,1);
  TH1D* h_BBgr_projy = new TH1D("h_BBgr_projy", "; projected y from track;", 100,-0.25,0.25);

  TH1D* h_BBgr_projx_clusx_diff =new TH1D("h_BBgr_projx_clusx_diff",";projx - 0.71*clus_xmean;",400,-0.4,0.4);
  TH1D* h_BBgr_projx_clusx_diff_allclus =new TH1D("h_BBgr_projx_clusx_diff_allclus",";projx - 0.71*clus_xmean;",400,-0.4,0.4);
  TH1D* h_BBgr_projx_clusx_diff_allclus_straight =new TH1D("h_BBgr_projx_clusx_diff_allclus_straight",";projx - clus_xmean;",400,-0.4,0.4);
  TH2D* h_BBtr_p_ph = new TH2D("h_BBtr_p_ph",";bb.tr.ph[0];bb.tr.p[0]",600,-0.3,0.3,500,0,5);
  TH2D* h_BBtr_p_th = new TH2D("h_BBtr_p_th",";bb.tr.th[0];bb.tr.p[0]",600,-0.3,0.3,500,0,5);
  TH2D* h_BBtr_p_th_cut_mirror2 = new TH2D("h_BBtr_p_th_cut_mirror2",";bb.tr.p[0];bb.tr.th[0]",250,0,5,300,-0.3,0.3);
  TH2D* h_BBtr_p_th_cut_mirror3 = new TH2D("h_BBtr_p_th_cut_mirror3",";bb.tr.p[0];bb.tr.th[0]",250,0,5,300,-0.3,0.3);

  TH1D* h_BBgr_projx_clusx_diff_clus =new TH1D("h_BBgr_projx_clusx_diff_clus",";projx - 0.71*clus_xmean;",400,-0.4,0.4);
 
  TH2D* h_BBgr_clusxdiff_trth_allclus = new TH2D("h_BBgr_clusxdiff_trth_allclus","; bb.tr.th[0] ; GRINCH dx", 600,-0.3,0.3,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp = new TH2D("h_BBgr_clusxdiff_trp","; bb.tr.p[0] ; projx - 0.71*clus_xmean", 500,0,5,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp_allclus = new TH2D("h_BBgr_clusxdiff_trp_allclus","; bb.tr.p[0] ; grinch dx", 500,0,5,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trph_allclus = new TH2D("h_BBgr_clusxdiff_trph_allclus","; bb.tr.ph[0] ; grinch dx", 400,-0.2,0.2,1000,-0.4, 0.4);

  TH2D* h_BBgr_clusxdiff_trth_clus = new TH2D("h_BBgr_clusxdiff_trth_clus","; bb.tr.th[0] ; GRINCH dx", 600,-0.3,0.3,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp_clus = new TH2D("h_BBgr_clusxdiff_trp_clus","; bb.tr.p[0] ; grinch dx", 500,0,5,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trph_clus = new TH2D("h_BBgr_clusxdiff_trph_clus","; bb.tr.ph[0] ; grinch dx", 400,-0.2,0.2,1000,-0.4, 0.4);


  TH1D* h_BBgr_projx_clusx_diff_allclus_mirror1 =new TH1D("h_BBgr_projx_clusx_diff_allclus_mirror1","mirror 1 ;GRINCH 'dx' ;",200,-0.4,0.4);
  TH2D* h_BBgr_clusxdiff_trth_allclus_mirror1 = new TH2D("h_BBgr_clusxdiff_trth_allclus_mirror1","mirror1 ; bb.tr.th[0] ; GRINCH 'dx' ", 300,-0.3,0.3,200,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp_mirror1 = new TH2D("h_BBgr_clusxdiff_trp_mirror1","mirror1; bb.tr.p[0] ; GRINCH 'dx' ", 500,0,5,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp_allclus_mirror1 = new TH2D("h_BBgr_clusxdiff_trp_allclus_mirror1","mirror1; bb.tr.p[0] ; GRINCH 'dx' ", 500,0,5,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trph_allclus_mirror1 = new TH2D("h_BBgr_clusxdiff_trph_allclus_mirror1","mirror1; bb.tr.ph[0] ; GRINCH 'dx' ", 400,-0.2,0.2,1000,-0.4, 0.4);
  TH1D* h_BBgr_projx_clusx_diff_allclus_mirror2 =new TH1D("h_BBgr_projx_clusx_diff_allclus_mirror2","mirror2;GRINCH 'dx';",200,-0.4,0.4);
  TH2D* h_BBgr_clusxdiff_trth_allclus_mirror2 = new TH2D("h_BBgr_clusxdiff_trth_allclus_mirror2","mirror2; bb.tr.th[0] ;GRINCH 'dx' ", 300,-0.3,0.3,200,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp_mirror2 = new TH2D("h_BBgr_clusxdiff_trp_mirror2","mirror2 ; bb.tr.p[0] ; GRINCH 'dx' ", 500,0,5,200,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp_allclus_mirror2 = new TH2D("h_BBgr_clusxdiff_trp_allclus_mirror2","mirror2; bb.tr.p[0] ; GRINCH 'dx' ", 500,0,5,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trph_allclus_mirror2 = new TH2D("h_BBgr_clusxdiff_trph_allclus_mirror2","mirror2; bb.tr.ph[0] ; GRINCH 'dx' ", 400,-0.2,0.2,1000,-0.4, 0.4);
  TH1D* h_BBgr_projx_clusx_diff_allclus_mirror3 =new TH1D("h_BBgr_projx_clusx_diff_allclus_mirror3","mirror 3 ;GRINCH 'dx' ;",200,-0.4,0.4);
  TH2D* h_BBgr_clusxdiff_trth_allclus_mirror3 = new TH2D("h_BBgr_clusxdiff_trth_allclus_mirror3","mirror3 ; bb.tr.th[0] ; GRINCH 'dx' ", 300,-0.3,0.3,200,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp_mirror3 = new TH2D("h_BBgr_clusxdiff_trp_mirror3","mirror3; bb.tr.p[0] ; GRINCH 'dx' ", 500,0,5,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp_allclus_mirror3 = new TH2D("h_BBgr_clusxdiff_trp_allclus_mirror3","mirror3; bb.tr.p[0] ; GRINCH 'dx' ", 500,0,5,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trph_allclus_mirror3 = new TH2D("h_BBgr_clusxdiff_trph_allclus_mirror3","mirror3; bb.tr.ph[0] ; GRINCH 'dx' ", 400,-0.2,0.2,1000,-0.4, 0.4);
  TH1D* h_BBgr_projx_clusx_diff_allclus_mirror4 =new TH1D("h_BBgr_projx_clusx_diff_allclus_mirror4","mirror 4 ;GRINCH 'dx' ;",200,-0.4,0.4);
  TH2D* h_BBgr_clusxdiff_trth_allclus_mirror4 = new TH2D("h_BBgr_clusxdiff_trth_allclus_mirror4","mirror4 ; bb.tr.th[0] ; GRINCH 'dx' ", 300,-0.3,0.3,200,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp_mirror4 = new TH2D("h_BBgr_clusxdiff_trp_mirror4","mirror4; bb.tr.p[0] ; GRINCH 'dx' ", 500,0,5,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trp_allclus_mirror4 = new TH2D("h_BBgr_clusxdiff_trp_allclus_mirror4","mirror4; bb.tr.p[0] ; GRINCH 'dx' ", 500,0,5,400,-0.4, 0.4);
  TH2D* h_BBgr_clusxdiff_trph_allclus_mirror4 = new TH2D("h_BBgr_clusxdiff_trph_allclus_mirror4","mirror4; bb.tr.ph[0] ; GRINCH 'dx' ", 400,-0.2,0.2,1000,-0.4, 0.4);

  TH2D* h_BBgr_clusxdiff_mirrornum_clus = new TH2D("h_BBgr_clusxdiff_mirrornum_clus",";mirror number ; GRINCH dx",6,-1,5,200,-0.4,0.4);
  TH2D* h_BBgr_clusxdiff_mirrornum_allclus = new TH2D("h_BBgr_clusxdiff_mirrornum_allclus",";mirror number ; GRINCH dx",6,-1,5,200,-0.4,0.4);

  TH1D* h_track_cut_dx = new TH1D("h_track_cut_dx","all mirrors ;clusx - (projx at window + expected dx);",50,-0.1,0.1);
  TH1D* h_track_cut_dx_mirror1 = new TH1D("h_track_cut_dx_mirror1","mirror1;clusx - (projx at window + expected dx);",50,-0.1,0.1);
  TH1D* h_track_cut_dx_mirror2 = new TH1D("h_track_cut_dx_mirror2","mirror2;clusx - (projx at window + expected dx);",50,-0.1,0.1);
  TH1D* h_track_cut_dx_mirror3 = new TH1D("h_track_cut_dx_mirror3","mirror3;clusx - (projx at window + expected dx);",50,-0.1,0.1);
  TH1D* h_track_cut_dx_mirror4 = new TH1D("h_track_cut_dx_mirror4","mirror4;clusx - (projx at window + expected dx);",50,-0.1,0.1);
  TH1D* h_track_cut_dy = new TH1D("h_track_cut_dy","all mirrors ;clusy - (projy at window + expected dy);",50,-0.1,0.1);
  TH1D* h_track_cut_dy_mirror1 = new TH1D("h_track_cut_dy_mirror1","mirror1;clusy - (projy at window + expected dy);",50,-0.1,0.1);
  TH1D* h_track_cut_dy_mirror2 = new TH1D("h_track_cut_dy_mirror2","mirror2;clusy - (projy at window + expected dy);",50,-0.1,0.1);
  TH1D* h_track_cut_dy_mirror3 = new TH1D("h_track_cut_dy_mirror3","mirror3;clusy - (projy at window + expected dy);",50,-0.1,0.1);
  TH1D* h_track_cut_dy_mirror4 = new TH1D("h_track_cut_dy_mirror4","mirror4;clusy - (projy at window + expected dy);",50,-0.1,0.1);

  TH2D* h_mirror_projection = new TH2D("h_mirror_projection",";projx + 0.66*BBtr_th[0] ;GRINCH dx",200,-1,1,90,-0.3,0.3);

  TH2D* h_BBgr_projx_dx = new TH2D("h_BBgr_projx_dx", "; projected x from track; actual dx", 100,-1,1, 400,-0.4,0.4);
  TH2D* h_BBgr_projx_dx_clus = new TH2D("h_BBgr_projx_dx_clus", "; projected x from track; actual dx", 100,-1,1, 400,-0.4,0.4);

  TH1D* h_BBgr_projy_clusy_diff_allclus = new TH1D("h_BBgr_projy_clusy_diff_allclus",";cluster y - projected y from track ",90,-0.3,0.3);
  TH1D* h_BBgr_projy_clusy_diff_clus = new TH1D("h_BBgr_projy_clusy_diff_clus",";cluster y - projected y from track ",90,-0.3,0.3);
  TH2D* h_BBgr_clusydiff_trph_allclus = new TH2D("h_BBgr_clusydiff_trph_allclus","; bb.tr.ph[0] ; GRINCH 'dy' ", 90,-0.15,0.15,90,-0.3, 0.3);
  TH2D* h_BBgr_clusydiff_trph_clus = new TH2D("h_BBgr_clusydiff_trph_clus","; bb.tr.ph[0] ; GRINCH 'dy' ", 90,-0.15,0.15,90,-0.3, 0.3);
  TH1D* h_BBgr_projy_clusy_diff_allclus_mirror2 = new TH1D("h_BBgr_projy_clusy_diff_allclus_mirror2","mirror2;cluster y - projected y from track ",90,-0.3,0.3);
  TH2D* h_BBgr_clusydiff_trph_allclus_mirror2 = new TH2D("h_BBgr_clusydiff_trph_allclus_mirror2"," mirror2; bb.tr.ph[0] ; GRINCH 'dy' ", 90,-0.15,0.15,90,-0.3, 0.3);
  TH1D* h_BBgr_projy_clusy_diff_allclus_mirror3 = new TH1D("h_BBgr_projy_clusy_diff_allclus_mirror3","mirror3;cluster y - projected y from track ",90,-0.3,0.3);
  TH2D* h_BBgr_clusydiff_trph_allclus_mirror3 = new TH2D("h_BBgr_clusydiff_trph_allclus_mirror3"," mirror3; bb.tr.ph[0] ; GRINCH 'dy' ", 90,-0.15,0.15,90,-0.3, 0.3);
  TH1D* h_BBgr_projy_clusy_diff_allclus_mirror1 = new TH1D("h_BBgr_projy_clusy_diff_allclus_mirror","mirror1;cluster y - projected y from track ",90,-0.3,0.3);
  TH2D* h_BBgr_clusydiff_trph_allclus_mirror1 = new TH2D("h_BBgr_clusydiff_trph_allclus_mirror1"," mirror1; bb.tr.ph[0] ; GRINCH 'dy' ", 90,-0.15,0.15,90,-0.3, 0.3);
  TH1D* h_BBgr_projy_clusy_diff_allclus_mirror4 = new TH1D("h_BBgr_projy_clusy_diff_allclus_mirror4","mirror4;cluster y - projected y from track ",90,-0.3,0.3);
  TH2D* h_BBgr_clusydiff_trph_allclus_mirror4 = new TH2D("h_BBgr_clusydiff_trph_allclus_mirror4"," mirror4; bb.tr.ph[0] ; GRINCH 'dy' ", 90,-0.15,0.15,90,-0.3, 0.3);


  // GRINCH all cluster variables 
  TH1D* h_BBgr_allclus_adc = new TH1D("h_BBgr_allclus_adc", ";bb.grinch_tdc.allclus.adc;",200,0,200);
  TH1D* h_BBgr_allclus_size =  new TH1D("h_BBgr_allclus_size", ";bb.grinch_tdc.allclus.size;",20,0,20);
  TH1D* h_BBgr_allclus_tmean =  new TH1D("h_BBgr_allclus_tmean", ";bb.grinch_tdc.allclus.t_mean;",60,-30,30);
  TH1D* h_BBgr_allclus_trms =  new TH1D("h_BBgr_allclus_trms", ";bb.grinch_tdc.allclus.t_rms;",50,0,5);
  TH1D* h_BBgr_allclus_tot_mean =  new TH1D("h_BBgr_allclus_tot_mean", ";bb.grinch_tdc.allclus.tot_mean;",50,0,5);

  TH1D* h_BBgr_allclus_trackindex =  new TH1D("h_BBgr_allclus_trackindex", ";bb.grinch_tdc.allclus.trackindex;",5,-1,4);
  TH2D* h_BBgr_allclus_size_clusindex =  new TH2D("h_BBgr_allclus_size_clusindex",";cluster index; cluster size", 20,0,20,20,0,20);
  TH1D* h_BBgr_allclus_xmean =  new TH1D("h_BBgr_allclus_xmean", ";bb.grinch_tdc.allclus.xmean;",200,-1,1);
  TH1D* h_BBgr_allclus_ymean =  new TH1D("h_BBgr_allclus_ymean", ";bb.grinch_tdc.allclus.ymean;",100,-0.25,0.25);
  TH1D* h_BBgr_allclus_Nclusters = new TH1D("h_BBgr_allclus_Nclusters", ";Number of clusters in event;", 20,0,20);
  TH2D* h_BBgr_allclus_xmean_projx = new TH2D("h_BBgr_allclus_xmean_projx", "; projected x at grinch window from track ;cluster x position",200,-1,1,200,-1,1);
  TH2D* h_BBgr_allclus_xmean_tmean =  new TH2D("h_BBgr_allclus_xmean_tmean", ";xmean;tmean",200,-1,1,60,-30,30);
  TH2D* h_BBgr_clus_xmean_tmean =  new TH2D("h_BBgr_clus_xmean_tmean", ";xmean;tmean",200,-1,1,60,-30,30);
  TH2D* h_BBgr_allclus_mirror2 = new TH2D("h_BBgr_allclus_mirror2", "Mirror 2; projected y at grinch window from track ;cluster y position",100,-0.25,0.25,100,-0.25,0.25);
  TH2D* h_BBgr_allclus_mirror1 = new TH2D("h_BBgr_allclus_mirror1", "Mirror 1; projected y at grinch window from track ;cluster y position",100,-0.25,0.25,100,-0.25,0.25);
  TH2D* h_BBgr_allclus_mirror3 = new TH2D("h_BBgr_allclus_mirror3", "Mirror 3; projected y at grinch window from track ;cluster y position",100,-0.25,0.25,100,-0.25,0.25);
  TH2D* h_BBgr_allclus_mirror4 = new TH2D("h_BBgr_allclus_mirror4", "Mirror 4; projected y at grinch window from track ;cluster y position",100,-0.25,0.25,100,-0.25,0.25);
 
  TH2D* h_BBgr_clus_size_ps = new TH2D("h_BBgr_clus_size_ps",";preshower e ; cluster size", 200,0,3 ,15,0,15);
  TH2D* h_BBgr_clus_size_ps_trackmatch = new TH2D("h_BBgr_clus_size_ps_trackmatch",";preshower e ; cluster size", 200,0,3 ,15,0,15);
  TH2D* h_BBgr_clus_size_ps_notrackmatch = new TH2D("h_BBgr_clus_size_ps_notrackmatch",";preshower e ; cluster size", 200,0,3 ,15,0,15);
  TH2D* h_BBgr_clus_adc_ps = new TH2D("h_BBgr_clus_adc_ps","; preshower e ; cluster ToT sum", 200,0,3, 200,0,200);
  TH2D* h_BBgr_clus_adc_ps_trackmatch = new TH2D("h_BBgr_clus_adc_ps_trackmatch","; preshower e ; cluster ToT sum", 200,0,3, 200,0,200);
  TH2D* h_BBgr_clus_adc_ps_notrackmatch = new TH2D("h_BBgr_clus_adc_ps_notrackmatch","; preshower e ; cluster ToT sum", 200,0,3, 200,0,200);
  TH2D* h_BBgr_clus_size_adc =  new TH2D("h_BBgr_clus_size_adc",";clus ToT sum; clus size;",200,0,200,15,0,15);
  TH2D* h_BBgr_allclus_adc_ps = new TH2D("h_BBgr_allclus_adc_ps","; preshower e ; cluster ToT sum", 200,0,3, 200,0,200);


  TH1D* h_BBgr_hit_time =  new TH1D("h_BBgr_hit_time", ";bb.grinch_tdc.hit.time -bb.hodotdc.clus.tmean[0] ;",60,-30,30);
  TH1D* h_BBgr_hit_time_nohodo =  new TH1D("h_BBgr_hit_time_nohodo", ";bb.grinch_tdc.hit.time;",60,-30,30);
  TH1D* h_BBgr_hit_time_pscut =  new TH1D("h_BBgr_hit_time_pscut", ";bb.grinch_tdc.hit.time -bb.hodotdc.clus.tmean[0] ;",60,-30,30);
  TH1D* h_BBgr_hit_time_tw_poly = new TH1D("h_BBgr_hit_time_tw_poly", "LE with poly tw corrrection", 60,-30,30);
  TH1D* h_BBgr_hit_time_tw_trad = new TH1D("h_BBgr_hit_time_tw_trad", "LE with trad tw corrrection", 60,-30,30);
  
  TH1D* h_BBgr_hit_numhits = new TH1D("h_BBgr_hit_numhits", ";bb.grinch_tdc.hit.time;",100,0,100);
  TH1D* h_BBgr_hit_amp =  new TH1D("h_BBgr_hit_amp", ";bb.grinch_tdc.hit.amp;",50,0,50);
  TH1D* h_BBgr_hit_clusindex =  new TH1D("h_BBgr_hit_clusindex"," ;bb.grinch_tdc.hit.clusindex;", 21,-1,20);
  TH2D* h_BBgr_hit_clusindex_time =  new TH2D("h_BBgr_hit_clusindex_time"," ;bb.grinch_tdc.hit.clusindex; bb.grinch_tdc.hit.time ", 21,-1,20, 60,-30,30);
  TH2D* h_BBgr_hit_time_amp =  new TH2D("h_BBgr_hit_time_amp","all PMTs with hits;hit branch tot; hit branch time -hodo_tmean ;",50,0,50,60,-30,30);
  TH2D* h_BBgr_hit_time_amp_tw_poly =  new TH2D("h_BBgr_hit_time_amp_tw_poly","all PMTs with hits;hit branch tot; hit branch time -hodo_tmean with poly tw correction;",50,0,50,60,-30,30);
  TH2D* h_BBgr_hit_time_amp_tw_trad=  new TH2D("h_BBgr_hit_time_amp_tw_trad","all PMTs with hits;hit branch tot; hit branch time -hodo_tmean with trad tw correction;",50,0,50,60,-30,30);

  TH2D* h_BBgr_hit_time_amp_nohodo =  new TH2D("h_BBgr_hit_time_amp_nohodo","all PMTs with hits;hit branch tot; hit branch time ;",50,0,50,60,-30,30);
  TH2D* h_BBgr_hit_time_amp_250 =  new TH2D("h_BBgr_hit_time_amp_250","PMT 250;hit branch tot; hit branch time ;",50,0,50,60,-30,30);
  TH2D  *h_BBgr_hit_time_elemID = new TH2D("h_BBgr_hit_time_elemID", "LE vs PMT for each PMT; PMT number; LE - hodo_tmean", 510,0,510, 60,-30,30 );
  TH2D  *h_BBgr_hit_time_elemID_nohodo = new TH2D("h_BBgr_hit_time_elemID_nohodo", "LE vs PMT for each PMT; PMT number; LE", 510,0,510, 60,-30,30 );
  TH2D  *h_BBgr_hit_time_elemID_pscut = new TH2D("h_BBgr_hit_time_elemID_pscut", "LE vs PMT for each PMT; PMT number; LE - hodo_tmean", 510,0,510, 60,-30,30 );

  TH2D  *h_BBgr_hit_time_elemID_bestcluster = new TH2D("h_BBgr_hit_time_elemID_bestcluster", " Best Cluster Only ; PMT number; LE - hodo_tmean", 510,0,510, 60,-30,30 );

  TH1D* h_BBgr_tdc_tdc =  new TH1D("h_BBgr_tdc_tdc","tdc_tdc",60,-30,30);
 
  // Set long int to keep track of total entries
  cout<<"Loading branches. Hold tight! "<<endl;
  Long64_t Nevents = C->GetEntries();
  UInt_t run_number = 0;
 
  cout<<"Entries: "<<Nevents<<endl;
 
  Double_t projx;
  Double_t projy;
  Double_t projx_clusx_diff;
  Double_t projx_clusx_diff_allclus;
  Double_t projx_clusx_diff_clus;
  Double_t projx_clusx_diff_allclus_mirror1;
  Double_t projx_clusx_diff_allclus_mirror2;
  Double_t projx_clusx_diff_allclus_mirror3;
  Double_t projx_clusx_diff_allclus_mirror4;
  Double_t projx_clusx_diff_allclus_straight;
  Double_t dx_expected; 
  Double_t dy_expected;
  Double_t track_cut_dx;
  Double_t track_cut_dy;
  Double_t projy_clusy_diff_allclus;
  Double_t projy_clusy_diff_clus;
  Double_t projy_clusy_diff_allclus_mirror2;
  Double_t projy_clusy_diff_allclus_mirror3;
  Double_t projy_clusy_diff_allclus_mirror1;
  Double_t projy_clusy_diff_allclus_mirror4;
  Double_t mirror_projection;

  
  // Looking at various cuts via 2D arrays of histograms. 
  const int nBins_ps = 6;
  const int nBins_gr_size = 6;
  double ps_cuts[nBins_ps] = {0, 0.16, 0.18, 0.20, 0.22, 0.24};
  int gr_size_cuts[nBins_gr_size] = {0,1,2,3,4,5};
  std::vector<double> gr_tot_cuts = {15,20,25,30,35,40,45,50,55,60};

  /// creating arrays of histograms to look at the various cuts. 
  TH1D* w2_array_cut[nBins_gr_size][nBins_ps]; // [row, col] [gr,ps]
  TH1D* w2_array_anticut[nBins_gr_size][nBins_ps]; // [row, col] [gr,ps]
  TH2D* h_array_BBgr_clus_size_ps_e[nBins_gr_size][nBins_ps]; // [row, col] [gr,ps]
  TH1D*  h_array_BBps_grcut[nBins_gr_size][nBins_ps]; // [row, col] [gr,ps]
  TH1D* h_array_BBgr_pscut[nBins_gr_size][nBins_ps]; // [row, col] [gr,ps]
  TH2D* h_array_anticut_BBgr_clus_size_ps_e[nBins_gr_size][nBins_ps]; // [row, col] [gr,ps]
  TH1D*  h_array_anticut_BBps_grcut[nBins_gr_size][nBins_ps]; // [row, col] [gr,ps]
  TH1D* h_array_anticut_BBgr_pscut[nBins_gr_size][nBins_ps]; // [row, col] [gr,ps]

  // loop over to initialize histograms
  for (int gr_cnt = 0 ; gr_cnt <nBins_gr_size ;gr_cnt++){ // row
    for (int ps_cnt = 0; ps_cnt <nBins_ps; ps_cnt++){ //col 
      // Create unique names for the histograms
      std::string histName1 = Form("h_W2_cut_grsize_%d_ps_e_0p%d",gr_size_cuts[gr_cnt], static_cast<int>(ps_cuts[ps_cnt] * 100));
      std::string histName2 = Form("h_W2_anticut_grsize_%d_ps_e_0p%d",gr_size_cuts[gr_cnt], static_cast<int>(ps_cuts[ps_cnt] * 100));
      std::string histName3 = Form("h_BBgr_clus_size_%d_ps_e_0p%d",gr_size_cuts[gr_cnt], static_cast<int>(ps_cuts[ps_cnt] * 100));
      std::string histName4 = Form("h_BBps_grcut_grsize_%d_ps_e_0p%d",gr_size_cuts[gr_cnt], static_cast<int>(ps_cuts[ps_cnt] * 100));
      std::string histName5 = Form("h_BBgr_pscut_grsize_%d_ps_e_0p%d",gr_size_cuts[gr_cnt], static_cast<int>(ps_cuts[ps_cnt] * 100));
      std::string histName6 = Form("h_anticut_BBgr_clus_size_%d_ps_e_0p%d",gr_size_cuts[gr_cnt], static_cast<int>(ps_cuts[ps_cnt] * 100));
      std::string histName7 = Form("h_anticut_BBps_grcut_grsize_%d_ps_e_0p%d",gr_size_cuts[gr_cnt], static_cast<int>(ps_cuts[ps_cnt] * 100));
      std::string histName8 = Form("h_anticut_BBgr_pscut_grsize_%d_ps_e_0p%d",gr_size_cuts[gr_cnt], static_cast<int>(ps_cuts[ps_cnt] * 100));
      // fill the histograms
      w2_array_cut[gr_cnt][ps_cnt] = new TH1D(histName1.c_str(), histName1.c_str(), 200, 0, 2);      
      w2_array_anticut[gr_cnt][ps_cnt] = new TH1D(histName2.c_str(), histName2.c_str(), 200, 0, 2);      
      h_array_BBgr_clus_size_ps_e[gr_cnt][ps_cnt] =  new TH2D(histName3.c_str()," ;preshower energy; grinch cluster size",200,0,2,40,0,40);
      h_array_BBps_grcut[gr_cnt][ps_cnt] =  new TH1D(histName4.c_str(), histName4.c_str(), 200, 0, 2);
      h_array_BBgr_pscut[gr_cnt][ps_cnt] = new TH1D(histName5.c_str(), histName5.c_str(), 40, 0, 40);
      h_array_anticut_BBgr_clus_size_ps_e[gr_cnt][ps_cnt] =  new TH2D(histName6.c_str()," ;preshower energy; grinch cluster size",200,0,2,40,0,40);
      h_array_anticut_BBps_grcut[gr_cnt][ps_cnt] =  new TH1D(histName7.c_str(), histName7.c_str(), 200, 0, 2);
      h_array_anticut_BBgr_pscut[gr_cnt][ps_cnt] = new TH1D(histName8.c_str(), histName8.c_str(), 40, 0, 40);
    }
  }

  std::vector <double> numerator_no_trackmatching={0,0,0,0,0,0,0,0};
  std::vector <double> denominator_no_trackmatching={0,0,0,0,0,0,0,0};
  std::vector <double> numerator={0,0,0,0,0,0,0,0};
  std::vector <double> denominator={0,0,0,0,0,0,0,0};

  
  
  double GR_adc_cut = 25;

  Int_t max = 0;
  if (entries_input == -1){
    max = Nevents;
  }
  else{
    max = entries_input;
  }
  if (max > Nevents){ max = Nevents;}
  cout<<"max = "<<max<<endl;

  Int_t cout_cnt = 0;
  Int_t cout_cnt2 =0;

  Bool_t in_order = true;
  Double_t tempsize = 100;

  // Loop over events
  for(Long64_t nevent = 0; nevent<max; nevent++){
   
    C->GetEntry(nevent); 
    if (nevent%50000==0) cout << " Entry = " << nevent << endl;
    if ( BBtr_n !=1) continue; // doing a global cut here for now instead of the elist so things load faster
    if(kineW2>2) continue; 
    if(abs(BBtr_vz[0])>0.08) continue;

    projx = BBtr_x[0] +BBtr_th[0]*0.48;//0.48
    projy = BBtr_y[0]+BBtr_ph[0]*0.48;
    projx_clusx_diff = projx - 0.71*BBgr_clus_xmean;

    mirror_projection = projx + 0.66*BBtr_th[0]; // grinch window to mirrors is 66cm. BBtr_x[0] + 1.14*BBtr_th[0]
      
	      
    
      
    //if (projx < -0.7 || projx > 0.8 || projy<-0.21 || projy >0.16) continue;
    h_BBgr_projx ->Fill(projx);
    h_BBgr_projy ->Fill(projy);
    h_BBtr_x ->Fill(BBtr_x[0]);
    h_BBtr_th ->Fill(BBtr_th[0]);
    h_BBtr_y ->Fill(BBtr_y[0]);
    h_BBtr_ph ->Fill(BBtr_ph[0]);
    h_BBtr_p ->Fill(BBtr_p[0]); 
    h_BBtr_p_th ->Fill(BBtr_th[0],BBtr_p[0]);
    h_BBtr_p_ph ->Fill(BBtr_ph[0],BBtr_p[0]);

    if (projx > GR_mirror2_min && projx < GR_mirror2_max)//mirror 2
      {
	h_BBtr_p_th_cut_mirror2 ->Fill(BBtr_p[0],BBtr_th[0]);
      }

    if (projx > GR_mirror3_min && projx < GR_mirror3_max)//mirror 3
      {
	h_BBtr_p_th_cut_mirror3 ->Fill(BBtr_p[0],BBtr_th[0]);
      }

    // if (BBtr_x[0]>0.1 && BBtr_x[0]<0.2)
    // 	 {
    // 	    h_BBtr_p_th_cut ->Fill(BBtr_th[0],BBtr_p[0]);
    // 	 }
      

    // if (projx < -0.3 ||  projx > 0 || projy <-0.15 || projy>0.1) continue; //TESTING LOOKING AT ONE MIRROR. NEED TO COMMENT OUT

    if ( BBgr_clus_size>=1 && BBgr_clus_adc>0)
      {
	// h_BBgr_clus_xmean_projx ->Fill(projx, BBgr_clus_xmean);
	if (BBgr_clus_trackindex != -1)
	  {
	    h_BBgr_clus_xmean_projx_trackmatch ->Fill(projx, BBgr_clus_xmean);
	    // h_BBgr_projx_clusx_diff -> Fill(projx_clusx_diff);
	    //h_BBgr_clusxdiff_trp->Fill(BBtr_p[0],projx_clusx_diff);
	  }
      }



    // Looking at RF time
    Double_t rftime;
    Double_t bbtrigger;
    Double_t diff;
    Double_t difftrig;
    Double_t rf_bbtrig_gr_diff;
    Double_t rf_bbtrig_diff;
    if (BBgr_clus_size !=0 && BBgr_clus_trackindex != -1 && BBgr_clus_adc > GR_cut)
      {
	for (Int_t ihit = 0; ihit< Ndata_bb_tdctrig_tdcelemID ; ihit++)
	  {
	    if (bb_tdctrig_tdcelemID[ihit] ==4 )
	      {
		rftime = bb_tdctrig_tdc[ihit];// 
		h_rftime ->Fill(rftime);
		h_rftime_BBgr_clus_tmean ->Fill(rftime, BBgr_clus_tmean);
		h_rftime_BBgr_clus_tmean_hodosub ->Fill(rftime, BBgr_clus_tmean - HODOtmean[0]);
		if (rftime > 0){
		  diff =  rftime - (BBgr_clus_tmean - HODOtmean[0] ); 
		  h_rftime_diff ->Fill(diff);
		}
	      }
	    if (bb_tdctrig_tdcelemID[ihit] ==5 )
	      {
		bbtrigger = bb_tdctrig_tdc[ihit];
		h_bbtrigger ->Fill(bbtrigger);
		h_bbtrigger_BBgr_clus_tmean ->Fill(bbtrigger, BBgr_clus_tmean);
		h_bbtrigger_BBgr_clus_tmean_hodosub ->Fill(bbtrigger, BBgr_clus_tmean - HODOtmean[0]);
		difftrig =  bbtrigger - (BBgr_clus_tmean - HODOtmean[0] ); 
		h_bbtrigger_diff ->Fill(difftrig);
	      }
	  } 
	if (rftime > 0 && bbtrigger >0)
	  {
	    rf_bbtrig_diff = bbtrigger-rftime;
	    h_rf_bbtrig_diff ->Fill(rf_bbtrig_diff);
	    rf_bbtrig_gr_diff = rf_bbtrig_diff - (BBgr_clus_tmean - HODOtmean[0]);
	    h_rf_bbtrig_gr_diff->Fill(rf_bbtrig_gr_diff );
	    h_rf_bbtrig_gr_tmean ->Fill(rf_bbtrig_diff, BBgr_clus_tmean - HODOtmean[0]);
	  }
      }
     


    // General Histos 
    h_W2 ->Fill(kineW2);
    h_HCAL_e ->Fill(HCAL_e);
    h_BBps_e ->Fill(BBps_e);
    h_BBsh_e ->Fill(BBsh_e);
    h_BBsh_ps_e ->Fill(BBsh_e,BBps_e); 
    h_BBps_sh_e ->Fill(BBsh_e + BBps_e);
    //h_BBps_sh_e ->Fill(BBps_e,BBsh_e); //why is this 2d printing in the terminal but not the projx xmean one lol
    Double_t e_over_p  = (BBps_e + BBsh_e)/BBtr_p[0];
    h_BB_e_p -> Fill(e_over_p);

    // GRINCH histos 
    h_BBgr_allclus_Nclusters ->Fill(Nclusters);

    // h_BBgr_clus_size ->Fill(BBgr_clus_size);  

    h_BBgr_clus_trackindex ->Fill(BBgr_clus_trackindex);  
    if(true) //BBgr_clus_trackindex !=-1 && BBgr_clus_size>0
      {
	projx_clusx_diff_clus = BBgr_clus_xmean - projx;
	h_BBgr_clus_xmean_projx ->Fill(projx, BBgr_clus_xmean);
	h_BBgr_projx_clusx_diff_clus -> Fill(projx_clusx_diff_clus);
	h_BBgr_clusxdiff_trp_clus->Fill(BBtr_p[0],projx_clusx_diff_clus);
	h_BBgr_clusxdiff_trth_clus->Fill(BBtr_th[0],projx_clusx_diff_clus);
	h_BBgr_clusxdiff_trph_clus->Fill(BBtr_ph[0],projx_clusx_diff_clus);
	h_BBgr_clus_xmean_tmean ->Fill(BBgr_clus_xmean,BBgr_clus_tmean);
	h_BBgr_projx_dx_clus ->Fill(projx,projx_clusx_diff_clus);
	h_BBgr_clusxdiff_mirrornum_clus ->Fill(BBgr_clus_mirrorindex, projx_clusx_diff_clus);



	projy_clusy_diff_clus = (BBgr_clus_ymean - projy);
	h_BBgr_projy_clusy_diff_clus ->Fill(projy_clusy_diff_clus);
	h_BBgr_clusydiff_trph_clus->Fill(BBtr_ph[0],projy_clusy_diff_clus);

	h_BBgr_clus_adc ->Fill(BBgr_clus_adc); 
	h_BBgr_clus_size ->Fill(BBgr_clus_size);  
	  
	//h_BBgr_clus_tmean ->Fill(BBgr_clus_tmean - HODOtmean[0]);  
	//h_BBgr_clus_tmean_nohodo ->Fill(BBgr_clus_tmean);  
	h_BBgr_clus_trms ->Fill(BBgr_clus_trms);  
	h_BBgr_clus_tot_mean ->Fill(BBgr_clus_tot_mean);  
	h_BBgr_clus_tot_mean_size ->Fill(BBgr_clus_size, BBgr_clus_tot_mean);  
	h_BBgr_clus_ymean ->Fill(BBgr_clus_ymean);
	h_BBgr_clus_xmean ->Fill(BBgr_clus_xmean);
	h_BBgr_clus_tmean_tot ->Fill(BBgr_clus_tot_mean, BBgr_clus_tmean);

	h_BBgr_clus_size_ps ->Fill( BBps_e, BBgr_clus_size);
	h_BBgr_clus_adc_ps ->Fill( BBps_e, BBgr_clus_adc);
	h_BBgr_clus_size_adc ->Fill(BBgr_clus_adc,BBgr_clus_size);
      }

    if (BBps_e < PS_cut) 
      {
	h_W2_ps_anticut ->Fill(kineW2);
	if(BBgr_clus_trackindex !=-1){
	  h_BBgr_clus_size_psanticut ->Fill(BBgr_clus_size);
	}
	else{
	  h_BBgr_clus_size_psanticut_notrackmatch ->Fill(BBgr_clus_size);
	}
      }
    if ((BBgr_clus_size <=2 && BBgr_clus_trackindex != -1)) // was BBgr_clus_adc < GR_cut || BBgr_clus_trackindex == -1 
      {
	h_W2_gr_anticut ->Fill(kineW2);
	h_BBps_granticut->Fill(BBps_e);
	h_BB_e_p_granticut -> Fill(e_over_p);
      }
    if (BBps_e <PS_cut && (BBgr_clus_adc <GR_cut || BBgr_clus_trackindex == -1 ))
      {
	h_W2_gr_ps_anticut ->Fill(kineW2);
      }
    if (BBps_e >PS_cut && (BBgr_clus_adc >GR_cut && BBgr_clus_trackindex != -1 ))
      {
	h_W2_gr_ps_cut ->Fill(kineW2);
      }
    if(BBgr_clus_size >= 3 && BBgr_clus_trackindex !=-1&& BBps_e >0) // was BBgr_clus_adc > GR_cut
      {
	h_W2_gr_cut ->Fill(kineW2);
	//h_BBgr_clus_tmean ->Fill(BBgr_clus_tmean - HODOtmean[0]);  
	//h_BBgr_clus_tmean_nohodo ->Fill(BBgr_clus_tmean);  
	h_BBps_grcut ->Fill(BBps_e);
	h_BB_e_p_grcut -> Fill(e_over_p);
      }

    if(BBgr_clus_adc >=20 && BBgr_clus_trackindex !=-1&& BBps_e >0)
      {
	h_BBps_grcut_totsum ->Fill(BBps_e);
	h_BB_e_p_grcut_totsum -> Fill(e_over_p);
      }

    if ((BBgr_clus_adc <20 && BBgr_clus_trackindex != -1))
      {
	h_BBps_granticut_totsum ->Fill(BBps_e);
	h_BB_e_p_granticut_totsum -> Fill(e_over_p);
      }

    if(BBgr_clus_trackindex == -1)
      {
	h_BBps_nogrtrack ->Fill(BBps_e);
	h_BBgr_clus_size_notrackmatch ->Fill(BBgr_clus_size);
	h_BBgr_clus_size_ps_notrackmatch->Fill( BBps_e, BBgr_clus_size);
	h_BBgr_clus_adc_ps_notrackmatch->Fill( BBps_e, BBgr_clus_adc);
      }

    if(BBgr_clus_trackindex !=-1)
      {
	h_BBgr_clus_size_ps_trackmatch ->Fill( BBps_e, BBgr_clus_size);
	h_BBgr_clus_adc_ps_trackmatch->Fill( BBps_e, BBgr_clus_adc);
      }

    if (BBps_e >PS_cut)
      {
	h_W2_ps_cut ->Fill(kineW2);
	if(BBgr_clus_trackindex !=-1){
	  h_BBgr_clus_size_pscut ->Fill(BBgr_clus_size);
	}
	else{
	  h_BBgr_clus_size_pscut_notrackmatch ->Fill(BBgr_clus_size);
	}
      }
     
      
    if (BBgr_clus_size > 1)
      {
	h_BBgr_clus_tmean ->Fill(BBgr_clus_tmean - HODOtmean[0]);  
	h_BBgr_clus_tmean_nohodo ->Fill(BBgr_clus_tmean);  
      }
    if (BBgr_clus_size == 1)
      {
	//h_BBgr_clus_tmean ->Fill(BBgr_clus_tmean - HODOtmean[0]);  
	//h_BBgr_clus_tmean_nohodo ->Fill(BBgr_clus_tmean);  
      }
    if (BBgr_clus_size == 0)
      {
	// h_BBgr_clus_adc ->Fill(BBgr_clus_adc);
	// h_BBgr_clus_tmean ->Fill(BBgr_clus_tmean - HODOtmean[0]);  
	//h_BBgr_clus_tmean_nohodo ->Fill(BBgr_clus_tmean);  
      }


    if(BBps_e > PS_cut && BBps_e+BBsh_e >=1.3 && HCAL_e >0.1 && abs(e_over_p -1) <0.2 )
      {
	h_W2_elastic ->Fill(kineW2);

	if(abs(kineW2-0.92) < 0.2)
	  {
	    h_BBgr_elas_trackindex  ->Fill(BBgr_clus_trackindex);
	  }
	if(BBgr_clus_trackindex == 0)
	  {
	    h_W2_elastic_trcut ->Fill(kineW2);
	  }
      
      }




    // Mirror Cuts for best cluster
    if (projx > -0.4 && projx < -0.1) //projx > -0.4 && projx < -0.1
      {
	h_BBgr_clus_mirror2 ->Fill(projy, BBgr_clus_ymean);

      }
    if (projx < -0.55)
      {
	h_BBgr_clus_mirror1 ->Fill(projy, BBgr_clus_ymean);
      }
    if (projx > 0.1 && projx < 0.4)
      {
	h_BBgr_clus_mirror3 ->Fill(projy, BBgr_clus_ymean);
      }
    if (projx > 0.6)
      {
	h_BBgr_clus_mirror4 ->Fill(projy, BBgr_clus_ymean);
      }
      
    Int_t maxadc= 0;
    Int_t maxadcindex = 0;

    tempsize = 100;

    for (Int_t i = 0; i < Nclusters; i ++) //loop over multiple clusters
      {
	h_BBgr_allclus_trackindex ->Fill(BBgr_allclus_trackindex[i]);  
	h_BBgr_allclus_adc ->Fill(BBgr_allclus_adc[i]);
	h_BBgr_allclus_size ->Fill(BBgr_allclus_size[i]);
	h_BBgr_allclus_tmean ->Fill(BBgr_allclus_tmean[i] - HODOtmean[0]);  
	h_BBgr_allclus_trms ->Fill(BBgr_allclus_trms[i]);  
	h_BBgr_allclus_tot_mean ->Fill(BBgr_allclus_tot_mean[i]);  
	h_BBgr_allclus_ymean ->Fill(BBgr_allclus_ymean[i]);
	h_BBgr_allclus_size_clusindex ->Fill(i, BBgr_allclus_size[i]);

	  

	//  h_BBgr_allclus_adc_ps ->Fill( BBps_e, BBgr_allclus_adc[i]);
	if( BBgr_allclus_size[i]>0 && BBgr_allclus_tot_mean[i]>0) // BBgr_allclus_size[i]>2
	  {
	    projx_clusx_diff_allclus = BBgr_allclus_xmean[i] - projx;
	    h_BBgr_allclus_xmean_projx ->Fill(projx, BBgr_allclus_xmean[i]);
	    h_BBgr_allclus_xmean ->Fill(BBgr_allclus_xmean[i]);
	    h_BBgr_projx_clusx_diff_allclus -> Fill(projx_clusx_diff_allclus);
	    h_BBgr_clusxdiff_trp_allclus->Fill(BBtr_p[0],projx_clusx_diff_allclus);
	    h_BBgr_clusxdiff_trth_allclus->Fill(BBtr_th[0],projx_clusx_diff_allclus);
	    h_BBgr_clusxdiff_trph_allclus->Fill(BBtr_ph[0],projx_clusx_diff_allclus);
	    h_BBgr_allclus_xmean_tmean ->Fill(BBgr_allclus_xmean[i],BBgr_allclus_tmean[i]);
	    h_BBgr_projx_dx ->Fill(projx,projx_clusx_diff_allclus);

	    h_BBgr_allclus_adc_ps ->Fill( BBps_e, BBgr_allclus_adc[i]);

	    projy_clusy_diff_allclus = (BBgr_allclus_ymean[i] - projy);
	    h_BBgr_projy_clusy_diff_allclus ->Fill(projy_clusy_diff_allclus);
	    h_BBgr_clusydiff_trph_allclus->Fill(BBtr_ph[0],projy_clusy_diff_allclus);
	      
	    h_BBgr_clusxdiff_mirrornum_allclus ->Fill(BBgr_allclus_mirrorindex[i], projx_clusx_diff_allclus);

	    // // trying to find mirror projections 
	    // if(projy > -0.025 && projy < 0.025 && BBtr_ph[0]> -0.01 && BBtr_ph[0]< 0.01){// look at a small slice in y
	    // 	mirror_projection_test = projx + 0.66*BBtr_th[0]; // grinch window to mirrors is 66cm 
	    // 	h_mirror_projection_test ->Fill(mirror_projection_test, projx_clusx_diff_allclus);
	    // }
	    
	    h_mirror_projection  ->Fill(mirror_projection, projx_clusx_diff_allclus);

	      
	    // Mirror Cuts for every cluster
	    if (mirror_projection > mirror2_min && mirror_projection < mirror2_max  )//mirror 2
	      {
		//projx_clusx_diff_allclus_mirror2 = (BBgr_allclus_xmean[i]- GR_clusx_projx_intercept_mirror2)/GR_clusx_projx_slope_mirror2 - projx; // x = (y-b)/m //probably set to slope 1 and intercept zero in the config file.
		h_BBgr_allclus_mirror2 ->Fill(projy, BBgr_allclus_ymean[i]);
		h_BBgr_projx_clusx_diff_allclus_mirror2 -> Fill(projx_clusx_diff_allclus);
		h_BBgr_clusxdiff_trp_allclus_mirror2->Fill(BBtr_p[0],projx_clusx_diff_allclus);
		h_BBgr_clusxdiff_trth_allclus_mirror2->Fill(BBtr_th[0],projx_clusx_diff_allclus);
		h_BBgr_clusxdiff_trph_allclus_mirror2->Fill(BBtr_ph[0],projx_clusx_diff_allclus);
		h_BBgr_projx_clusx_diff_allclus -> Fill(projx_clusx_diff_allclus);//
		h_BBgr_clusxdiff_trth_allclus->Fill(BBtr_th[0],projx_clusx_diff_allclus);
		 
		dx_expected = GR_dx_th_slope_mirror2*BBtr_th[0] +GR_dx_th_intercept_mirror2; // expected linear fit
		track_cut_dx = BBgr_allclus_xmean[i] - (projx+ dx_expected);
		//// clusx - (projx at window + expected dx)
		h_track_cut_dx->Fill(track_cut_dx);//should be around zero!
		h_track_cut_dx_mirror2 ->Fill(track_cut_dx);//should be around zero!

		// projy_clusy_diff_allclus = (BBgr_allclus_ymean[i] - projy);
		h_BBgr_projy_clusy_diff_allclus_mirror2 ->Fill(projy_clusy_diff_allclus);
		h_BBgr_clusydiff_trph_allclus_mirror2->Fill(BBtr_ph[0],projy_clusy_diff_allclus);

		dy_expected = GR_dy_ph_par0_mirror2 + GR_dy_ph_par1_mirror2*BBtr_ph[0] + GR_dy_ph_par2_mirror2*pow(BBtr_ph[0],2) + GR_dy_ph_par3_mirror2*pow(BBtr_ph[0],3);// par0 + par1*ph + par2*ph^2 +par3*ph^3
		track_cut_dy = BBgr_allclus_ymean[i] - (projy+ dy_expected);
		//// clusy - (projy at window + expected dy)
		h_track_cut_dy ->Fill(track_cut_dy);//should be around zero!
		h_track_cut_dy_mirror2 ->Fill(track_cut_dy);//should be around zero!

	      }
	    if (mirror_projection < mirror1_max )// mirror 1. really only works in sbs9
	      {
		//projx_clusx_diff_allclus_mirror1 = BBgr_allclus_xmean[i] - projx;		  
		h_BBgr_allclus_mirror1 ->Fill(projy, BBgr_allclus_ymean[i]);
		h_BBgr_clusxdiff_trp_allclus_mirror1->Fill(BBtr_p[0],projx_clusx_diff_allclus);
		h_BBgr_clusxdiff_trth_allclus_mirror1->Fill(BBtr_th[0],projx_clusx_diff_allclus);
		h_BBgr_clusxdiff_trph_allclus_mirror1->Fill(BBtr_ph[0],projx_clusx_diff_allclus);
		h_BBgr_projx_clusx_diff_allclus -> Fill(projx_clusx_diff_allclus);//
		h_BBgr_clusxdiff_trth_allclus->Fill(BBtr_th[0],projx_clusx_diff_allclus);

		dx_expected = GR_dx_th_slope_mirror1*BBtr_th[0] +GR_dx_th_intercept_mirror1;
		track_cut_dx = BBgr_allclus_xmean[i] - (projx+ dx_expected);
		//// clusx - (projx at window + expected dx)
		h_track_cut_dx ->Fill(track_cut_dx);//should be around zero!
		h_track_cut_dx_mirror1 ->Fill(track_cut_dx);//should be around zero!

		//projy_clusy_diff_allclus = (BBgr_allclus_ymean[i] - projy);
		h_BBgr_projy_clusy_diff_allclus_mirror1 ->Fill(projy_clusy_diff_allclus);
		h_BBgr_clusydiff_trph_allclus_mirror1->Fill(BBtr_ph[0],projy_clusy_diff_allclus);

		dy_expected = GR_dy_ph_par0_mirror1 + GR_dy_ph_par1_mirror1*BBtr_ph[0] + GR_dy_ph_par2_mirror1*pow(BBtr_ph[0],2) + GR_dy_ph_par3_mirror1*pow(BBtr_ph[0],3);// par0 + par1*ph + par2*ph^2 +par3*ph^3
		track_cut_dy = BBgr_allclus_ymean[i] - (projy+ dy_expected);
		//// clusy - (projy at window + expected dy)
		h_track_cut_dy ->Fill(track_cut_dy);//should be around zero!
		h_track_cut_dy_mirror1 ->Fill(track_cut_dy);//should be around zero!

	      }

	    if (mirror_projection > mirror3_min && mirror_projection < mirror3_max)//mirror 3
	      {
		//projx_clusx_diff_allclus_mirror3 = (BBgr_allclus_xmean[i]- GR_clusx_projx_intercept_mirror3)/GR_clusx_projx_slope_mirror3 - projx; // x = (y-b)/m
		h_BBgr_allclus_mirror3 ->Fill(projy, BBgr_allclus_ymean[i]);
		h_BBgr_projx_clusx_diff_allclus_mirror3 -> Fill(projx_clusx_diff_allclus);
		h_BBgr_clusxdiff_trp_allclus_mirror3->Fill(BBtr_p[0],projx_clusx_diff_allclus);
		h_BBgr_clusxdiff_trth_allclus_mirror3->Fill(BBtr_th[0],projx_clusx_diff_allclus);
		h_BBgr_clusxdiff_trph_allclus_mirror3->Fill(BBtr_ph[0],projx_clusx_diff_allclus);
		h_BBgr_projx_clusx_diff_allclus -> Fill(projx_clusx_diff_allclus);//
		h_BBgr_clusxdiff_trth_allclus->Fill(BBtr_th[0],projx_clusx_diff_allclus);

		dx_expected = GR_dx_th_slope_mirror3*BBtr_th[0] +GR_dx_th_intercept_mirror3;
		track_cut_dx = BBgr_allclus_xmean[i] - (projx + dx_expected);
		//// clusx - (projx at window + expected dx)
		h_track_cut_dx ->Fill(track_cut_dx);//should be around zero!
		h_track_cut_dx_mirror3 ->Fill(track_cut_dx);//should be around zero!

		//projy_clusy_diff_allclus= (BBgr_allclus_ymean[i] - projy);
		h_BBgr_projy_clusy_diff_allclus_mirror3 ->Fill(projy_clusy_diff_allclus);
		h_BBgr_clusydiff_trph_allclus_mirror3->Fill(BBtr_ph[0],projy_clusy_diff_allclus);

		dy_expected = GR_dy_ph_par0_mirror3 + GR_dy_ph_par1_mirror3*BBtr_ph[0] + GR_dy_ph_par2_mirror3*pow(BBtr_ph[0],2) + GR_dy_ph_par3_mirror3*pow(BBtr_ph[0],3);// par0 + par1*ph + par2*ph^2 +par3*ph^3
		track_cut_dy = BBgr_allclus_ymean[i] - (projy+ dy_expected);
		//// clusy - (projy at window + expected dy)
		h_track_cut_dy ->Fill(track_cut_dy);//should be around zero!
		h_track_cut_dy_mirror3 ->Fill(track_cut_dy);//should be around zero!

	      }
	    if (mirror_projection > mirror4_min)//mirror4 //really only works in sbs8
	      {
		//projx_clusx_diff_allclus_mirror4 = BBgr_allclus_xmean[i] - projx;		  
		h_BBgr_allclus_mirror4 ->Fill(projy, BBgr_allclus_ymean[i]);
		h_BBgr_clusxdiff_trp_allclus_mirror4->Fill(BBtr_p[0],projx_clusx_diff_allclus);
		h_BBgr_clusxdiff_trth_allclus_mirror4->Fill(BBtr_th[0],projx_clusx_diff_allclus);
		h_BBgr_clusxdiff_trph_allclus_mirror4->Fill(BBtr_ph[0],projx_clusx_diff_allclus);
		h_BBgr_projx_clusx_diff_allclus -> Fill(projx_clusx_diff_allclus);//
		h_BBgr_clusxdiff_trth_allclus->Fill(BBtr_th[0],projx_clusx_diff_allclus);

		dx_expected = GR_dx_th_slope_mirror4*BBtr_th[0] +GR_dx_th_intercept_mirror4;
		track_cut_dx = BBgr_allclus_xmean[i] - (projx+ dx_expected);
		//// clusx - (projx at window + expected dx)
		h_track_cut_dx ->Fill(track_cut_dx);//should be around zero!
		h_track_cut_dx_mirror4 ->Fill(track_cut_dx);//should be around zero!

		// projy_clusy_diff_allclus = (BBgr_allclus_ymean[i] - projy);
		h_BBgr_projy_clusy_diff_allclus_mirror4 ->Fill(projy_clusy_diff_allclus);
		h_BBgr_clusydiff_trph_allclus_mirror4->Fill(BBtr_ph[0],projy_clusy_diff_allclus);

		dy_expected = GR_dy_ph_par0_mirror4 + GR_dy_ph_par1_mirror4*BBtr_ph[0] + GR_dy_ph_par2_mirror4*pow(BBtr_ph[0],2) + GR_dy_ph_par3_mirror4*pow(BBtr_ph[0],3);// par0 + par1*ph + par2*ph^2 +par3*ph^3
		track_cut_dy = BBgr_allclus_ymean[i] - (projy+ dy_expected);
		//// clusy - (projy at window + expected dy)
		h_track_cut_dy ->Fill(track_cut_dy);//should be around zero!
		h_track_cut_dy_mirror4 ->Fill(track_cut_dy);//should be around zero!
	      }


	  }
	  
	if(BBgr_allclus_adc[i] > maxadc)
	  {
	    maxadc = BBgr_allclus_adc[i];
	    maxadcindex = i;
	  }

      }// end loop over Nclusters

    if (BBgr_allclus_adc[maxadcindex] > GR_cut)
      {
	h_W2_allclus_cut ->Fill(kineW2);
      }
    else
      {
	h_W2_allclus_anticut ->Fill(kineW2);
      }

    if(Nclusters ==0)
      {
	h_W2_allclus_anticut ->Fill(kineW2);
      }

 
    if (BBgr_allclus_adc[maxadcindex] > GR_cut && BBps_e > PS_cut)
      {
	h_W2_allclus_ps_cut ->Fill(kineW2);
      }
    else if (BBgr_allclus_adc[maxadcindex]< GR_cut && BBps_e < PS_cut)
      {
	h_W2_allclus_ps_anticut ->Fill(kineW2);
      }
    else if(Nclusters ==0 && BBps_e < PS_cut)
      {
	h_W2_allclus_ps_anticut ->Fill(kineW2);
      }


  

  
    // // testing a simple eff
    // // [0] is the electron detection eff on a cluster size cut
    // // [1] is the electron detection eff on an adc cut
    // // [2] is the pion rejection eff on a size cut
    // // [3] is the pion rejection eff on an adc cut
    // // [4] is the pion rejection eff with a tight pion cut and cluster size cut
    // // [5] is the pion rejection eff with a tight pion cut and adc cut

    // // NOT CHECKING TRACK MATHCING 
    // // [0] is the electron detection eff on a cluster size cut 
    // // [1] is the electron detection eff on an adc cut
    // // [2] is the pion rejection eff on a size cut
    // // [3] is the pion rejection eff on an adc cut
    // // [4] is the pion rejection eff on a tight pion cut and cluster size cut
    // // [5] is the pion rejection eff on a tight pion cut and an adc cut


    bool tight_electron_cut = (BBps_e >0.45 && BBps_e < 0.55 && BBsh_e >0.9 && BBsh_e < 1.05 && abs(e_over_p-0.98)<0.16);

    bool tight_pion_cut = (BBps_e >= 0.08 && BBps_e <= 0.12 && BBsh_e >= 1.1 && BBsh_e <=1.3 && e_over_p<0.8);
    
    
    // [0] is the electron detection eff on a cluster size cut
    if(BBps_e > 0.2){
      denominator[0] = denominator[0] + 1;
      if (BBgr_clus_size >= 3 && BBgr_clus_trackindex !=-1)
	{
	  numerator[0] = numerator[0] +1;
	}
    }


    // [1] is the electron detection eff on an adc cut
    if(BBps_e > 0.2){
      denominator[1] = denominator[1] + 1;
      if (BBgr_clus_adc >  GR_adc_cut && BBgr_clus_trackindex !=-1)
	{
	  numerator[1] = numerator[1] +1;
	}
    }

 // [2] is the electron detection eff on a cluster size cut
    if(tight_electron_cut){
      denominator[2] = denominator[2] + 1;
      if (BBgr_clus_size >= 3 && BBgr_clus_trackindex !=-1)
	{
	  numerator[2] = numerator[2] +1;
	}
    }



    // [1] is the electron detection eff on an adc cut
    if(tight_electron_cut){
      denominator[3] = denominator[3] + 1;
      if (BBgr_clus_adc >  GR_adc_cut && BBgr_clus_trackindex !=-1)
	{
	  numerator[3] = numerator[3] +1;
	}
    }
    
      
  
    // [2] is the pion rejection eff on a size cut
    if (BBps_e <=0.2){
      denominator[4] = denominator[4] + 1;
      if (BBgr_clus_size >= 3 && BBgr_clus_trackindex !=-1)
	{
	  numerator[4] = numerator[4] +1;
	}
    }
      
    // denominator.push_back(0);
    // numerator.push_back(0);
      
    // [3] is the pion rejection eff on an adc cut
    if (BBps_e <=0.2){
      denominator[5] = denominator[5] + 1;
      if (BBgr_clus_adc >GR_adc_cut && BBgr_clus_trackindex !=-1)
	{
	  numerator[5] = numerator[5] +1;
	}
    }

    // [4] is the pion rejection eff on a cluster size cut
    // Tight "pion" cut 
    if (tight_pion_cut )
      {
	denominator[6] = denominator[6] + 1;
	if (BBgr_clus_size >= 3 && BBgr_clus_trackindex !=-1)
	  {
	    numerator[6] = numerator[6] +1;
	  }
      }

    // [5] is the pion rejection eff on an adc cut
    // Tight "pion" cut 
    if (tight_pion_cut)
      {
	denominator[7] = denominator[7] + 1;
	if (BBgr_clus_adc > GR_adc_cut && BBgr_clus_trackindex !=-1)
	  {
	    numerator[7] = numerator[7] +1;
	  }
      }

  

    // NOT checking for track match 
    // [0] is the electron detection eff on a cluster size cut
    if(BBps_e > 0.2){
      denominator_no_trackmatching[0] = denominator_no_trackmatching[0] + 1;
      if (BBgr_clus_size >= 3)
	{
	  numerator_no_trackmatching[0] = numerator_no_trackmatching[0] +1;
	}
    }

    // denominator.push_back(0);
    // numerator.push_back(0);

    // NOT checking for track match 
    // [1] is the electron detection eff on an adc cut
    if(BBps_e > 0.2){
      denominator_no_trackmatching[1] = denominator_no_trackmatching[1] + 1;
      if (BBgr_clus_adc >GR_adc_cut )
	{
	  numerator_no_trackmatching[1] = numerator_no_trackmatching[1] +1;
	}
    }

    

    // [2] is the electron detection eff on a cluster size cut
    if(tight_electron_cut){
      denominator_no_trackmatching[2] = denominator_no_trackmatching[2] + 1;
      if (BBgr_clus_size >= 3 )
	{
	  numerator_no_trackmatching[2] = numerator_no_trackmatching[2] +1;
	}
    }



    // [1] is the electron detection eff on an adc cut
    if(tight_electron_cut){
      denominator_no_trackmatching[3] = denominator_no_trackmatching[3] + 1;
      if (BBgr_clus_adc >  GR_adc_cut )
	{
	  numerator_no_trackmatching[3] = numerator_no_trackmatching[3] +1;
	}
    }
    

    
    // denominator.push_back(0);
    // numerator.push_back(0);
    
    // NOT checking for track match 
    // [2] is the pion rejection eff on a size cut
    if (BBps_e <=0.2){
      denominator_no_trackmatching[4] = denominator_no_trackmatching[4] + 1;
      if (BBgr_clus_size >= 3 )
	{
	  numerator_no_trackmatching[4] = numerator_no_trackmatching[4] +1;
	}
    }
      
    // denominator.push_back(0);
    // numerator.push_back(0);

    // NOT checking for track match 
    // [3] is the pion rejection eff on an adc cut
    if (BBps_e <=0.2){
      denominator_no_trackmatching[5] = denominator_no_trackmatching[5] + 1;
      if (BBgr_clus_adc >GR_adc_cut )
	{
	  numerator_no_trackmatching[5] = numerator_no_trackmatching[5] +1;
	}
    }


    // Tight "pion" cut
    // NOT checking for track match 
    // [4] is the pion rejection eff on a size cut
    if (tight_pion_cut)
      {
	denominator_no_trackmatching[6] = denominator_no_trackmatching[6] + 1;
	if (BBgr_clus_size >= 3)
	  {
	    numerator_no_trackmatching[6] = numerator_no_trackmatching[6] +1;
	  }
      }

    // Tight "pion" cut
    // NOT checking for track match 
    // [5] is the pion rejection eff on an adc cut
    if (tight_pion_cut)
      {
	denominator_no_trackmatching[7] = denominator_no_trackmatching[7] + 1;
	if (BBgr_clus_adc > GR_adc_cut)
	  {
	    numerator_no_trackmatching[7] = numerator_no_trackmatching[7] +1;
	  }
      }



  

    // trying out the 2d array of cuts on preshower and grinch
    for (int gr_cnt = 0 ; gr_cnt <nBins_gr_size ;gr_cnt++){ // row
      for (int ps_cnt = 0; ps_cnt <nBins_ps; ps_cnt++){ //col 
	if (BBps_e >= ps_cuts[ps_cnt] && BBgr_clus_size >= gr_size_cuts[gr_cnt] && BBgr_clus_trackindex !=-1)
	  {
	    w2_array_cut[gr_cnt][ps_cnt] ->Fill(kineW2);
	    h_array_BBgr_clus_size_ps_e[gr_cnt][ps_cnt] ->Fill(BBps_e, BBgr_clus_size);
	    h_array_BBgr_pscut[gr_cnt][ps_cnt] ->Fill(BBgr_clus_size);
	    h_array_BBps_grcut[gr_cnt][ps_cnt]->Fill(BBps_e);
	  }
	else if (BBgr_clus_trackindex !=-1){// need to think about track matching and how that plays together
	  w2_array_anticut[gr_cnt][ps_cnt] ->Fill(kineW2);
	  h_array_anticut_BBgr_clus_size_ps_e[gr_cnt][ps_cnt] ->Fill(BBps_e, BBgr_clus_size);
	  h_array_anticut_BBgr_pscut[gr_cnt][ps_cnt] ->Fill(BBgr_clus_size);
	  h_array_anticut_BBps_grcut[gr_cnt][ps_cnt]->Fill(BBps_e);

	}

      }// end loop over gr_cnt (row)
    }//end loop over ps_cnt (col)

 
    //SBS 9
    //Double_t tw_poly = -0.326;
    //Double_t tw_poly_intercept = 5.29;
    //SBS 8 
    Double_t tw_poly = -0.317;
    Double_t tw_poly_intercept = 5.24;

    Double_t tw_trad = 14.41;
    Double_t tw_trad_intercept = -2.78;



    // loop over the grinch hit brances
    h_BBgr_hit_numhits ->Fill(hitsGR);
    for (Int_t i = 0; i < hitsGR; i++)
      {
	if (BBgr_hit_pmtnum[i] >= 510 || BBgr_hit_amp[i]==0 || BBgr_hit_pmtnum[i]==133) continue;
	Double_t time_hodocor = BBgr_hit_time[i]  - HODOtmean[0];
	Double_t time_poly = time_hodocor - BBgr_hit_amp[i]*tw_poly - tw_poly_intercept;
	Double_t time_trad = time_hodocor - tw_trad * pow(BBgr_hit_amp[i],-0.5) - tw_trad_intercept;
	h_BBgr_hit_time ->Fill(time_hodocor);
	h_BBgr_hit_time_tw_poly ->Fill(time_poly);
	h_BBgr_hit_time_tw_trad ->Fill(time_trad);
	h_BBgr_hit_time_nohodo ->Fill(BBgr_hit_time[i]);
	h_BBgr_hit_amp ->Fill(BBgr_hit_amp[i]);
	h_BBgr_hit_time_amp ->Fill(BBgr_hit_amp[i],time_hodocor);
	h_BBgr_hit_time_amp_tw_poly ->Fill(BBgr_hit_amp[i],time_poly );
	h_BBgr_hit_time_amp_tw_trad ->Fill(BBgr_hit_amp[i],time_trad) ;
	h_BBgr_hit_time_amp_nohodo ->Fill(BBgr_hit_amp[i],BBgr_hit_time[i] );
	h_BBgr_hit_time_elemID ->Fill(BBgr_hit_pmtnum[i],BBgr_hit_time[i] - HODOtmean[0] );
	h_BBgr_hit_time_elemID_nohodo ->Fill(BBgr_hit_pmtnum[i],BBgr_hit_time[i] );	  
	h_BBgr_hit_clusindex ->Fill(BBgr_hit_clustindex[i]);
	h_BBgr_hit_clusindex_time ->Fill(BBgr_hit_clustindex[i], BBgr_hit_time[i]);
	// h_BBgr_hit_clusindex_size ->Fill(BBgr_hit_clustindex[i], BBgr_hit_time[i]);
	if (BBps_e > 0.2)
	  {
	    h_BBgr_hit_time_pscut ->Fill(time_hodocor);
	    h_BBgr_hit_time_elemID_pscut ->Fill(BBgr_hit_pmtnum[i],BBgr_hit_time[i] - HODOtmean[0] );
	  }

       
	if(BBgr_hit_clustindex[i] == BBgr_bestcluster)
	  {
	    h_BBgr_hit_time_elemID_bestcluster -> Fill(BBgr_hit_pmtnum[i],BBgr_hit_time[i] - HODOtmean[0] );
	  }


	if (BBgr_hit_pmtnum[i] == 250)
	  {
	    h_BBgr_hit_time_amp_250 ->Fill(BBgr_hit_amp[i], BBgr_hit_time[i]);
	  }
	if (cout_cnt < 50)
	  {
	    // cout<< "grinch hit time: "<<BBgr_hit_time[i]<<" hodo time: "<< HODOtmean[0]<<endl;
	    cout_cnt++;
	  }
      }//end loop over hits
      
      


    // for (Int_t i = 0 ; i < Ndata_BBgr_tdcelemID; i++)
    // 	{
    // 	  if (BBgr_tdc_tdcelemID[i]<510)
    // 	    {
    // 	      h_BBgr_tdc_tdc ->Fill(BBgr_tdc_tdc[i]);
    // 	      if (cout_cnt2 <50)
    // 		{
    // 		  cout<<"tdc_tdc " << BBgr_tdc_tdc[i] <<endl;
    // 		  cout_cnt2 ++;
    // 		}
    // 	    }
    // 	}
      
   
    
  }//end for loop over Nevents

  cout<<"----------------------------------------------------"<<endl;
  cout<<"Utilizing track-matching: "<<endl;
  cout<<"Electron Detection Eff: clus size >= 3"<<endl;
  cout<<"           = "<<numerator[0]<<" / "<<denominator[0] <<" = " <<numerator[0]/denominator[0]<<endl;

  cout<<"Electron Detection Eff: clus adc > " <<GR_adc_cut<<endl;
  cout<<"           = "<<numerator[1]<<" / "<<denominator[1] <<" = " <<  numerator[1]/denominator[1]<<endl;

  cout<<"Electron Detection Eff with tight electron cut: clus size >= 3"<<endl;
  cout<<"           = "<<numerator[2]<<" / "<<denominator[2] <<" = " <<numerator[2]/denominator[2]<<endl;

  cout<<"Electron Detection Eff with tight electron cut: clus adc > " <<GR_adc_cut<<endl;
  cout<<"           = "<<numerator[3]<<" / "<<denominator[3] <<" = " <<  numerator[3]/denominator[3]<<endl;

  
  cout<<"Pion rejection eff: clus size >= 3"<<endl;
  cout<<"           = "<<denominator[4]-numerator[4]<<" / "<<denominator[4] <<" = " <<(denominator[4]-numerator[4])/denominator[4]<<endl;

  cout<<"Pion rejection eff: clus adc > "<<GR_adc_cut<<endl;
  cout<<"           = "<<denominator[5]-numerator[5]<<" / "<<denominator[5] <<" = " << (denominator[5]- numerator[5])/denominator[5]<<endl;

  cout<<"Pion rejection eff with tight pion spot cut: clus size >= 3"<<endl;
  cout<<"           = "<<denominator[6]-numerator[6]<<" / "<<denominator[6] <<" = " <<(denominator[6]-numerator[6])/denominator[6]<<endl;

  cout<<"Pion rejection eff with tight pion spot cut: clus adc > "<<GR_adc_cut<<endl;
  cout<<"           = "<<denominator[7]-numerator[7]<<" / "<<denominator[7] <<" = " <<  (denominator[7]-numerator[7])/denominator[7]<<endl;
  cout<<"----------------------------------------------------"<<endl;

    
  cout<<"----------------------------------------------------"<<endl;
  cout<<"Not Explicitly Utilizing track-matching: "<<endl;
  cout<<"Electron Detection Eff: clus size >= 3"<<endl;
  cout<<"           = "<<numerator_no_trackmatching[0]<<" / "<<denominator_no_trackmatching[0] <<" = " <<numerator_no_trackmatching[0]/denominator_no_trackmatching[0]<<endl;

  cout<<"Electron Detection Eff: clus adc > 20"<<endl;
  cout<<"           = "<<numerator_no_trackmatching[1]<<" / "<<denominator_no_trackmatching[1] <<" = " <<  numerator_no_trackmatching[1]/denominator_no_trackmatching[1]<<endl;

  cout<<"Electron Detection Eff with tight electron cut: clus size >= 3"<<endl;
  cout<<"           = "<<numerator_no_trackmatching[2]<<" / "<<denominator_no_trackmatching[2] <<" = " <<numerator_no_trackmatching[2]/denominator_no_trackmatching[2]<<endl;

  cout<<"Electron Detection Eff with tight electron cut: clus adc > " <<GR_adc_cut<<endl;
  cout<<"           = "<<numerator_no_trackmatching[3]<<" / "<<denominator_no_trackmatching[3] <<" = " <<  numerator_no_trackmatching[3]/denominator_no_trackmatching[3]<<endl;

  cout<<"Pion rejection eff: clus size >= 3"<<endl;
  cout<<"           = "<<denominator_no_trackmatching[4]-numerator_no_trackmatching[4]<<" / "<<denominator_no_trackmatching[4] <<" = " <<(denominator_no_trackmatching[4]-numerator_no_trackmatching[4])/denominator_no_trackmatching[4]<<endl;

  cout<<"Pion rejection eff: clus adc > " <<GR_adc_cut<<endl;
  cout<<"           = "<<denominator_no_trackmatching[5]-numerator_no_trackmatching[5]<<" / "<<denominator_no_trackmatching[5] <<" = " <<  (denominator_no_trackmatching[5]-numerator_no_trackmatching[5])/denominator_no_trackmatching[5]<<endl;


  cout<<"Pion rejection eff with tight pion spot cut: clus size >= 3"<<endl;
  cout<<"           = "<<denominator_no_trackmatching[6]-numerator_no_trackmatching[6]<<" / "<<denominator_no_trackmatching[6] <<" = " <<(denominator_no_trackmatching[6]-numerator_no_trackmatching[6])/denominator_no_trackmatching[6]<<endl;

  cout<<"Pion rejection eff with tight pion spot cut: clus adc > " <<GR_adc_cut<<endl;
  cout<<"           = "<<denominator_no_trackmatching[7]-numerator_no_trackmatching[7]<<" / "<<denominator_no_trackmatching[7] <<" = " <<  (denominator_no_trackmatching[7]-numerator_no_trackmatching[7])/denominator_no_trackmatching[7]<<endl;
  cout<<"----------------------------------------------------"<<endl;

    

  TCanvas *c1 = new TCanvas("c1","c1",1000,500);

  THStack *stack = new THStack("stack", "preshower energy");
  h_BBps_nogrtrack ->SetLineColorAlpha(kRed,0.35);
  h_BBps_nogrtrack ->SetFillColorAlpha(kRed,0.35);
  h_BBps_granticut ->SetLineColorAlpha(kBlue,0.35);
  h_BBps_granticut ->SetFillColorAlpha(kBlue,0.35);
  h_BBps_grcut ->SetLineColorAlpha(kGreen+2,0.35);
  h_BBps_grcut ->SetFillColorAlpha(kGreen+2,0.35);
  h_BBps_e ->SetLineColorAlpha(kBlack,1);
  stack->Add(h_BBps_nogrtrack);
  stack ->Add(h_BBps_e);
  stack->Add(h_BBps_grcut);
  stack->Add(h_BBps_granticut);
  stack ->Draw("nostack hist");

  //Add a legend
  auto legend = new TLegend();
  legend->AddEntry(h_BBps_e,"bb.ps.e","l");
  legend->AddEntry(h_BBps_grcut,"gr clus size >=3 && track-matched","f");
  legend->AddEntry(h_BBps_granticut,"gr clus size <=2 && track-matched","f");
  legend->AddEntry(h_BBps_nogrtrack,"no track match","f");
  legend->Draw();


  TCanvas *c2 = new TCanvas("c2","c2",1000,500);
  c2->Divide(2,1);
  c2->cd(1);

  THStack *stack2 = new THStack("stack2", "GRINCH cluster size: ps.e>0.2");
  h_BBgr_clus_size ->SetLineColor(kBlack);
  h_BBgr_clus_size_pscut ->SetLineColorAlpha(kGreen+1,0.35);
  h_BBgr_clus_size_pscut ->SetFillColorAlpha(kGreen+1,0.35);
  h_BBgr_clus_size_pscut_notrackmatch ->SetLineColorAlpha(kMagenta,0.35);
  h_BBgr_clus_size_pscut_notrackmatch ->SetFillColorAlpha(kMagenta,0.35);
  stack2->Add(h_BBgr_clus_size);
  stack2->Add(h_BBgr_clus_size_pscut);
  stack2->Add(h_BBgr_clus_size_pscut_notrackmatch);
  stack2 ->Draw("nostack hist");
    
  //Add a legend
  auto legend2 = new TLegend();
  legend2->AddEntry( h_BBgr_clus_size,"gr clus size","l");
  legend2->AddEntry( h_BBgr_clus_size_pscut,"ps.e > 0.2 && gr track matched","f");
  legend2->AddEntry( h_BBgr_clus_size_pscut_notrackmatch,"ps.e > 0.2 && not gr track matched","f");
  legend2->Draw();

  c2->cd(2);

  THStack *stack3 = new THStack("stack3", "GRINCH cluster size: ps.e <=0.2");
  h_BBgr_clus_size_psanticut ->SetLineColorAlpha(kGreen,0.35);
  h_BBgr_clus_size_psanticut ->SetFillColorAlpha(kGreen,0.35);
  h_BBgr_clus_size_psanticut_notrackmatch ->SetLineColorAlpha(kViolet,0.35);
  h_BBgr_clus_size_psanticut_notrackmatch ->SetFillColorAlpha(kViolet,0.35);
  stack3->Add(h_BBgr_clus_size);
  stack3->Add(h_BBgr_clus_size_psanticut);
  stack3->Add(h_BBgr_clus_size_psanticut_notrackmatch);
  stack3->Draw("nostack hist");

  //Add a legend
  auto legend3 = new TLegend();
  legend3->AddEntry( h_BBgr_clus_size,"gr clus size","l");
  legend3->AddEntry( h_BBgr_clus_size_psanticut,"ps.e <= 0.2 && gr track matched","f");
  legend3->AddEntry( h_BBgr_clus_size_psanticut_notrackmatch,"ps.e <= 0.2 && not gr track matched","f");
  legend3->Draw();


  TCanvas *c3 = new TCanvas("c3","c3",1000,500);
  c3->Divide(2,1);
  c3->cd(1);
   
  THStack *stack4 = new THStack("stack4", "GRINCH cluster size: track matched");
  h_BBgr_clus_size_pscut ->SetLineColorAlpha(kGreen+1,0.35);
  h_BBgr_clus_size_pscut ->SetFillColorAlpha(kGreen+1,0.35);
  h_BBgr_clus_size_psanticut ->SetLineColorAlpha(kBlue,0.35);
  h_BBgr_clus_size_psanticut ->SetFillColorAlpha(kBlue,0.35);
  stack4->Add(h_BBgr_clus_size);
  stack4->Add(h_BBgr_clus_size_pscut);
  stack4->Add(h_BBgr_clus_size_psanticut);
  stack4->Draw("nostack hist");
  //Add a legend
  auto legend4 = new TLegend();
  legend4->AddEntry( h_BBgr_clus_size,"gr clus size","l");
  legend4->AddEntry( h_BBgr_clus_size_pscut,"ps.e > 0.2 && gr track matched","f");
  legend4->AddEntry( h_BBgr_clus_size_psanticut,"ps.e <= 0.2 && gr track matched","f");
  legend4->Draw();

  c3->cd(2);
  THStack *stack5 = new THStack("stack5", "GRINCH cluster size: Not track matched");
  h_BBgr_clus_size_psanticut_notrackmatch ->SetLineColorAlpha(kRed,0.35);
  h_BBgr_clus_size_psanticut_notrackmatch ->SetFillColorAlpha(kRed,0.35);
  h_BBgr_clus_size_pscut_notrackmatch ->SetLineColorAlpha(kMagenta,0.35);
  h_BBgr_clus_size_pscut_notrackmatch ->SetFillColorAlpha(kMagenta,0.35);
  stack5 ->Add(h_BBgr_clus_size);
  stack5->Add(h_BBgr_clus_size_pscut_notrackmatch);
  stack5->Add(h_BBgr_clus_size_psanticut_notrackmatch);
  stack5->Draw("nostack hist");
    
  auto legend5 = new TLegend();
  legend5->AddEntry( h_BBgr_clus_size,"gr clus size","l");
  legend5->AddEntry( h_BBgr_clus_size_pscut_notrackmatch,"ps.e > 0.2 && Not gr track matched","f");
  legend5->AddEntry( h_BBgr_clus_size_psanticut_notrackmatch,"ps.e <= 0.2 && Not gr track matched","f");
  legend5->Draw();

  TCanvas *c4 = new TCanvas("c4","c4",1000,500);

  THStack *stack6 = new THStack("stack6", "preshower energy");
  h_BBps_nogrtrack ->SetLineColorAlpha(kRed,0.35);
  h_BBps_nogrtrack ->SetFillColorAlpha(kRed,0.35);
  h_BBps_granticut_totsum ->SetLineColorAlpha(kBlue,0.35);
  h_BBps_granticut_totsum ->SetFillColorAlpha(kBlue,0.35);
  h_BBps_grcut_totsum ->SetLineColorAlpha(kGreen+2,0.35);
  h_BBps_grcut_totsum ->SetFillColorAlpha(kGreen+2,0.35);
  h_BBps_e ->SetLineColorAlpha(kBlack,1);
  stack6->Add(h_BBps_nogrtrack);
  stack6 ->Add(h_BBps_e);
  stack6->Add(h_BBps_grcut_totsum);
  stack6->Add(h_BBps_granticut_totsum);
  stack6->Draw("nostack hist");

  //Add a legend
  auto legend6 = new TLegend();
  legend6->AddEntry(h_BBps_e,"bb.ps.e","l");
  legend6->AddEntry(h_BBps_grcut_totsum,Form("gr clus tot >= %.0f && track-matched",GR_cut),"f");
  legend6->AddEntry(h_BBps_granticut_totsum, Form("gr clus tot < %.0f && track-matched",GR_cut),"f");
  legend6->AddEntry(h_BBps_nogrtrack,"no track match","f");
  legend6->Draw();


  cout<< "Looped over " << max<< " entries." <<endl;
  fout -> Write();
}
//end main
