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
#include <stack>
#include <TBox.h>

/////////////////////////////////////////////////////
//// cluster_finding_hit.C
//// Maria Satik
//// January 2022
//// msatnik@email.wm.edu
////
//// hoping to make a simpler cluster finding for the GRINCH with this program
////
//global variables
stack <Int_t> short_row_on_left_search;
stack <Int_t> short_row_on_right_search;
stack <Int_t> long_row_on_left_search;
stack <Int_t> long_row_on_right_search;
stack <Int_t> middle_search;
stack <Int_t> top_row_middle_search;
stack <Int_t> top_row_left_search;
stack <Int_t> top_row_right_search;
stack <Int_t> bottom_row_middle_search;
stack <Int_t> bottom_row_left_search;
stack <Int_t> bottom_row_right_search;
stack <Int_t> clusterstack;

stack <Int_t> neighbors_stack;
stack <Int_t> cluster_stack;

const Int_t N_ROW=60;

Int_t row_array[510];
Int_t col_array[510];
Int_t tdc_to_adc_map_array[510]; // Needs to be updated when ADCs are moved to different Nino cards. 

Double_t tr_x[10];
Double_t tr_y[10];
Double_t tr_n;
Double_t tr_vz[100];
Double_t tr_tg_th[100];
Double_t tr_tg_ph[100];
Double_t gem_track_nhits;
Double_t tr_p[100];
Double_t hcal_e;

Double_t ps_e;  
Double_t sh_e;
Double_t pmom;
Double_t xptar; 
Double_t yptar;   
Double_t ytar; 
Int_t GrinchNum;
Double_t tdcGID[1000]; 
Double_t tdcGLe[1000]; 
Double_t tdcGTot[1000]; 
Double_t tdcGTe[1000]; 
Double_t tdcGMult[1000]; 

Double_t hodo_tmean;

Int_t GrinchNumHit;
Double_t tdcGHitLe[1000];
Double_t tdcGHitID[1000];
Double_t tdcGHitTot[1000];
Double_t tdcGHitTe[1000];
//is there one for multiplicity? or maybe you just use the tdcGMult...

Int_t GrinchADCNum;
Double_t adcGID[1000]; // ADC Channel Numbern
Double_t adcGAtime[1000]; // The time the signal arrives in the ADC window
Double_t adcGAmp[1000]; // Amplitude of the ADC signal
Double_t grinch_adc[1000];// Integral of the ADC signal
Double_t adcGMult[1000]; // ADC multiplity: how many times the signal went above threshold in the time window
Double_t gTrigBits; // Should be "1" if it is a bbcal trigger, 16 if it is a GRINCH LED event (I think)

 
Int_t hit_flag_array[510]; // array to help keep track of which PMTs were hit during an event
Int_t mult_array[510]; // array for the multiplicities of each PMT during an event
Int_t sum_array[510]; // array to help with cluster finding: records how many adjacent PMTs have a hit for each PMT
Int_t root_index_array[510]; // An array to help us convert betweeen the PMT number and the index that root uses for an event
Int_t adc_root_index_array[510]; // same as above but for ADC signals 

//histograms
TH1F* h_ps_e = new TH1F("h_ps_e"," ; Total PreShower E ",300,0.0,2.0);
TH1F* h_TBB_e = new TH1F("h_TBB_e"," ; Total PRe+Shower E ",300,0.0,5.0);
TH1F* h_ratio_e = new TH1F("h_ratio_e"," ; Total PRe+Shower E/ track norm ",300,0.0,2.5);
TH1F* h_gTrigBits = new TH1F("h_gTrigBits", "; gTrigBits;",35,0,35);

TH1F *h_hodo_tmean = new TH1F("h_hodo_tmean",";bb.hodotdc.clus.tmean;",40,-20,20); 
TH1F* h_tr_x = new TH1F("h_tr_x", "bb.tr.x[0]", 200,-1,1);
TH1F* h_tr_y = new TH1F("h_tr_y", "bb.tr.y[0]", 200, -1,1);
TH1F* h_tr_p = new TH1F("h_tr_p","bb.tr.p[0]",200,0,5);
TH1F* h_tr_vz = new TH1F("h_tr_vz","bb.tr.vz[0])",200,-0.2,0.2);
TH1F* h_hcal_e = new TH1F("h_hcal_e","sbs.hcal.e" ,200,0,1);
TH1F* h_gem_track_nhits = new TH1F("h_gem_track_nhits","bb.gem.track.nhits",10,0,10);
TH1F* h_tr_tg_ph = new TH1F("h_tr_tg_ph", "bb.tr.tg_ph[0]",200,-0.1,0.1);
TH1F* h_tr_tg_th = new TH1F("h_tr_tg_th", "bb.tr.tg_th[0]",200,-0.3,0.3);
TH1F* h_tr_n = new TH1F("h_tr_n","bb.tr.n",10,0,10);




TH2F* h_grinch_cluster_center_vert_tr = new TH2F("h_grinch_cluster_center_vert_tr",";GRINCH cluster center vertical ; track x;", 120,0,60,250,-1,1);
TH2F* h_grinch_cluster_center_horiz_tr = new TH2F("h_grinch_cluster_center_horiz_tr",";GRINCH cluster center horizontal; track y;", 18,0,9,75,-0.3,0.3);
TH2F* h_grinch_cluster_center_horiz_tr_display = new TH2F("h_grinch_cluster_center_horiz_tr",";cluster center horizontal; bb.tr.y[0];", 18,-9,0,75,-0.3,0.3);

TH2F* h_grinch_CUT_cluster_center_vert_tr = new TH2F("h_grinch_CUT_cluster_center_vert_tr",";cluster center vertical ; bb.tr.x[0];", 120,0,60,250,-1,1);
TH2F* h_grinch_CUT_cluster_center_horiz_tr = new TH2F("h_grinch_CUT_cluster_center_horiz_tr",";cluster center horizontal; bb.tr.y[0];", 18,0,9,75,-0.3,0.3);

TH2F* h_grinch_cluster_track_xy = new TH2F("h_grinch_cluster_track_xy",";bb.tr.y[0];bb.tr.x[0];", 60,-0.3,0.3,200,-1,1);

TH2F* h_good_event_track_xy = new TH2F("h_good_event_track_xy",";bb.tr.y[0];bb.tr.x[0];",  60,-0.3,0.3,200,-1,1);

TH2F* h_grinch_eff_xy = new TH2F("h_grinch_eff_xy",";bb.tr.y[0];bb.tr.x[0];",  8,-0.2,0.2,28,-0.7,0.7);

TH2F* h_grinch_cluster_track_subtraction = new TH2F("h_grinch_cluster_track_subtraction",";clusterx - tr_y[0]; clustery -tr_x[0];", 500,-1, 9,500,-1,60);
TH1F* h_grinch_cluster_center_vert_tr_subtraction =  new TH1F("h_grinch_cluster_center_vert_tr_subtraction", "clustery - tr_x[0]", 500,-1,60);
TH1F* h_grinch_cluster_center_horiz_tr_subtraction = new TH1F("h_grinch_cluster_center_horiz_tr_subtraction ", "clusterx - tr_y[0]", 500, -1,9);


TH2F* h_grinch_CUT_cluster_center_display = new TH2F("h_grinch_CUT_cluster_center_display", "; horizontal ; vertical ;",16,0,8,60,-60,0);


TH1F* h_grinch_cluster_size = new TH1F("h_grinch_cluster_size" ,"; Cluster Size ",12,3,15);
TH1F* h_grinch_cluster_size_good = new TH1F("h_grinch_cluster_size_good" ,"; Cluster Size (shower cut) ",12,3,15);
TH1F* h_grinch_cluster_size_bad = new TH1F("h_grinch_cluster_size_bad" ,"; Cluster Size (outside shower cut)",12,3,15);
TH2F* h_grinch_cluster_center = new TH2F("h_grinch_cluster_center", "; horizontal ; vertical ;",16,0,8,60,0,60);
TH2F* h_grinch_cluster_center_display = new TH2F("h_grinch_cluster_center_display", "; horizontal ; vertical ;",16,0,8,60,-60,0);



TH2F* h_grinch_cluster_center_spreadcut = new TH2F("h_grinch_cluster_center_spreadcut",";horizontal; vertical;",16,0,8,60,-60,0);
TH2F* test_coord_histo = new TH2F(" test_coord_histo", "; horizontal ; vertical ;",18,0,9,60,0,60);
TH1F* h_grinch_pmt_good_hit = new TH1F("h_grinch_pmt_good_hit",";PMT ;",510,0,510);
TH2F* h_grinch_cluster_spread = new TH2F("h_grinch_cluster_spread",";cluster width;  cluster height;",8,1,9,59,1,60);
TH1F* h_grinch_cluster_size_cuts = new TH1F("h_grinch_cluster_size_cuts" ,"; Cluster Size with cuts",12,3,15);
TH2F* h_grinch_cluster_center_bad_sh_e  = new TH2F("h_grinch_cluster_center_bad_sh_e ", "; horizontal ; vertical ;",16,0,8,60,-60,0);
TH2F* h_grinch_cluster_center_good_sh_e  = new TH2F("h_grinch_cluster_center_good_sh_e ", "; horizontal ; vertical ;",16,0,8,60,-60,0);

TH2F* h_grinch_cluster_center_3 = new TH2F("h_grinch_cluster_center_3", "; horizontal ; vertical ;",16,0,8,60,-60,0);
TH2F* h_grinch_cluster_center_4 = new TH2F("h_grinch_cluster_center_4", "; horizontal ; vertical ;",16,0,8,60,-60,0);
TH2F* h_grinch_cluster_center_5 = new TH2F("h_grinch_cluster_center_5", "; horizontal ; vertical ;",16,0,8,60,-60,0);
TH2F* h_grinch_cluster_center_6 = new TH2F("h_grinch_cluster_center_6", "; horizontal ; vertical ;",16,0,8,60,-60,0);
TH2F* h_grinch_cluster_center_7 = new TH2F("h_grinch_cluster_center_7", "; horizontal ; vertical ;",16,0,8,60,-60,0);



TH1F* h_grinch_cluster_elem = new TH1F("h_grinch_cluster_elem",";GRINCH TDC elemID;", 510,0,510);
TH2F* h_grinch_cluster_le_elem = new TH2F("h_grinch_cluster_le_elem","; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,200,800,1000);

TH2F* h_grinch_cluster_le_elem_hodo = new TH2F("h_grinch_cluster_le_elem_hodo","; GRINCH TDC elemID ; GRINCH TDC LE (ns) - hodo tmean ",510,0,510,200,800,1000);

TH2F* h_grinch_cluster_te_elem = new TH2F("h_grinch_cluster_te_elem"," ; GRINCH TDC elemID ; GRINCH TDC TE (ns)",510,0,510,200,800,1000);
TH2F* h_grinch_cluster_tot_elem = new TH2F("h_grinch_cluster_tot_elem"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) ",510,0,510,110,-10,100);
TH2F* h_grinch_cluster_mult_elem = new TH2F("h_grinch_cluster_mult_elem"," ; GRINCH TDC elemID ; GRINCH TDC Mult ",510,0,510,10,0,10);
TH2F* h_grinch_cluster_adcMult_elem = new TH2F("h_grinch_cluster_adcMult_elem"," ; GRINCH ADC elem id ; GRINCH ADC Mult. ",63,0,63,10,0,10);
TH2F* h_grinch_cluster_atime_elem = new TH2F("h_grinch_cluster_atime_elem"," ; GRINCH ADC elemID ; GRINCH ADC time (ns) ",63,0,63,500,0,2000);
TH2F* h_grinch_cluster_amp_elem = new TH2F("h_grinch_cluster_amp_elem"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) ",63,0,63,1000,-100,900);    
TH2F* h_grinch_cluster_int_elem = new TH2F("h_grinch_cluster_int_elem","; GRINCH ADC elemID ; GRINCH ADC channel ",64,0,63,400,0,199);//maria
TH2F* h_grinch_cluster_adc_tdc =new TH2F("h_grinch_cluster_adc_tdc","; GRINCH TDC ToT ; GRINCH ADC AMP ",60,0,60,200,0,200);
TH2F* h_grinch_cluster_adcInt_tdc =new TH2F("h_grinch_cluster_adcInt_tdc","; GRINCH TDC ToT ; GRINCH ADC INT ",60,0,60,150,0,150);
TH2F* h_grinch_cluster_int_amp =new TH2F("h_grinch_cluster_int_amp","; GRINCH ADC AMP ; GRINCH ADC INT ",200,0,200,150,0,150);
TH2F* h_grinch_cluster_adc_tdc_test =new TH2F("h_grinch_cluster_adc_tdc_test","; GRINCH TDC ToT ; GRINCH ADC AMP ",60,0,60,200,0,200);

TH2F* h_grinch_cluster_bad_le_elem = new TH2F("h_grinch_cluster_bad_le_elem","; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,600,600,1200);
TH2F* h_grinch_cluster_bad_te_elem = new TH2F("h_grinch_cluster_bad_te_elem"," ; GRINCH TDC elemID ; GRINCH TDC TE (ns)",510,0,510,600,600,1200);
TH2F* h_grinch_cluster_bad_tot_elem = new TH2F("h_grinch_cluster_bad_tot_elem"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) ",510,0,510,110,-10,100);
TH2F* h_grinch_cluster_bad_mult_elem = new TH2F("h_grinch_cluster_bad_mult_elem"," ; GRINCH TDC elemID ; GRINCH TDC Mult ",510,0,510,10,0,10);
TH2F* h_grinch_cluster_bad_adcMult_elem = new TH2F("h_grinch_cluster_bad_adcMult_elem"," ; GRINCH ADC elem id ; GRINCH ADC Mult. ",63,0,63,10,0,10);
TH2F* h_grinch_cluster_bad_atime_elem = new TH2F("h_grinch_cluster_bad_atime_elem"," ; GRINCH ADC elemID ; GRINCH ADC time (ns) ",63,0,63,500,0,2000);
TH2F* h_grinch_cluster_bad_amp_elem = new TH2F("h_grinch_cluster_bad_amp_elem"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) ",63,0,63,500,0,500);    
TH2F* h_grinch_cluster_bad_int_elem = new TH2F("h_grinch_cluster_bad_int_elem","; GRINCH ADC elemID ; GRINCH ADC channel ",64,0,63,510,-10,500);//maria
TH2F* h_grinch_cluster_bad_adc_tdc =new TH2F("h_grinch_cluster_bad_adc_tdc","; GRINCH TDC ToT ; GRINCH ADC AMP ",60,0,60,200,0,200);
TH2F* h_grinch_cluster_bad_adcInt_tdc =new TH2F("h_grinch_cluster_bad_adcInt_tdc","; GRINCH TDC ToT ; GRINCH ADC INT ",60,0,60,150,0,150);
TH2F* h_grinch_cluster_bad_int_amp =new TH2F("h_grinch_cluster_bad_int_amp","; GRINCH ADC AMP ; GRINCH ADC INT ",200,0,200,150,0,150);
TH2F* h_grinch_cluster_bad_adc_tdc_test =new TH2F("h_grinch_cluster_bad_adc_tdc_test","; GRINCH TDC ToT ; GRINCH ADC AMP ",60,0,60,200,0,200);

TGraph* background_graph = new TGraph();
TGraph* rate_graph = new TGraph();
TGraphErrors* offset_graph = new TGraphErrors();
TGraphErrors* mean_graph = new TGraphErrors();
TGraphErrors* sigma_graph = new TGraphErrors();
TGraphErrors* sigma_squared_graph = new TGraphErrors();


TH2F* h_grinch_le_elem = new TH2F("h_grinch_le_elem"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,2600,0,2600);
TH2F* h_grinch_tot_elem = new TH2F("h_grinch_tot_elem"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) ",510,0,510,101,0,100);

TH2F* h_grinch_le_elem_hodo = new TH2F("h_grinch_le_elem_hodo"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) - hodo tmean ",510,0,510,2600,0,2600);


TH2F* h_grinch_mult_elem = new TH2F("h_grinch_mult_elem"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,10,0,10);
TH1F* h_grinch_elem = new TH1F("h_grinch_elem",";GRINCH elem ID;",510,0,510);
TH2F* h_grinch_amp_elem = new TH2F("h_grinch_amp_elem"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) ",63,0,63,1000,-100,1900);    

TH2F* h_grinch_background_elem = new TH2F("h_grinch_background_elem","; GRINCH TDC elemID; Background",510,0,510,600,0,600);

TH2F* h_grinch_rate_elem = new TH2F("h_grinch_rate_elem"," ; GRINCH TDC elemID; calc. bkgnd rate",510,0,510,10000,0,1000000);

TH1F* h_grinch_hitslistlength = new TH1F("h_grinch_hitslistlength",";total fires in the multi-hit tdc;", 500,0,500);

//TH2F* h_grinch_hits_elem = new TH2F("h_grinch_hits_elem", " ; ; ")

TH1F* h_grinch_le_all = new TH1F("h_grinch_le_all","; GRINCH LE ALL ;", 2500,0,2500);

TH1F* h_grinch_tot_all =new TH1F("h_grinch_tot_all","; GRINCH TOT ALL ;", 101,0,100);

TH2F* h_grinch_hit_le_elem = new TH2F("h_grinch_hit_le_elem"," ; GRINCH TDC elemID ; GRINCH TDC HIT LE (ns) ",510,0,510,2700,-100,2600);
TH1F* h_grinch_hit_le_all = new TH1F("h_grinch_hit_le_all","; GRINCH LE HIT ALL ;", 2500,0,2500);

TH1F* h_grinch_cluster_le_all = new TH1F("h_grinch_cluster_le_all","; GRINCH CLUSTER LE ALL ;", 2500,0,2500);
TH1F* h_grinch_cluster_tot_all = new TH1F("h_grinch_cluster_tot_all","; GRINCH CLUSTER TOT ALL ;", 101,0,100);

TH1F* h_grinch_cluster_tot[511];
TH1F* h_grinch_cluster_le[511];
TH1F* h_grinch_cluster_amp[64];//maria 
TH1F* h_grinch_amp[64];

Int_t sh_ps_hit_cnt = 0; //counter for the number of events that pass the shower cut
Int_t sh_hit_no_grinch_hit_cnt = 0; //counter for the number of events that pass the shower cut but do not have a grinch cluster
Int_t grinch_hit_sh_hit_cnt = 0; // counter for the number of event that pass the shower cut and have a grinch cluster
Int_t grinch_hit_no_sh_hit_cnt = 0;// counter for the number of events that have a grinch cluster but do not pass the shower cut
Int_t no_sh_no_grinch_cnt = 0;// counter for the number of events that do not pass the shower cut and do not have a grinch cluster
Int_t cluster_in_tr_cut_cnt = 0;
Int_t tr_cut_cnt = 0; 

Int_t cluster_in_tr_cut_cnt_2D[100][50] = {0};
Int_t tr_cut_cnt_2D[100][50] = {0};
Double_t nsteps_trx = 28;
Double_t nsteps_try = 8;
Double_t stepsize = 0.05;
Double_t min_trx = -0.7;
Double_t min_try = -0.2;
Double_t trx_steps_array[100] = {0};
Double_t try_steps_array[100] = {0};



TH1F* h_grinch_sh_ps_hit_cnt;
TH1F* h_grinch_sh_hit_no_grinch_hit_cnt;
TH1F* h_grinch_grinch_hit_sh_hit_cnt;
TH1F* h_grinch_grinch_hit_no_sh_hit_cnt;
TH1F* h_grinch_no_sh_no_grinch_cnt; 
TH1F* h_grinch_efficiency_percent = new TH1F("h_grinch_efficiency_percent" ,"; percent efficiency",100,0,100);

//function declarations
void make_row_col_array(); // this is probably already in the database and isn't needed
void Fill_Search_Stacks();
Bool_t Is_NOT_In_Stack(stack <Int_t> checkstack, Int_t PMT);
void Print_Stack(stack <Int_t> inputstack);
stack <Int_t> Clear_Stack();
Int_t Sum_Adjacent_Hits(Int_t PMT, stack <Int_t>& neighbors_stack);
stack <Int_t> Add_to_Stack(stack <Int_t> inputstack, Int_t PMT);
void Calculate_PMT_Coord(double& horiz, double& vert, Int_t PMT); // May be unecessary: already in database?
void Calculate_Cluster_Center(stack <Int_t> inputstack, double& horiz, double& vert);
void Find_Cluster_Extrema(stack <Int_t> inputstack, Int_t& horizspread, Int_t& vertspread);
void Fill_Cluster_Histos(stack <Int_t> inputstack);
void Fill_Bad_Cluster_Histos(stack <Int_t> inputstack);
void make_tdc_to_adc_map_array(); // Probably uncessary: should add to database instead
Int_t Stack_Size(stack <Int_t> inputstack);
Double_t offset_gaus(Double_t *x, Double_t *par);
void Palette1();



void cluster_finding_hit(TString basename="",Int_t nrun=2043,TString configfilename="run_list_grinch.txt", Bool_t show_track = kFALSE){ //MAIN
  if (basename=="") {
    cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
    cin >> basename;
  }
  TString fullname = Form("%s_%d",basename.Data(),nrun);
  // gStyle->SetPalette(kBlueRedYellow);
  Palette1();
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
  TString inputroot;
  inputroot="rootfiles/"+fullname+".root";
  TString outputhist;
  outputhist= "hist/"+fullname+"_grinch_hist.root";
  cout<<" writing to file "<<outputhist<<endl;
  TObjArray HList(0);
  TString Chainroot;
  configfilename="RunLists/"+configfilename;
  ifstream configfile(configfilename.Data());
  TString currentline;
  TChain *fchain = new TChain("T");
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ) {
    if( !currentline.BeginsWith("#") ){
      cout << " add file : " << currentline << endl;
      fchain->Add(currentline);
    }   
  } 
  //  Int_t hcalNum;
  //fchain->SetBranchAddress("Ndata.sbs.hcal.tdcelemID",&hcalNum) ;
  
 
  fchain->SetBranchAddress("bb.ps.e",&ps_e) ; // preshower energy
  fchain->SetBranchAddress("bb.sh.e",&sh_e) ;  // shower energy
  fchain->SetBranchAddress("BB.gold.p",&pmom) ; //? 
  fchain->SetBranchAddress("BB.gold.th",&xptar) ; //?
  fchain->SetBranchAddress("BB.gold.ph",&yptar) ; //?
  fchain->SetBranchAddress("BB.gold.y",&ytar) ; //?

  fchain ->SetBranchAddress("bb.tr.x",&tr_x);
  fchain ->SetBranchAddress("bb.tr.y",&tr_y);
  fchain ->SetBranchAddress("bb.tr.n", &tr_n);
  fchain ->SetBranchAddress("bb.tr.vz", &tr_vz);
  fchain ->SetBranchAddress("bb.tr.tg_th",&tr_tg_th);
  fchain ->SetBranchAddress("bb.tr.tg_ph",&tr_tg_ph);
  fchain ->SetBranchAddress("bb.gem.track.nhits",&gem_track_nhits);
  fchain ->SetBranchAddress("bb.tr.p",&tr_p);
  fchain ->SetBranchAddress("sbs.hcal.e", &hcal_e);

 

  fchain ->SetBranchAddress("bb.hodotdc.clus.tmean",&hodo_tmean); 

  
  // fchain->SetBranchAddress("Ndata.bb.grinch_tdc.tdcelemID",&GrinchNum) ;//number of PMTs with a signal in an event (I think)  
  // fchain->SetBranchAddress("bb.grinch_tdc.tdcelemID",&tdcGID) ; // The PMT number 
  // fchain->SetBranchAddress("bb.grinch_tdc.tdc",&tdcGLe) ; // leading edge
  // fchain -> SetBranchAddress("bb.grinch_tdc.tdc_te",&tdcGTe) ; // trailing edge
  // fchain->SetBranchAddress("bb.grinch_tdc.tdc_tot",&tdcGTot) ; // tdc time over threshold (te-le)
  
  fchain->SetBranchAddress("bb.grinch_tdc.tdc_mult",&tdcGMult) ; //tdc multiplicity
  
  ////I need to go back to a-onl and see how I was looking at all the hits and get that set up in my ifarm   
  ////for when we are looking at all the TDC hits in the window.  grinch_tdc->SetStoreRawHits(kTRUE) in replay_GRINCH
  fchain->SetBranchAddress("bb.grinch_tdc.hits.TDCelemID", &tdcGID); //pmt number
  fchain->SetBranchAddress("Ndata.bb.grinch_tdc.hits.TDCelemID", &GrinchNum); 
  fchain->SetBranchAddress("bb.grinch_tdc.hits.t",&tdcGLe); //leading edge for all hits in the window
  fchain->SetBranchAddress("bb.grinch_tdc.hits.t_te",&tdcGTe);
  fchain->SetBranchAddress("bb.grinch_tdc.hits.t_tot",&tdcGTot);


  fchain->SetBranchAddress("Ndata.bb.grinch_adc.adcelemID",&GrinchADCNum) ; //number of PMTs with ADC channels with a signal in an event. 
  fchain->SetBranchAddress("bb.grinch_adc.adcelemID",&adcGID) ; // The ADC channel 
  fchain->SetBranchAddress("bb.grinch_adc.a_time",&adcGAtime) ; // the time the first(?)ADC signal went over threshold in the window  
  fchain->SetBranchAddress("bb.grinch_adc.a_amp_p",&adcGAmp) ; // pedestal-subtracted amplitude
  fchain->SetBranchAddress("bb.grinch_adc.a_mult",&adcGMult) ;// ADC multiplicity 
  fchain->SetBranchAddress("bb.grinch_adc.a_p",&grinch_adc) ; // ADC Integral (I think)

  fchain->SetBranchAddress("g.trigbits",&gTrigBits); // The type of trigger the event was. (1 is the bbcal trig, 16 is the grinch LED) (May not be working)
  
  
 
  //Making histos for the individual channels 
  for (Int_t ig=0;ig<511;ig++) {    
    h_grinch_cluster_tot[ig] = new TH1F(Form("h_grinch_cluster_tot_%d",ig),Form(" ; GRINCH TDC ToT CLUSTER PMT %d ; ",ig),50,0,50);
    HList.Add(h_grinch_cluster_tot[ig]);
      h_grinch_cluster_le[ig] = new TH1F(Form("h_grinch_cluster_le_%d",ig),Form(" ; GRINCH TDC LE CLUSTER PMT %d ; ",ig),50,875,925);
    HList.Add(h_grinch_cluster_le[ig]);
  }
 
  for (Int_t i=0;i<64;i++){   
    h_grinch_cluster_amp[i] = new TH1F(Form("h_grinch_cluster_amp_%d",i),Form(" ; ADC Amp CLUSTER %d ; ",i), 501, -0.5, 500.5); 
    h_grinch_amp[i] = new TH1F(Form("h_grinch_amp_%d",i),Form(" ; ADC Amp %d ; ",i), 501, -0.5, 500.5);
    HList.Add(h_grinch_cluster_amp[i]);
  }
    
 
  //// functions that initialize a few things 
  Fill_Search_Stacks(); //sets up cluster finding stacks
  make_row_col_array(); //sets up an array to quickly get the row/col from PMT number
  make_tdc_to_adc_map_array();//sets up arrays to quickly convert between TDC number and ADC number 


  
  Double_t value_x = 0;
  Double_t value_y = 0;

  for (Int_t stepcnt_x = 0; stepcnt_x < nsteps_trx+1; stepcnt_x ++)
    {
      value_x = min_trx + stepsize*stepcnt_x;
      if (value_x >-0.0001 && value_x <0.0001){value_x = 0;}//// it was saving zero as 1.22 E-16, so I hardcoded it to zero. 
      trx_steps_array[stepcnt_x] =value_x;
      cout<<"trx_steps_array["<<stepcnt_x<<"] = "<<value_x<<endl;
    }
  
  for (Int_t stepcnt_y = 0; stepcnt_y <nsteps_try+1; stepcnt_y ++)
    {
      value_y = min_try + stepsize*stepcnt_y;
      try_steps_array[stepcnt_y] = value_y;
      cout<< "try_steps_array["<<stepcnt_y <<"] = "<<value_y<<endl;
    }
  
  
  
  /*
  for (Int_t pmt = 0; pmt<510;pmt++)
    {
      Int_t adc_chan = tdc_to_adc_map_array[pmt];
      if(adc_chan != -1)
	{
	  cout << "ADC Chan "<<adc_chan<< " is on PMT "<<pmt<<endl;
	}
    }
  */

  //// adding histograms to the "HList" so that we can look at them later using macros
  HList.Add(h_ps_e);  
  HList.Add(h_TBB_e); 
  HList.Add(h_ratio_e); 
 
  HList.Add(h_gTrigBits);

  HList.Add(h_grinch_cluster_size);
  HList.Add(h_grinch_cluster_size_good);
  HList.Add(h_grinch_cluster_size_bad);
  HList.Add(h_grinch_cluster_center);
  HList.Add(h_grinch_cluster_center_display);
  HList.Add(test_coord_histo); 
  HList.Add(h_grinch_pmt_good_hit);
  HList.Add(h_grinch_cluster_spread);
  HList.Add(h_grinch_cluster_center_spreadcut);
  HList.Add(h_grinch_cluster_size_cuts);
  HList.Add(h_grinch_cluster_center_bad_sh_e);
  HList.Add(h_grinch_cluster_center_good_sh_e);
  HList.Add(h_grinch_cluster_adc_tdc_test);

  HList.Add(h_grinch_cluster_center_3);
  HList.Add(h_grinch_cluster_center_4);
  HList.Add(h_grinch_cluster_center_5);
  HList.Add(h_grinch_cluster_center_6);
  HList.Add(h_grinch_cluster_center_7);

  HList.Add(h_grinch_cluster_le_elem);
  HList.Add(h_grinch_cluster_te_elem);
  HList.Add(h_grinch_cluster_tot_elem);
  HList.Add(h_grinch_cluster_mult_elem);
  HList.Add(h_grinch_cluster_amp_elem);
  HList.Add(h_grinch_cluster_int_elem);
  HList.Add(h_grinch_cluster_atime_elem);
  HList.Add(h_grinch_cluster_adcMult_elem);
  HList.Add(h_grinch_cluster_adc_tdc);
  HList.Add(h_grinch_cluster_adcInt_tdc);
  HList.Add(h_grinch_cluster_int_amp);

  HList.Add(h_grinch_cluster_bad_le_elem);
  HList.Add(h_grinch_cluster_bad_te_elem);
  HList.Add(h_grinch_cluster_bad_tot_elem);
  HList.Add(h_grinch_cluster_bad_mult_elem);
  HList.Add(h_grinch_cluster_bad_amp_elem);
  HList.Add(h_grinch_cluster_bad_int_elem);
  HList.Add(h_grinch_cluster_bad_atime_elem);
  HList.Add(h_grinch_cluster_bad_adcMult_elem);
  HList.Add(h_grinch_cluster_bad_adc_tdc);
  HList.Add(h_grinch_cluster_bad_adcInt_tdc);
  HList.Add(h_grinch_cluster_bad_int_amp);
  HList.Add(h_grinch_cluster_bad_adc_tdc_test);

  HList.Add(h_grinch_sh_ps_hit_cnt);
  HList.Add(h_grinch_sh_hit_no_grinch_hit_cnt);
  HList.Add(h_grinch_grinch_hit_sh_hit_cnt);
  HList.Add(h_grinch_grinch_hit_no_sh_hit_cnt);
  HList.Add(h_grinch_no_sh_no_grinch_cnt);
  HList.Add(h_grinch_efficiency_percent);

  HList.Add(h_grinch_le_elem);
  HList.Add(h_grinch_le_all); 

  HList.Add(h_grinch_cluster_le_all); 

  //HList.Add(h_grinch_hit_le_elem);
  // HList.Add(h_grinch_hit_le_all);


  TH1F* h_grinch_hit_le[511];
  TH1F* h_grinch_hit_tot[511];
  for (Int_t ig=0;ig<511;ig++) {
    h_grinch_hit_le[ig] = new TH1F(Form("h_grinch_hit_le_%d",ig),Form(" ; GRINCH TDC HIT LE (ns) PMT %d  ; ",ig),2600,0,2600);
    h_grinch_hit_tot[ig] = new TH1F(Form("h_grinch_hit_tot_%d",ig),Form(" ; GRINCH TDC HIT TOT (ns) PMT %d  ; ",ig),50,0,50);
    HList.Add(h_grinch_hit_le[ig]);  
  }

  TH1F* h_grinch_le[511];
  TH1F* h_grinch_tot[511];
  TH1F* h_grinch_hodo_le[510];
   for (Int_t ig=0;ig<511;ig++) {
    h_grinch_le[ig] = new TH1F(Form("h_grinch_le_%d",ig),Form(" ; GRINCH TDC LE (ns) PMT %d  ; ",ig),2600,0,2600);
    h_grinch_tot[ig] = new TH1F(Form("h_grinch_tot_%d",ig),Form(" ; GRINCH TDC ToT (ns) PMT %d  ; ",ig),50,0,50);
    HList.Add(h_grinch_le[ig]);    

    h_grinch_hodo_le[ig] = new TH1F(Form("h_grinch_hodo_le_%d",ig),Form(" ; GRINCH TDC LE (ns) -hodo_tmean  PMT %d  ; ",ig),2600,0,2600);
  }



  Long64_t nentries = fchain->GetEntries();
  cout<<"nentries: "<<nentries<<endl;

  Int_t temp_cnt=0;

  Int_t zero_count = 0; 
  Int_t clustercnt = 0;
  Int_t max = nentries; // make this lower if you don't want to analyze all of the entries
  if (max > nentries){ max = nentries;}
  cout<<"max = "<<max<<endl;

  Int_t hit_sum = 0; 

  
  //// filling out a histogram to show the x-y coordinates of each PMT
  double testx=0;
  double testy=0;  
  for (Int_t pmtnum = 0; pmtnum<510;pmtnum++)
    {
     Calculate_PMT_Coord(testx, testy, pmtnum);
     test_coord_histo ->Fill(testx,testy);
    }
 
 

  //// Loop over the number of entries /////////
  for (int i = 0; i < max ; i++) {//nentries
    fchain->GetEntry(i);
    if (i%5000==0) cout << " Entry = " << i << endl;
    
    Double_t tot_e = ps_e+sh_e;
    Double_t rat = tot_e/pmom;
    h_ps_e->Fill(ps_e);
    h_TBB_e->Fill(tot_e);
    h_ratio_e->Fill(rat);

    h_hodo_tmean ->Fill(hodo_tmean);

    h_tr_x ->Fill(tr_x[0]);
    h_tr_y ->Fill(tr_y[0]);

    h_tr_n ->Fill(tr_n);
    h_tr_vz ->Fill(tr_vz[0]);
    h_tr_tg_th ->Fill(tr_tg_th[0]);
    h_tr_tg_ph ->Fill(tr_tg_ph[0]);
    h_tr_p ->Fill(tr_p[0]);
    h_hcal_e ->Fill(hcal_e);
    h_gem_track_nhits ->Fill(gem_track_nhits);
    
    
    h_gTrigBits ->Fill(gTrigBits);

    if(gTrigBits!=1)// gTrigBits == 1 is the bbcal trigger, == 16 is the grinch LED. 
      ///Only want to process bbcal trigger events.
      {
    	continue; //breaks out and goes to next entry 
      }

    Bool_t sh_flag = kFALSE; //flag for when an event passes the shower cut
    Bool_t g_cluster_flag =kFALSE;//flag for when an event has a grinch cluster

    //// initialize some arrays
    for(Int_t n = 0; n< 64; n++) //initilalize root index array
      {
	adc_root_index_array[n]=0;
      }

    
    //if (ps_e > 0.2 && abs(rat - 1.15) < 0.1 && abs(tot_e -3.15) < 0.5)// shower cut. abs(tot_e - 3.25)< 0.5 && ps_e > 0.2 && abs(rat - 1.2) < 0.2
    // if(tr_n ==1 && abs(tr_vz[0])<0.08 && abs(tr_tg_th[0])<0.15 && abs(tr_tg_ph[0])<0.3 && gem_track_nhits > 3 && tr_p[0] >3.0 && tr_p[0]<4 && hcal_e > 0.025 && ps_e >0.22) //SBS 8 cuts 
        if(tr_n ==1 && abs(tr_vz[0])<0.08 && gem_track_nhits > 3 && tr_p[0] >1.4 && tr_p[0]<2.0 && hcal_e > 0.025 && ps_e >0.22) //SBS 9 cuts && abs(tr_tg_th[0])<0.15 && abs(tr_tg_ph[0])<0.3 &&
       {
	////
	sh_flag = kTRUE; //throw up a flag that this is a good electron hit 
	sh_ps_hit_cnt++;
	////
	
	 h_good_event_track_xy->Fill(tr_y[0], tr_x[0]);  

	if(tr_x[0] <- 0.12 && tr_x[0] > -0.14 && tr_y[0] > -0.13 && tr_y[0] <-0.11){
	  tr_cut_cnt ++;
	}

	for (Int_t stepcnt_x = 0; stepcnt_x < nsteps_trx; stepcnt_x ++)
	  {
	    for (Int_t stepcnt_y = 0; stepcnt_y < nsteps_try+1; stepcnt_y ++)
	      {
		if( tr_x[0] > trx_steps_array[stepcnt_x] && tr_x[0] <  trx_steps_array[stepcnt_x+1] && tr_y[0] > try_steps_array[stepcnt_y] && tr_y[0] <try_steps_array[stepcnt_y+1]){
		  tr_cut_cnt_2D[stepcnt_x][stepcnt_y] ++;
		}		
	      }
	  }
	
       }

    //// loop over the PMTs with ADC's that had a signal for this event
    for (Int_t ig=0;ig<GrinchADCNum;ig++) { //since PMTs with no signal are not written to the root file, this "ig" is NOT the ADC channel
      Int_t gindex_adc=adcGID[ig]; // we need to specifically ask which ADC channel this signal is for         

      //should have if statment on atime or on a good tdc hit
      h_grinch_amp_elem ->Fill(gindex_adc,adcGAmp[ig]);
      h_grinch_amp[gindex_adc] ->Fill(adcGAmp[ig]);

      //// Map the "gindex_adc" (which is the ADC Chan number) to the "ig index" which is what the branch needs.
      adc_root_index_array[gindex_adc]=ig;
      //// need this in order to go back after finding clusters to look at ADC values and whatnot.          
    }// end loop over ADC 

    for (int n = 0; n<510;n++) //initialize some more arrays
      {
	hit_flag_array[n]=0;
	mult_array[n]=0;
	root_index_array[n]=0;
      }

    /*

    //// ----------------------------------------------------------------------
    ////Attempting to read out the multiple hits on the TDC.    
    for (Int_t igHit = 0; igHit < GrinchNumHit; igHit++)
      {	
	Int_t gindexHit = tdcGHitID[igHit];	

	if(gindexHit<511 && gindexHit>-1){
	    h_grinch_hit_le_all ->Fill(tdcGHitLe[igHit]);
	    h_grinch_hit_le_elem -> Fill(gindexHit, tdcGHitLe[igHit]);
	    h_grinch_hit_le[gindexHit] ->Fill(tdcGHitLe[igHit]);
	}
      }
    //// -------------------------------------------------------------------------
    */
    Double_t grinch_minus_hodo_time;
    Double_t goodhit=0;
    hit_sum = hit_sum + GrinchNum; //adding up all of the hits from the multi-hit TDC to normalize later
    h_grinch_hitslistlength ->Fill(GrinchNum);
    //// Loop over the PMTs that had a signal for this event 
    for (Int_t ig=0;ig<GrinchNum;ig++) { //since PMTs with no signal are not written to the root file, this "ig" is NOT the PMT number
      Int_t gindex=tdcGID[ig]; // we need to specifically ask which PMT this signal is for
      h_grinch_elem->Fill(gindex);
      if (gindex <511 && gindex >-1) {
	h_grinch_le_elem->Fill(tdcGID[ig],tdcGLe[ig] - hodo_tmean);
	h_grinch_le[gindex]->Fill(tdcGLe[ig] - hodo_tmean); // 
	h_grinch_hodo_le[gindex]->Fill(tdcGLe[ig] -hodo_tmean); //
	
	h_grinch_le_all ->Fill(tdcGLe[ig]-hodo_tmean);
	h_grinch_tot_elem ->Fill(tdcGID[ig],tdcGTot[ig]);
	h_grinch_tot_all ->Fill(tdcGTot[ig]);
	h_grinch_tot[gindex] ->Fill(tdcGTot[ig]);
	
	grinch_minus_hodo_time = tdcGLe[ig] - hodo_tmean;
	h_grinch_le_elem_hodo ->Fill(tdcGID[ig], grinch_minus_hodo_time);


	// note for eric: I belive this cut is similar to the one in SBSGRINCH. So it would be like we are starting here. 
	if (abs(tdcGLe[ig]-900)<25 && sh_flag){ // if the event is within the good timing cut for the grinch 
	  goodhit++;
	  ////////////////////////////
	  hit_flag_array[gindex] = 1; // mark that this PMT has a good hit
	  mult_array[gindex] = tdcGMult[ig]; // save it's multiplicity to the array.
	  h_grinch_mult_elem ->Fill(gindex, tdcGMult[ig]);
	  root_index_array[gindex]=ig; // Map the "gindex" (which is the PMT number) to the "ig index" which is what the branch needs.
	                               // I want to be able to go back after finding clusters to look at ADC and TDC and whatnot on those PMTs.   
	  h_grinch_pmt_good_hit ->Fill(gindex);
	  ////////////////////////////
	}
     
      }         
     
    }//end loop over TDC PMTs
  

    clusterstack = Clear_Stack(); //clearing the stack we are using for cluster finding for the new event
    
    
   
    if (temp_cnt <max){// temporary loop here to only do cluster finding on a certain number of events for debugging. 
      for(int m = 0; m<510;m++) //fill with 0 
	{
	  sum_array[m]=0;
	}

      //// Start of Cluster Finding
      
      cluster_stack = Clear_Stack();
      Int_t top;
      for (int pmtcntr = 9; pmtcntr <=500 ; pmtcntr ++) // algorithm doesn't need the top row and bottom row.
	{
	  neighbors_stack = Clear_Stack();
	  sum_array[pmtcntr] = Sum_Adjacent_Hits(pmtcntr,neighbors_stack); // retuns the number PMTs with hits adjacent to the PMT we are looking at	 
	
	  if (sum_array[pmtcntr] != 0)
	    {
	      //cout<<"sum_array [PMT] "<<pmtcntr<< " = "<<sum_array[pmtcntr]<<endl;
	    }
	  if (sum_array[pmtcntr] >= 2) // if the PMT has 2 or more neighboring PMTs with hits
	    {
	      //cout<< "Possible cluster around PMT "<<  pmtcntr << ": sum = " <<sum_array[pmtcntr]<<endl;	      
	      clusterstack.push(pmtcntr); // put that PMT in the stack to be analyzed to get cluster size
	      //Print_Stack(neighbors_stack);
	      while(!neighbors_stack.empty())
		{
		  top = neighbors_stack.top();
		  if(Is_NOT_In_Stack(cluster_stack, top ) )
		    {
		      cluster_stack.push(top);
		    }
		  neighbors_stack.pop();
		}
	    }
	}

     

      Int_t ClusterSize = Stack_Size(cluster_stack);
      
      temp_cnt ++; 

      

      g_cluster_flag =kFALSE;
      if(ClusterSize >= 3)
	{
	  //cout<<endl;
	  // cout<<"EVENT "<<i<<endl;
	  //Print_Stack(cluster_stack);
	  //cout<<"cluster_stack"<<endl;
	  //cout<< "cluster size = "<< ClusterSize<<endl;

	  clustercnt++;
	  g_cluster_flag = kTRUE;

	  h_grinch_cluster_size ->Fill(ClusterSize);
	
	
	  //cout<<"clustercnt = "<<clustercnt<<endl;
	  
      	 	  
	  /*
	    stack <Int_t> copystack = cluster_stack;
	    Int_t pmttest=0;
	    double x=0;
	    double y=0;
	    Int_t test = 0;
	    while(!copystack.empty())
	    {
	    //pmttest = testcnt;
	    //Calculate_PMT_Coord(x,y,PMT);
	    test = copystack.top();
	    copystack.pop();	   
	    Calculate_PMT_Coord(x, y, test);	  
	    cout<< "PMT "<<test << " is at "<< x << ", "<<y<<endl;
	    //cout<<"PMT "<<PMT<<" is at x = "<< x << " y = "<<y<<endl;
	    }
	  */
	  
	  
	  
	  double clusterx = 0;
	  double clustery = 0;  
	  
	  Calculate_Cluster_Center(cluster_stack,clusterx,clustery); // clusterx and clustery are reference parameters so the function updates them.
	  //// cluster_stack is a global variable that was filled in Calculate_Cluster_Size. I should prob change it to a referene param so it's not so confusing 
	  //cout<<"cluster center is at "<<clusterx<<", "<<clustery<<endl;
	  h_grinch_cluster_center ->Fill(clusterx,clustery);
	  h_grinch_cluster_center_display ->Fill(clusterx,-clustery);

	  h_grinch_cluster_center_horiz_tr->Fill(clusterx, tr_y[0]); //horizontal
	  h_grinch_cluster_center_horiz_tr_display->Fill(-clusterx, tr_y[0]); //horizontal with neg clusterx
	  h_grinch_cluster_center_vert_tr ->Fill(clustery, tr_x[0]); //vertical 

	  h_grinch_cluster_center_vert_tr_subtraction ->Fill( clustery - tr_x[0]);
	  h_grinch_cluster_center_horiz_tr_subtraction -> Fill(clusterx - tr_y[0]);
	  h_grinch_cluster_track_subtraction ->Fill (clusterx - tr_y[0] , clustery-tr_x[0]); 
	  
	
	  h_grinch_cluster_track_xy ->Fill(tr_y[0], tr_x[0]);  
	  
	  if(tr_x[0] <- 0.12 && tr_x[0] > -0.14 && tr_y[0] > -0.13 && tr_y[0] <-0.11)
	    {//looking at just a small region in the tr_x tr_y plane
	      cluster_in_tr_cut_cnt ++;
	      h_grinch_CUT_cluster_center_horiz_tr->Fill(clusterx, tr_y[0]); //horizontal
	      h_grinch_CUT_cluster_center_vert_tr ->Fill(clustery, tr_x[0]); //vertical 
	      h_grinch_CUT_cluster_center_display  ->Fill(clusterx,-clustery);
	    }

	  for (Int_t stepcnt_x = 0; stepcnt_x < nsteps_trx; stepcnt_x ++)
	    {
	      for (Int_t stepcnt_y = 0; stepcnt_y < nsteps_try+1; stepcnt_y ++)
		{
		  if( tr_x[0] >= trx_steps_array[stepcnt_x] && tr_x[0] <  trx_steps_array[stepcnt_x+1] && tr_y[0] >= try_steps_array[stepcnt_y] && tr_y[0] <try_steps_array[stepcnt_y+1]){
		   cluster_in_tr_cut_cnt_2D[stepcnt_x][stepcnt_y] ++;
		  }		
		}
	    }


	  Fill_Cluster_Histos(cluster_stack); // filling ADC and TDC histograms for only the PMTs in the cluster. This does not necessarily have to be done in a sperate function like this. 
	  
	  /*
	  if(!sh_flag) //filling different histos based on if there was a good hit in the shower or not
	    {
	      h_grinch_cluster_center_bad_sh_e ->Fill(clusterx,-clustery);
	      //Fill_Bad_Cluster_Histos(cluster_stack);
	      h_grinch_cluster_size_bad ->Fill(ClusterSize);
	    }
	  if(sh_flag)
	    {
	       h_grinch_cluster_center_good_sh_e ->Fill(clusterx,-clustery);
	       // Fill_Cluster_Histos(cluster_stack);
	       h_grinch_cluster_size_good ->Fill(ClusterSize);
	    }
	  */

	  if(ClusterSize == 3)
	    {
	      h_grinch_cluster_center_3 -> Fill(clusterx,-clustery);
	    }

          if(ClusterSize == 4)
	    {
	      h_grinch_cluster_center_4 -> Fill(clusterx,-clustery);
	    }

	  if(ClusterSize == 5)
	    {
	      h_grinch_cluster_center_5 -> Fill(clusterx,-clustery);
	    }

	  if(ClusterSize == 6)
	    {
	      h_grinch_cluster_center_6 -> Fill(clusterx,-clustery);
	    }
	  if(ClusterSize == 7)
	    {
	      h_grinch_cluster_center_7 -> Fill(clusterx,-clustery);
	    }

	  

	  Int_t vertspread = 0;
	  Int_t horizspread = 0; 

	  Find_Cluster_Extrema(cluster_stack,horizspread,vertspread); // function measures the height and width of the cluster
	  h_grinch_cluster_spread->Fill(horizspread,vertspread);
	  
          if(vertspread <=6 && horizspread <=5) //if the cluster is wider or taller than this, it is most likely multiple clusters. 
	    {
	      h_grinch_cluster_center_spreadcut ->Fill(clusterx,-clustery);
	      h_grinch_cluster_size_cuts ->Fill(ClusterSize);
	    }
	  
	  //cout<<"cluster is "<<horizspread<<" wide, "<< vertspread<<" tall"<<endl;


	  //cout<< "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " <<endl;
	}	
	
      //// Asking if the cluster had a good shower hit or if there was a good shower hit and no cluster
      if(sh_flag && g_cluster_flag)
	{
	  grinch_hit_sh_hit_cnt++;
	}
      if(sh_flag && !g_cluster_flag)
	{
	  sh_hit_no_grinch_hit_cnt++;
	}
      if(!sh_flag && g_cluster_flag)
	{
	  grinch_hit_no_sh_hit_cnt++;
	}
      if(!sh_flag && !g_cluster_flag)
	{
	  no_sh_no_grinch_cnt++;
	}
      
      
    }// end test if statment
       

  }//end event loop

  
  Double_t background_sum[510] = {0};
  Double_t background_per_bin[510];
  Double_t rate[510] = {0};
  Double_t histo_int;
  const Double_t ns = 0.000000001;

  Double_t pmt[510];

  Double_t bincontent;
  Int_t histo_entries;

  TF1* background_fit[510];
  Double_t background_from_fit[510];
 

  //This only works with all the TDC hits. Although it does work for giving the offset gaussian fit a place to start maybe? 

  for (Int_t i =0;i<510;i++)
    {
      for (Int_t bin = 500; bin<=700 ;bin++)
	{ 
	  bincontent =  h_grinch_le[i]->GetBinContent(bin);
	  background_sum[i]= background_sum[i] + bincontent;
	  histo_entries = h_grinch_le[i]->GetEntries();
	}
      background_fit[i] = new TF1(Form("background_fit_%d",i),"pol0",500,700);
      h_grinch_le[i] ->Fit(background_fit[i],"RQ0");
      background_from_fit[i] = background_fit[i] ->GetParameter(0);
      pmt[i] = i;
     
      background_per_bin[i] = background_sum[i]/200;

      // h_grinch_background_elem -> Fill(i,background_per_bin[i]);
      background_graph ->AddPoint(i, background_from_fit[i]);
      rate[i] = background_from_fit[i] / ( ns * max) /1000;
      h_grinch_rate_elem ->Fill(i, rate[i]);
      //cout<< "rate "<< i << ": "<<rate[i]<<endl;
    }
  rate_graph = new TGraph(510,pmt,rate);

  rate_graph->SetTitle(Form("GRINCH background rate run %d",nrun)); 
  rate_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
  rate_graph ->GetYaxis() -> SetTitle("Calculated Background rate (KHz)");
 

  

  

  TF1* fit[510];
  Double_t mean[510];
  Double_t mean_error[510];
  Double_t sigma[510];
  Double_t sigma_error[510];
  Double_t offset[510];
  Double_t offset_error[510];
  Double_t sigma_squared[510];
  Double_t sigma_squared_error[510];

  TF1* fit_cluster[510];
  Double_t mean_cluster[510];
  Double_t sigma_cluster[510];
  Double_t offset_cluster[510];

   
  

    
  for (Int_t i = 0;i<510; i++)
    {
      //fit[i] = new TF1(Form("fit_%d",i),"gaus",800,1000);     
      fit[i] = new TF1(Form("fit_%d",i),offset_gaus,800,1000,4);
      fit[i] ->SetParameters(2*background_per_bin[i],910,5,background_per_bin[i]);//need to do some stats on the histo to give these params?
      fit[i] ->SetParNames("Amp.","Mean","Sigma","Offset");
      TFitResultPtr r = h_grinch_le[i] ->Fit(fit[i],"RQ0");
      Int_t fitStatus = r;
      if (fitStatus != 0)
	{
	  cout<<"fit error for PMT "<< i <<" ,"<< " Status = "<<r<<endl;

	  mean[i] = -10;
      	  mean_error[i] = 0;
      	  sigma[i]= -10;
      	  sigma_error[i]= 0;
      	  offset[i]= -10;
      	  offset_error[i]=0;
	  // sigma_squared[i] = -10;
	  //sigma_squared_error[i] = 0;
	  continue;	  
	}

      mean[i] = fit[i] ->GetParameter(1);
      mean_error[i] = fit[i] ->GetParError(1);      
      sigma[i] = abs(fit[i] ->GetParameter(2));
      sigma_error[i] = fit[i] ->GetParError(2);
      offset[i] = fit[i] ->GetParameter(3);
      offset_error[i] = fit[i] ->GetParError(3);
      
      //sigma[i] = sqrt(sigma_squared[i]);      
      //sigma_error[i] = 0.5 * 1/sqrt(sigma_squared[i]) * sigma_squared_error[i];
      

      
      // if(mean[i] <0 || mean[i]>1200 || sigma[i]>300 || sigma_error[i] >100)
      // 	{
      // 	  cout<<"bad fit in PMT "<< i <<endl;
      // 	  mean[i] = -10;
      // 	  mean_error[i] = 0;
      // 	  sigma[i]= -10;
      // 	  sigma_error[i]= 0;
      // 	  offset[i]= -10;
      // 	  offset_error[i]=0;
      // 	}
      

      
      // //fit_cluster[i] = new TF1(Form("fit_cluster_%d",i),"gaus",875,925);
      // fit_cluster[i] = new TF1(Form("fit_cluster_%d",i),offset_gaus,875,925,4);
      // fit_cluster[i] ->SetParameters(20,900,10,1);
      // fit_cluster[i] ->SetParNames("Amp.","Mean","Sigma","Offset");
      // h_grinch_cluster_le[i] ->Fit(fit_cluster[i],"RQ0");
      // mean_cluster[i] = fit_cluster[i] ->GetParameter(1);
      // sigma_cluster[i] = abs(fit_cluster[i] ->GetParameter(2));
      // offset_cluster[i] = fit_cluster[i] ->GetParameter(3);
      
    }

  Double_t rate_from_offset_250=0;
  Double_t rate_from_offset_253=0;
  Double_t rate_from_offset_244=0;
  Double_t rate_from_offset_165=0;
  Double_t rate_from_offset_168=0;
  Double_t rate_from_offset_162=0;
  Double_t rate_from_offset_420=0;
  Double_t rate_from_offset_423=0;
  Double_t rate_from_offset_434=0;


  rate_from_offset_250 = offset[250] / (ns*max)/ 1000;
  rate_from_offset_253 = offset[253] / (ns*max)/ 1000;
  rate_from_offset_244 = offset[244] / (ns*max)/ 1000;
  rate_from_offset_165 = offset[165] / (ns*max)/ 1000;
  rate_from_offset_168 = offset[168] / (ns*max)/ 1000;
  rate_from_offset_162 = offset[162] / (ns*max)/ 1000;
  rate_from_offset_420 = offset[420] / (ns*max)/ 1000;
  rate_from_offset_423 = offset[423] / (ns*max)/ 1000;
  rate_from_offset_434 = offset[434] / (ns*max)/ 1000;


 



  cout<< "-------------------------------------------------------"<<endl;
  cout<< "rates from the offset fits for run "<<nrun<<endl;
  cout<<"PMT 162: "<<rate_from_offset_162 << " KHz,   PMT 165: "<<rate_from_offset_165 <<" KHz,   PMT 168: "<<rate_from_offset_168<<" KHz"<<endl;
  cout<<"PMT 244: " <<rate_from_offset_244<< " KHz,   PMT 250: "<<rate_from_offset_250 <<" KHz,   PMT 253: "<<rate_from_offset_253<<" KHz"<< endl;
  cout<<"PMT 434: " <<rate_from_offset_434<< " KHz,   PMT 420: "<<rate_from_offset_420 <<" KHz,   PMT 423: "<<rate_from_offset_423<<" KHz"<< endl;
  cout<<"--------------------------------------------------------"<<endl;



  



  Double_t zero[510] = {0};
  
  mean_graph = new TGraphErrors(510,pmt,mean,zero,mean_error);
  sigma_graph = new TGraphErrors(510,pmt,sigma,zero,sigma_error);
  offset_graph = new TGraphErrors(510,pmt,offset,zero,offset_error);
  //sigma_squared_graph = new TGraphErrors(510,pmt,sigma_squared,zero,sigma_squared_error);




  mean_graph->SetTitle(Form("GRINCH LE Mean run %d",nrun)); 
  mean_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
  mean_graph ->GetYaxis() -> SetTitle("Mean from Gausian fit");
  sigma_graph->SetTitle(Form("GRINCH LE Sigma run %d",nrun)); 
  sigma_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
  sigma_graph ->GetYaxis() -> SetTitle("Sigma from Gausian fit");
  offset_graph->SetTitle(Form("GRINCH LE Offset run %d",nrun)); 
  offset_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
  offset_graph ->GetYaxis() -> SetTitle("Offset from Gausian fit");

  // sigma_squared_graph->SetTitle(Form("GRINCH LE Sigma Squared run %d",nrun)); 
  //sigma_squared_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
  // sigma_squared_graph ->GetYaxis() -> SetTitle("Sigma Squared from Gausian fit");
   
 
  

  /*
  for (Int_t i=0 ;i<510; i++)
    {
      cout<<"PMT "<<i<< ": "<< "mean = " << mean[i] << ", mean_error = " << mean_error[i] << endl;
      cout<< " sigma = "<< sigma[i]<< ", sigma_error = "<< sigma_error[i]<< endl;
      cout<< " sigma squared = "<< sigma_squared[i]<< ", sigma_squared_error = "<< sigma_squared_error[i]<< endl;
      cout<< " offset = "<<offset[i]<<", offset_error = "<< offset_error[i]<< endl; 
      
      cout<<"-"<<endl;
    }
  */
  
  


  /*
  TF1* test_fit = new TF1("test_fit",offset_gaus,800,1000,4);
  test_fit ->SetParameters(20,900,10);
  test_fit ->SetParNames("Amp.","Mean","Sigma Squared","Offset");
  h_grinch_le[446] ->Fit("test_fit");
  */
  
  
  
  
  
  


  cout<< "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<" written to file "<<outputhist<<endl;
  cout<<"cluster analysis on "<<max<<" events"<<endl;
    cout<<clustercnt<<" clusters found"<<endl;
  cout<<"hit_sum = "<< hit_sum << endl;

  cout<<"good shower energy on "<< sh_ps_hit_cnt<<" events" <<endl;
  cout<<"good shower energy and NO good grinch cluster on "<< sh_hit_no_grinch_hit_cnt<<" events" <<endl;
  cout<<"good shower energy and good grinch cluster on "<< grinch_hit_sh_hit_cnt <<" events" <<endl;
  //cout<<"NO good shower energy and good grinch cluster on "<< grinch_hit_no_sh_hit_cnt <<" events" <<endl;
  //cout<<"NO good shower energy and NO grinch cluster on "<< no_sh_no_grinch_cnt <<" events" <<endl;

  //h_grinch_sh_hit_no_grinch_hit_cnt ->Fill(sh_hit_no_grinch_hit_cnt);
  //h_grinch_grinch_hit_sh_hit_cnt ->Fill(grinch_hit_sh_hit_cnt);
  // h_grinch_grinch_hit_no_sh_hit_cnt ->Fill(grinch_hit_no_sh_hit_cnt);
  //h_grinch_no_sh_no_grinch_cnt ->Fill(no_sh_no_grinch_cnt);
  // h_grinch_sh_ps_hit_cnt ->Fill(sh_ps_hit_cnt);

  Double_t grinch_hit_double = grinch_hit_sh_hit_cnt;
  Double_t sh_hit_double = sh_ps_hit_cnt;
  Double_t percent = 100 * grinch_hit_double/sh_hit_double; 

  Double_t cluster_in_tr_cut_cnt_double = cluster_in_tr_cut_cnt;
  Double_t tr_cut_cnt_double = tr_cut_cnt ; 
  Double_t tr_cut_percent = 100 * cluster_in_tr_cut_cnt_double/tr_cut_cnt_double;


  h_grinch_efficiency_percent->Fill(percent);

  cout<<"Effiency on Good Hits: "<< grinch_hit_sh_hit_cnt<< "/" << sh_ps_hit_cnt <<" = " << percent <<" %"<<endl; 

  cout<<"Effiency on the small track region we cut on for a test: " << cluster_in_tr_cut_cnt<<" / "<< tr_cut_cnt << " = "<< tr_cut_percent << " %"<<endl;

  Double_t cluster_in_tr_cut_cnt_double_2D[100][50] = {0};
  Double_t tr_cut_cnt_double_2D[100][50] = {0};
  Double_t tr_cut_percent_2D[100][50] = {0};
 
  for (Int_t stepcnt_x = 0; stepcnt_x < nsteps_trx; stepcnt_x ++)
    {
      for (Int_t stepcnt_y = 0; stepcnt_y < nsteps_try; stepcnt_y ++)
	{
	  Double_t percent = 0;
	  cluster_in_tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] = cluster_in_tr_cut_cnt_2D[stepcnt_x][stepcnt_y];
	  tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] = tr_cut_cnt_2D[stepcnt_x][stepcnt_y];
	  if(tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] == 0 && cluster_in_tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] == 0){
	    percent = 0;
	     h_grinch_eff_xy ->SetBinContent(stepcnt_y+1,stepcnt_x+1, 1);
	  }
	  else  if( tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] != 0 && cluster_in_tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] ==0)
	    {
	     percent = 0;
	      h_grinch_eff_xy ->SetBinContent(stepcnt_y+1,stepcnt_x+1, percent);
	    }
	  else{
	    percent = 100 * cluster_in_tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y]/tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y];
	     h_grinch_eff_xy ->SetBinContent(stepcnt_y+1,stepcnt_x+1, percent);
	  }	  
	  tr_cut_percent_2D[stepcnt_x][stepcnt_y] = percent;	 
	  cout<<"["<<trx_steps_array[stepcnt_x]<<"] ["<<try_steps_array[stepcnt_y]<<"]       "<< cluster_in_tr_cut_cnt_2D[stepcnt_x][stepcnt_y] <<" / "<< tr_cut_cnt_2D[stepcnt_x][stepcnt_y]  << " = "<< tr_cut_percent_2D[stepcnt_x][stepcnt_y]  << " %"<<endl;
	}		
    }
 
  h_grinch_eff_xy ->SetContour(256);

  TH2F* test2D = new TH2F("test2D",";horiz;vert;",  8,-0.2,0.2,28,-0.7,0.7);
  Int_t testvalue=1;

  for (Int_t testcnt = 1; testcnt <= 28; testcnt++)
    {
      for(Int_t testcnt2 = 1 ;testcnt2 <= 8 ;testcnt2 ++)
	{
	  test2D ->SetBinContent(testcnt2, testcnt, testvalue);
	  testvalue++;
	}
    }
  
  test2D ->SetContour(256);

  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  

  TFile hsimc(outputhist,"recreate");
  HList.Write();

}// end main

Double_t offset_gaus(Double_t *x, Double_t *par)
{
  Double_t height = par[0];
  Double_t mean = par[1];
  //Double_t sigma_squared = par[2];
  Double_t sigma = par[2];
  Double_t offset = par[3];

  Double_t fit = height * exp(-0.5 * pow((x[0] - mean)/sigma,2)) +offset;
  return fit;

}







////
//// See how many adjacent PMTs also had a hit
////
Int_t Sum_Adjacent_Hits(Int_t PMT, stack <Int_t>& neighbors_stack)
{
  //// I'm considering changing this function so that is also fills a stack with the values of it's neighbors. That way, when we are summing the PMTs in a cluster, we don't have to also search nighbors. But it does work fine just how I have it. 
  Int_t row = row_array[PMT];
  Int_t col = col_array[PMT];
  Int_t hit_sum = 0;
  stack <Int_t> searchstack;
  Int_t neighbor;
 
  if(hit_flag_array[PMT] == 0)
    {
      return 0;
    }
  
  //cout<<"Sum_Adjacent_Hits : "<<PMT<<endl;

  if (row%2 == 0)//short row of 8
    {
      if (col == 0)// on the left edge
	{
	  searchstack = short_row_on_left_search;
	  //cout<<"short row on left"<<endl;
	}
      else if (col == 7) // on the right edge
	{
	  searchstack = short_row_on_right_search;
	  //cout<<"short row on right"<<endl;
	}
      else
	{
	  //cout<<"short row in middle"<<endl;
	  searchstack = middle_search;
	}
    }
  else //long row of 9
    {
      if (col == 0 || col == 8)
	{
	  // don't sum anything. I don't want to cluster find on the long row edges. 
	  //cout<<"long row on edge"<<endl;
	  hit_sum=0;
	}
      else //not on an edge
	{
	  //cout<<"long row in the middle"<<endl;
	  searchstack = middle_search;
	}
    }

  //Print_Stack(searchstack);
  //cout<<"searchstack"<<endl;

  // stack <Int_t> neighborstack = Add_to_Stack(searchstack,PMT); //"PMT" to every element of the stack 
  // Print_Stack(neighborstack);
  // cout<<"neighborstack"<<endl;  
 
  neighbors_stack.push(PMT);
  while (!searchstack.empty())
    {
      neighbor = PMT-searchstack.top();
      hit_sum = hit_sum + hit_flag_array[neighbor];     
      if(hit_flag_array[neighbor] == 1)
	{
	  neighbors_stack.push(neighbor);
	}
    searchstack.pop();
    }
  
  //cout<< "hit_sum =" <<hit_sum <<endl;
  return hit_sum;

}//end Sum_Adjacent_Hits





void Fill_Search_Stacks()
{
  //// pre-filling stacks to make searching for neighbors and cluster finding cleaner in the other functions. 
  ////
  //// if [0] is the PMT we want to check the neighbors of, this is the PMT numbers of it's neighbors. 
  ////
  ////   [-9] [-8]
  //// [-1] [0] [+1]
  ////   [+8] [+9]
  ////
  ////
  //// So if we want to search around PMT 20, it would look like this:
  ////    [11] [12]
  ////  [19] [20] [21]
  ////    [28] [29]
  //// 
  //// If we are on the edge of the detector, we need to adjust slighty. 
  ////
  //// short_row_on_left : exclude [-1]
  //// short_row_on_right: exclude [+1]
  //// top_row_middle : exclue [-8], [-9]
  //// top_row_left : exclude [-1], [-8], [-9] this is only for PMT 0
  //// top_row_right: exclude [-9] ,[-8], [+1] this is only for PMT 7
  //// bottom_row_middle : exclude [+8], [+9] 
  //// bottom_row_left : exclude [-9], [-1], [+8], [+9] this is only for PMT 501
  //// bottom_row_right: exclude [-8], [+1], [+8], [+9] this is only for PMT 509
  ////

  Int_t search[6]= {-9,-8,9,8};
  for (Int_t counter = 0; counter < 4; counter ++)
    {
      short_row_on_left_search.push(search[counter]);
      short_row_on_right_search.push(search[counter]);
      middle_search.push(search[counter]);
    }
  short_row_on_left_search.push(1);
  middle_search.push(1);
  short_row_on_right_search.push(-1);
  middle_search.push(-1);
  
  long_row_on_left_search.push(-8);
  long_row_on_left_search.push(1);
  long_row_on_left_search.push(9);
  long_row_on_right_search.push(-9);
  long_row_on_right_search.push(-1);
  long_row_on_right_search.push(8);

  top_row_middle_search.push(-1);
  top_row_middle_search.push(1);
  top_row_middle_search.push(8);
  top_row_middle_search.push(9);

  top_row_left_search.push(8);
  top_row_left_search.push(1);
  top_row_left_search.push(9);

  top_row_right_search.push(-1);
  top_row_right_search.push(8);
  top_row_right_search.push(9);

  bottom_row_middle_search.push(-1);
  bottom_row_middle_search.push(1);
  bottom_row_middle_search.push(-8);
  bottom_row_middle_search.push(-9);  

  bottom_row_left_search.push(-8);
  bottom_row_left_search.push(1);

  bottom_row_right_search.push(-1);
  bottom_row_right_search.push(-9); 
}// end Fill_Search_Stacks


void Calculate_PMT_Coord(double& horiz, double& vert, Int_t PMT) // this is probably uncessary as these should already be in the database
{
  Int_t row = row_array[PMT];
  Int_t col = col_array[PMT]; 

  double temp_horiz = 0;
  
  if (row%2 == 0)//on a short row of 8
    {
      temp_horiz = col +0.5;
    }
  else //on a long row of 9
    {
      temp_horiz = col;
    }
  
  horiz = temp_horiz;
  vert = row;
  return;
}// end Calculate_PMT_Coord


void Calculate_Cluster_Center(stack <Int_t> inputstack, double& horiz, double& vert)
{//// Takes the average of the horizontal and vertical positions of the PMTs in the cluster. 
 //// The result is returned via reference parameter. 

  Double_t horiz_sum=0;
  Double_t vert_sum=0;
  Double_t counter=0;
  Double_t temp_horiz = 0;
  Double_t temp_vert = 0;
  Int_t PMT=0;

  while (!inputstack.empty())
    {
      PMT = inputstack.top();
      Calculate_PMT_Coord(temp_horiz,temp_vert,PMT);
      horiz_sum = horiz_sum + temp_horiz;
      vert_sum = vert_sum + temp_vert;
      counter = counter +1;
      inputstack.pop();
    }
  
  horiz = horiz_sum/counter;
  vert = vert_sum/counter;  
} //end calculate cluster center


void Find_Cluster_Extrema(stack <Int_t> inputstack, Int_t& horizspread, Int_t& vertspread)
{ //// finds the height and width of the cluster and returns them via reference parameter
  Int_t maxcol;
  Int_t mincol;
  Int_t maxrow;
  Int_t minrow;
  Int_t temprow = 0;
  Int_t tempcol = 0;
  Int_t pmt = 0;

  pmt = inputstack.top();
  maxrow = row_array[pmt];
  minrow = row_array[pmt];
  maxcol = col_array[pmt];
  mincol = col_array[pmt];
  inputstack.pop();

  while(!inputstack.empty())
    {
      pmt = inputstack.top();
      inputstack.pop();
      temprow = row_array[pmt];
      tempcol = col_array[pmt];      
      if(temprow > maxrow){maxrow = temprow;}
      if(temprow < minrow){minrow = temprow;}
      if(tempcol > maxcol){maxcol = tempcol;}
      if(tempcol < mincol){mincol = tempcol;}
    }
  //cout<<"minrow = "<<minrow<<endl;
  //cout <<"maxrow = "<<maxrow<<endl;
  //cout <<"mincol = "<<mincol<<endl;
  //cout <<"maxcol = "<<maxcol<<endl;
  vertspread = maxrow - minrow +1;
  horizspread = maxcol - mincol +1;  
    
  return;   
}// end Find_Cluster_Extrema



void Fill_Cluster_Histos(stack <Int_t> inputstack)
{ //// Fills the hitograms for the cluser in this function to keep the main loop clearer. 
  Int_t PMT;
  Int_t root_id;
  Int_t ADC_chan;
  Int_t ADC_root_id;
  Double_t grinch_minus_hodo_time;

  while (!inputstack.empty())
    {
      PMT = inputstack.top();
      inputstack.pop();
      root_id = root_index_array[PMT];//get the index that the branches need in order to look at the specific PMT for this event
      h_grinch_cluster_le_elem -> Fill(PMT,tdcGLe[root_id] - hodo_tmean);
      h_grinch_cluster_te_elem -> Fill(PMT,tdcGTe[root_id]);
      h_grinch_cluster_tot_elem -> Fill(PMT,tdcGTot[root_id]);
      h_grinch_cluster_mult_elem -> Fill(PMT,tdcGMult[root_id]);
      h_grinch_cluster_tot[PMT]->Fill(tdcGTot[root_id]);
      h_grinch_cluster_le[PMT]->Fill(tdcGLe[root_id]-hodo_tmean);
      h_grinch_cluster_le_all -> Fill(tdcGLe[root_id]-hodo_tmean);
      h_grinch_cluster_tot_all ->Fill(tdcGTot[root_id]);
      h_grinch_cluster_elem ->Fill(PMT);

      grinch_minus_hodo_time = tdcGLe[root_id] - hodo_tmean;
      h_grinch_cluster_le_elem_hodo ->Fill(PMT, grinch_minus_hodo_time);
      

      ////see if the PMT has a corresponding ADC channel 
      ADC_chan = tdc_to_adc_map_array[PMT];
      if(ADC_chan !=-1)//the PMT has a corresponding ADC channel. 
	{
	  ADC_root_id = adc_root_index_array[ADC_chan];
	  h_grinch_cluster_atime_elem ->Fill(ADC_chan, adcGAtime[ADC_root_id]);
	  h_grinch_cluster_amp_elem ->Fill(ADC_chan, adcGAmp[ADC_root_id]);
	  h_grinch_cluster_amp[ADC_chan] ->Fill(adcGAmp[ADC_root_id]);
	  h_grinch_cluster_int_elem ->Fill(ADC_chan, grinch_adc[ADC_root_id]);
	  h_grinch_cluster_adcMult_elem ->Fill(ADC_chan,adcGMult[ADC_root_id]); 
	  h_grinch_cluster_adc_tdc ->Fill(tdcGTot[root_id],adcGAmp[ADC_root_id]); 
	  h_grinch_cluster_adcInt_tdc ->Fill(tdcGTot[root_id],grinch_adc[ADC_root_id]); 
	  h_grinch_cluster_int_amp ->Fill(adcGAmp[ADC_root_id],grinch_adc[ADC_root_id]);//just testing the amp and integral are proportional to eachother 
	  if(ADC_chan == 23)
	    { ////looking at one specific ADC channel for a sanity check. 
	      h_grinch_cluster_adc_tdc_test ->Fill(tdcGTot[root_id],adcGAmp[ADC_root_id]);
	    }
	}
    }
}////end Fill_Cluster_Histos


void Fill_Bad_Cluster_Histos(stack <Int_t> inputstack)
{
  //// "Bad" here just means that the shower cut was not passed. I was curious what clusters that did not pass the shower counter looked like. Answer: they look the same to me. 
  Int_t PMT;
  Int_t root_id;
  Int_t ADC_chan;
  Int_t ADC_root_id;
  while (!inputstack.empty())
    {
      PMT = inputstack.top();
      inputstack.pop();
      root_id = root_index_array[PMT];//get the index that the branches need in order to look at the specific PMT for this event
      h_grinch_cluster_bad_le_elem -> Fill(PMT,tdcGLe[root_id]);
      h_grinch_cluster_bad_te_elem -> Fill(PMT,tdcGTe[root_id]);
      h_grinch_cluster_bad_tot_elem -> Fill(PMT,tdcGTot[root_id]);
      h_grinch_cluster_bad_mult_elem -> Fill(PMT,tdcGMult[root_id]);

      //see if the PMT has a corresponding ADC channel 
      ADC_chan = tdc_to_adc_map_array[PMT];
      if(ADC_chan !=-1)//the PMT has a corresponding ADC channel. 
	{
	  ADC_root_id = adc_root_index_array[ADC_chan];
	  h_grinch_cluster_bad_atime_elem ->Fill(ADC_chan, adcGAtime[ADC_root_id]);
	  h_grinch_cluster_bad_amp_elem ->Fill(ADC_chan, adcGAmp[ADC_root_id]);
	  h_grinch_cluster_bad_int_elem ->Fill(ADC_chan, grinch_adc[ADC_root_id]);
	  h_grinch_cluster_bad_adcMult_elem ->Fill(ADC_chan,adcGMult[ADC_root_id]); 
	  h_grinch_cluster_bad_adc_tdc ->Fill(tdcGTot[root_id],adcGAmp[ADC_root_id]); 
	  h_grinch_cluster_bad_adcInt_tdc ->Fill(tdcGTot[root_id],grinch_adc[ADC_root_id]); 
	  h_grinch_cluster_bad_int_amp ->Fill(adcGAmp[ADC_root_id],grinch_adc[ADC_root_id]);//just testing the amp and integral are proportional to eachother 
	  if(ADC_chan == 23)
	    {
	      h_grinch_cluster_bad_adc_tdc_test ->Fill(tdcGTot[root_id],adcGAmp[ADC_root_id]);
	    }
	}
    }
}//end Fill_Bad_Cluster_Histos


void make_tdc_to_adc_map_array() //get the ADC chan number from the PMT number (if there is one)
// note for eric: 
// this is also something that should maybe go in the database. We probably aren't going to be changing the location of the arrays anytime soon. 
{
  //Fill the array so that it is -1 where the PMT index does not have a corresponding TDC channel, and the ADC channel number where it does. 

  for (Int_t i = 0; i<510; i ++)//initialize all to -1
    {           
      tdc_to_adc_map_array[i]= -1;
    }
  
  Int_t ADC_cnt = 0; //ADCs are on PMTs 192 to 255 as of 2/1/2022. 255 is the LED trigger. Leaving that out for now. 
  for (Int_t PMT = 192; PMT <= 254; PMT ++)
    {
      tdc_to_adc_map_array[PMT] = ADC_cnt;
      ADC_cnt++;
    } 
}//end make_tdc_to_adc_map_array


void make_row_col_array()// making it easy to convert from PMT num to row, col
// this is also probably something that is in the database or should go in the database
{
  Int_t pmt = 0;
  Int_t N_COL;
  for (Int_t i = 0;i<N_ROW;i++)
    {
      N_COL = 8 + i%2;
      for (Int_t j=0 ; j< N_COL ; j++)
	{
	  row_array[pmt] = i;
	  col_array[pmt] = j;
	  pmt ++;
	}
    }
  /*
    for (Int_t r = 0; r< 510; r++)
    {
    cout<< "PMT "<<r<<": ROW "<< row_array[r] <<", COL"<< col_array[r]<<endl;
    }
  */
}// end make_row_col_array


Bool_t Is_NOT_In_Stack(stack <Int_t> checkstack, Int_t PMT) //simple function to loop through and see if a number is not in the stack. 
{
  while (!checkstack.empty())
    {
      if (PMT == checkstack.top())
	{
	  return kFALSE; // the PMT is already in the stack
	}
      checkstack.pop();
    }
  return kTRUE; // the PMT is not already in the stack
}// end Is_NOT_In_Stack


Int_t Stack_Size(stack <Int_t> inputstack)
{ //// simple function to return the number of elements in a stack
  Int_t counter = 0;
  while(!inputstack.empty())
    {
      counter ++;
      inputstack.pop();
    }
  return counter;
}

void Print_Stack(stack <Int_t> inputstack)
{ //// function to print out the contents of a stack for debugging purposes 
  // cout<< "~~~~~~~"<<endl;
  cout<<endl;
  while (!inputstack.empty())
    {
      cout<<"| " <<inputstack.top()<< " |" <<endl;
      inputstack.pop();
    }
  
  cout<< "|_____|"<<endl;
}// end Print_Stack


stack <Int_t> Add_to_Stack(stack <Int_t> inputstack, Int_t PMT) 
//// This is a function to help fill up search stacks
//// The integer "PMT" is added to every element of the stack
{
  stack <Int_t> tempstack;
  Int_t result = 0;
  while (!inputstack.empty())
    {
      result = PMT + inputstack.top();
      tempstack.push(result);
      inputstack.pop();
    }
  return tempstack;
} //end Add_to_Stack


stack <Int_t> Clear_Stack()
{ // simple function to return an empty stack
  stack <Int_t> clearstack; //default constructor is an empty stack
  return clearstack;
}// end Clear_Stack




 // void Draw_TBox( Double_t tr_cut_percent_2D[100][50] )
 // {  
 //   TCanvas* C = new TCanvas;
 //   TBox* eff_box[100][50];
   
 //   h_grinch_cluster_track_xy ->Draw("colz");

 //   for (Int_t stepcnt_x = 0; stepcnt_x < nsteps_trx; stepcnt_x ++)
 //    {
 //      for (Int_t stepcnt_y = 0; stepcnt_y < nsteps_try; stepcnt_y ++)
 // 	{
 // 	  eff_box[stepcnt_x][stepcnt_y] = new TBox(try_steps_array[stepcnt_y] ,trx_steps_array[stepcnt_x] ,try_steps_array[stepcnt_y +1] ,trx_steps_array[stepcnt_x +1] );
 // 	  eff_box[stepcnt_x][stepcnt_y]  ->SetFillColorAlpha(kBlue+stepcnt_y +stepcnt_x,0.35);
 // 	  eff_box[stepcnt_x][stepcnt_y] ->Draw("same");	
 // 	}
 //    }
   
 
   
 //   TBox* testbox = new TBox(0,0,0.5,0.5);
 //   testbox ->SetLineColor(1);
 //   testbox ->SetLineWidth(5);
 //   testbox -> SetFillColorAlpha(kBlue,0.35);
 //   // testbox ->Draw("same");  
 // }


void Palette1()
{
  Int_t n =256;
  static Int_t colors[256];
  static Bool_t initialized = kFALSE;
  
  Double_t Red[5] =  {0.00, 0.00, 1.00, 1.00, 1.00};
  Double_t Green[5] =  {0.00, 0.00, 0.00, 1.00, 1.00}; 
  Double_t Blue[5] = {0.00,1.00, 0.00, 0.00, 1.00};
  Double_t Length[5] = {0.00,0.8,0.9, 0.95, 1.00};


  if (!initialized)
    {
      Int_t FI = TColor::CreateGradientColorTable(5,Length,Red,Green,Blue,n);
      for (int i = 0;i<n;i++) colors [i] = FI+i;
      initialized = kTRUE;
      return;
    }
  gStyle ->SetPalette(n,colors);
}
