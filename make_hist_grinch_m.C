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

/////////////////////////////////////////////////////
//// make_hist_grinch_m.C
//// Maria Satik
//// January 2022
//// msatnik@email.wm.edu
////
//// This program processes data for the grinch from the Generic Detector class and performs basic cluster finding on it. 
//// A type of dynamic list, called a stack, is mainly utilized for this.
//// 
////
///////////////////////////////////////////////////


using namespace std;


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
stack <Int_t> cluster_members_stack;

const Int_t N_ROW=60;

Int_t row_array[510];
Int_t col_array[510];
Int_t tdc_to_adc_map_array[510]; // Needs to be updated when ADCs are moved to different Nino cards. 

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
Int_t GrinchADCNum;
Double_t adcGID[1000]; // ADC Channel Numbern
Double_t adcGAtime[1000]; // The time the signal arrives in the ADC window
Double_t adcGAmp[1000]; // Amplitude of the ADC signal
Double_t grinch_adc[1000];// Integral of the ADC signal
Double_t adcGMult[1000]; // ADC multiplity: how many times the signal went above threshold in the time window
UInt_t fTrigBits; // Should be "1" if it is a bbcal trigger, 16 if it is a GRINCH LED event (I think)
UInt_t fEvtNum; // the event number
 
Int_t hit_flag_array[510]; // array to help keep track of which PMTs were hit during an event
Int_t mult_array[510]; // array for the multiplicities of each PMT during an event
Int_t sum_array[510]; // array to help with cluster finding: records how many adjacent PMTs have a hit for each PMT
Int_t root_index_array[510]; // An array to help us convert betweeen the PMT number and the index that root uses for an event
Int_t adc_root_index_array[510]; // same as above but for ADC signals 

//histograms
TH1F* h_ps_e = new TH1F("h_ps_e"," ; Total PreShower E ",300,0.0,2.0);
TH1F* h_TBB_e = new TH1F("h_TBB_e"," ; Total PRe+Shower E ",300,0.0,5.0);
TH1F* h_ratio_e = new TH1F("h_ratio_e"," ; Total PRe+Shower E/ track norm ",300,0.0,2.5);
TH2F* h_grinch_atime_elem = new TH2F("h_grinch_atime_elem"," ; GRINCH ADC elemID ; GRINCH ADC time (ns) ",63,0,63,500,0,2000);
TH2F* h_grinch_atime_elem_multcut = new TH2F("h_grinch_atime_elem_multcut"," ; GRINCH ADC elemID ; GRINCH ADC time (ns) multcut ",63,0,63,500,0,2000); 
TH2F* h_grinch_atime_elem_SHcut  = new TH2F("h_grinch_atime_elem_SHcut"," ; GRINCH ADC elemID ; GRINCH ADC time (ns) Shower Cut ",63,0,63,500,0,2000);
TH2F* h_grinch_amp_elem = new TH2F("h_grinch_amp_elem"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) ",63,0,63,500,0,500);    
TH2F* h_grinch_amp_elem_multcut = new TH2F("h_grinch_amp_elem_multcut"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) multcut ",63,0,63,600,0,600); 
TH2F* h_grinch_amp_elem_SHcut = new TH2F("h_grinch_amp_elem_SHcut"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) Shower Cut ",63,0,63,600,0,600); 
TH2F* h_grinch_amp_elem_multcut_timecut = new TH2F("h_grinch_amp_elem_multcut_timecut"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) multcut timecut ",63,0,63,600,0,600);
TH2F* h_grinch_adcMult_elem = new TH2F("h_grinch_adcMult_elem"," ; GRINCH ADC elem id ; GRINCH ADC Mult. ",63,0,63,10,0,10);
TH2F* h_grinch_adcMult_elem_SHcut = new TH2F("h_grinch_adcMult_elem_SHcut"," ; GRINCH ADC elem id ; GRINCH ADC Mult. Shower Cut ",63,0,63,10,0,10);
TH2F* h_grinch_le_elem = new TH2F("h_grinch_le_elem"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,2700,-100,2600);
TH2F* h_grinch_le_elem_multcut = new TH2F("h_grinch_le_elem_multcut"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) multcut",510,0,510,2700,-100,2600);
TH2F* h_grinch_le_elem_SHcut = new TH2F("h_grinch_le_elem_SHcut"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) Shower Cut",510,0,510,2700,-100,2600);
TH2F* h_grinch_te_elem_multcut = new TH2F("h_grinch_te_elem_multcut"," ; GRINCH TDC elemID ; GRINCH TDC TE (ns) multcut",510,0,510,2700,-100,2600);
TH2F* h_grinch_te_elem = new TH2F("h_grinch_te_elem"," ; GRINCH TDC elemID ; GRINCH TDC TE (ns)",510,0,510,2700,-100,2600);
TH2F* h_grinch_te_elem_SHcut = new TH2F("h_grinch_te_elem_SHcut"," ; GRINCH TDC elemID ; GRINCH TDC TE (ns) Shower Cut",510,0,510,2700,-100,2600);
TH1F* h_grinch_ghits = new TH1F("h_grinch_ghits"," ; GRINCH Good hits ; ",10,0,10);
TH1F* h_grinch_elem = new TH1F("h_grinch_elem"," ; GRINCH TDC elemID ; ",510,0,510);
TH1F* h_grinch_tot_all = new TH1F("h_grinch_tot_all"," ; GRINCH TDC Tot (no cut) ; ",50,0,50);
TH1F* h_grinch_tot_all_lecut = new TH1F("h_grinch_tot_all_lecut"," ; GRINCH TDC Tot (TDC LE cut) ; ",50,0,50);
TH2F* h_grinch_tot_elem = new TH2F("h_grinch_tot_elem"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) ",510,0,510,110,-10,100);
TH2F* h_grinch_tot_elem_multcut = new TH2F("h_grinch_tot_elem_multcut"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) multcut ",510,0,510,110,-10,100);
TH2F* h_grinch_tot_elem_SHcut = new TH2F("h_grinch_tot_elem_SHcut"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) Shower Cut ",510,0,510,110,-10,100);
TH2F* h_grinch_mult_elem = new TH2F("h_grinch_mult_elem"," ; GRINCH TDC elemID ; GRINCH TD.qC Mult ",510,0,510,15,0,15);
TH2F* h_grinch_mult_elem_SHcut = new TH2F("h_grinch_mult_elem_SHcut"," ; GRINCH TDC elemID ; GRINCH TDC Mult Shower Cut ",510,0,510,15,0,15);
TH2F* h_grinch_adc_chan_vs_num = new TH2F("h_grinch_adc_chan_vs_num","; GRINCH ADC elemID ; GRINCH ADC channel ",64,0,63,510,-10,500);//maria
TH2F* h_grinch_adc_chan_vs_num_SHcut = new TH2F("h_grinch_adc_chan_vs_num_SHcut","; GRINCH ADC elemID ; GRINCH ADC channel Shower Cut",64,0,63,510,-10,500);//maria
TH2F* h_grinch_adc_chan_vs_num_multcut = new TH2F("h_grinch_adc_chan_vs_num_multcut","; GRINCH ADC elemID ; GRINCH ADC channel multcut ",64,0,64,510,-10,500);//ma1ria
TH1F* h_EvtHdr_TrigBits = new TH1F("h_EvtHdr_TrigBits", "; fTrigBits;",35,0,35);

TH1F* h_grinch_cluster_size = new TH1F("h_grinch_cluster_size" ,"; Cluster Size ",12,3,15);
TH1F* h_grinch_cluster_size_good = new TH1F("h_grinch_cluster_size_good" ,"; Cluster Size (shower cut) ",12,3,15);
TH1F* h_grinch_cluster_size_bad = new TH1F("h_grinch_cluster_size_bad" ,"; Cluster Size (outside shower cut)",12,3,15);
TH2F* h_grinch_cluster_center = new TH2F("h_grinch_cluster_center", "; horizontal ; vertical ;",16,0,8,60,0,60);
TH2F* h_grinch_cluster_center_display = new TH2F("h_grinch_cluster_center_display", "; horizontal ; vertical ;",16,0,8,60,-60,0);
TH2F* h_grinch_cluster_center_spreadcut = new TH2F("h_grinch_cluster_center_spreadcut",";horizontal; vertical;",16,0,8,60,-60,0);
TH2F* test_coord_histo = new TH2F(" test_coord_histo", "; horizontal ; vertical ;",18,0,9,60,0,60);
TH1F* h_grinch_pmt_of_interest = new TH1F("h_grinch_pmt_of_interest",";PMT ;",510,0,510);
TH1F* h_grinch_pmt_good_hit = new TH1F("h_grinch_pmt_good_hit",";PMT ;",510,0,510);
TH2F* h_grinch_cluster_spread = new TH2F("h_grinch_cluster_spread",";cluster width;  cluster height;",8,1,9,59,1,60);
TH1F* h_grinch_cluster_size_cuts = new TH1F("h_grinch_cluster_size_cuts" ,"; Cluster Size with cuts",12,3,15);
TH2F* h_grinch_cluster_center_bad_sh_e  = new TH2F("h_grinch_cluster_center_bad_sh_e ", "; horizontal ; vertical ;",16,0,8,60,-60,0);
TH2F* h_grinch_cluster_center_good_sh_e  = new TH2F("h_grinch_cluster_center_good_sh_e ", "; horizontal ; vertical ;",16,0,8,60,-60,0);


TH2F* h_grinch_cluster_le_elem = new TH2F("h_grinch_cluster_le_elem","; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,600,600,1200);
TH2F* h_grinch_cluster_te_elem = new TH2F("h_grinch_cluster_te_elem"," ; GRINCH TDC elemID ; GRINCH TDC TE (ns)",510,0,510,600,600,1200);
TH2F* h_grinch_cluster_tot_elem = new TH2F("h_grinch_cluster_tot_elem"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) ",510,0,510,110,-10,100);
TH2F* h_grinch_cluster_mult_elem = new TH2F("h_grinch_cluster_mult_elem"," ; GRINCH TDC elemID ; GRINCH TDC Mult ",510,0,510,10,0,10);
TH2F* h_grinch_cluster_adcMult_elem = new TH2F("h_grinch_cluster_adcMult_elem"," ; GRINCH ADC elem id ; GRINCH ADC Mult. ",63,0,63,10,0,10);
TH2F* h_grinch_cluster_atime_elem = new TH2F("h_grinch_cluster_atime_elem"," ; GRINCH ADC elemID ; GRINCH ADC time (ns) ",63,0,63,500,0,2000);
TH2F* h_grinch_cluster_amp_elem = new TH2F("h_grinch_cluster_amp_elem"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) ",63,0,63,500,0,500);    
TH2F* h_grinch_cluster_int_elem = new TH2F("h_grinch_cluster_int_elem","; GRINCH ADC elemID ; GRINCH ADC channel ",64,0,63,510,-10,500);//maria
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

TH1F* h_grinch_le[511];
TH1F* h_grinch_tot[511];
TH1F* h_grinch_tot_lecut[511];
TH1F* h_grinch_cluster_tot[511];
TH1F* h_grinch_adc[64];//maria 
TH1F* h_grinch_amp_p[64];//maria 
TH1F* h_grinch_cluster_amp_p[64];//maria 

Int_t sh_ps_hit_cnt = 0; //counter for the number of events that pass the shower cut
Int_t sh_hit_no_grinch_hit_cnt = 0; //counter for the number of events that pass the shower cut but do not have a grinch cluster
Int_t grinch_hit_sh_hit_cnt = 0; // counter for the number of event that pass the shower cut and have a grinch cluster
Int_t grinch_hit_no_sh_hit_cnt = 0;// counter for the number of events that have a grinch cluster but do not pass the shower cut
Int_t no_sh_no_grinch_cnt = 0;// counter for the number of events that do not pass the shower cut and do not have a grinch cluster

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
Int_t Sum_Adjacent_Hits(Int_t PMT);
stack <Int_t> Add_to_Stack(stack <Int_t> inputstack, Int_t PMT);
Int_t Calculate_Cluster_Size( stack <Int_t> inputstack);
void Calculate_PMT_Coord(double& horiz, double& vert, Int_t PMT); // May be unecessary: already in database
void Calculate_Cluster_Center(stack <Int_t> inputstack, double& horiz, double& vert);
void Find_Cluster_Extrema(stack <Int_t> inputstack, Int_t& horizspread, Int_t& vertspread);
void Fill_Cluster_Histos(stack <Int_t> inputstack);
void Fill_Bad_Cluster_Histos(stack <Int_t> inputstack);
void make_tdc_to_adc_map_array(); // Probably uncessary: should add to database instead



void make_hist_grinch_m(TString basename="",Int_t nrun=2043,TString configfilename="run_list_grinch.txt", Bool_t show_track = kFALSE){ //MAIN
  if (basename=="") {
    cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
    cin >> basename;
  }
  TString fullname = Form("%s_%d",basename.Data(),nrun);
  gStyle->SetPalette(1,0);
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
  fchain->SetBranchAddress("Ndata.bb.grinch_tdc.tdcelemID",&GrinchNum) ;//number of PMTs with a signal in an event (I think)   
  fchain->SetBranchAddress("bb.grinch_tdc.tdcelemID",&tdcGID) ; // The PMT number 
  fchain->SetBranchAddress("bb.grinch_tdc.tdc",&tdcGLe) ; // leading edge
  fchain -> SetBranchAddress("bb.grinch_tdc.tdc_te",&tdcGTe) ; // trailing edge
  fchain->SetBranchAddress("bb.grinch_tdc.tdc_mult",&tdcGMult) ; //tdc multiplicity
  fchain->SetBranchAddress("bb.grinch_tdc.tdc_tot",&tdcGTot) ; // tdc time over threshold (te-le)
     
  fchain->SetBranchAddress("Ndata.bb.grinch_adc.adcelemID",&GrinchADCNum) ; //number of PMTs with ADC channels with a signal in an event. 
  fchain->SetBranchAddress("bb.grinch_adc.adcelemID",&adcGID) ; // The ADC channel 
  fchain->SetBranchAddress("bb.grinch_adc.a_time",&adcGAtime) ; // the time the first(?)ADC signal went over threshold in the window  
  fchain->SetBranchAddress("bb.grinch_adc.a_amp_p",&adcGAmp) ; // pedestal-subtracted amplitude
  fchain->SetBranchAddress("bb.grinch_adc.a_mult",&adcGMult) ;// ADC multiplicity 
  fchain->SetBranchAddress("bb.grinch_adc.a_p",&grinch_adc) ; // ADC Integral (I think)

  fchain->SetBranchAddress("fEvtHdr.fTrigBits",&fTrigBits); // The type of trigger the event was. (1 is the bbcal trig, 16 is the grinch LED) (May not be working)
  fchain->SetBranchAddress("fEvtHdr.fEvtNum", &fEvtNum); // The event number

 
  //Making histos for the individual channels 
  for (Int_t ig=0;ig<511;ig++) {
    h_grinch_le[ig] = new TH1F(Form("h_grinch_le_%d",ig),Form(" ; GRINCH TDC LE (ns) PMT %d  ; ",ig),25,875,925);
    HList.Add(h_grinch_le[ig]);
    h_grinch_tot[ig] = new TH1F(Form("h_grinch_tot_%d",ig),Form(" ; GRINCH TDC Tot PMT %d (no cut) ; ",ig),50,0,50);
    HList.Add(h_grinch_tot[ig]);
    h_grinch_tot_lecut[ig] = new TH1F(Form("h_grinch_tot_lecut_%d",ig),Form(" ; GRINCH TDC Tot PMT %d (LE cut) ; ",ig),50,0,50);
    HList.Add(h_grinch_tot_lecut[ig]);   
    h_grinch_cluster_tot[ig] = new TH1F(Form("h_grinch_cluster_tot_%d",ig),Form(" ; GRINCH TDC ToT CLUSTER PMT %d ; ",ig),50,0,50);
    HList.Add(h_grinch_cluster_tot[ig]);
  }
 
    
  for (Int_t i=0;i<64;i++){
    h_grinch_adc[i] = new TH1F(Form("h_grinch_adc_%d",i),Form(" ; GRINCH ADC Channels ADC %d ; ",i), 201, -0.5, 200.5); 
    HList.Add(h_grinch_adc[i]);
    h_grinch_amp_p[i] = new TH1F(Form("h_grinch_amp_p_%d",i),Form(" ; GRINCH ADC Amp (pedsub) %d ; ",i), 501, -0.5, 500.5); 
    HList.Add(h_grinch_amp_p[i]);
    h_grinch_cluster_amp_p[i] = new TH1F(Form("h_grinch_cluster_amp_p_%d",i),Form(" ; ADC Amp CLUSTER %d ; ",i), 501, -0.5, 500.5); 
    HList.Add(h_grinch_cluster_amp_p[i]);
  }

  //// functions that initialize a few things 
  Fill_Search_Stacks(); //sets up cluster finding stacks
  make_row_col_array(); //sets up an array to quickly get the row/col from PMT number
  make_tdc_to_adc_map_array();//sets up arrays to quickly convert between TDC number and ADC number 
  
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

  HList.Add(h_grinch_atime_elem);
  HList.Add(h_grinch_atime_elem_multcut);
  HList.Add(h_grinch_atime_elem_SHcut);

  HList.Add(h_grinch_amp_elem);     
  HList.Add(h_grinch_amp_elem_multcut);
  HList.Add(h_grinch_amp_elem_SHcut);

  HList.Add(h_grinch_amp_elem_multcut_timecut);
  HList.Add(h_grinch_adcMult_elem);
  HList.Add(h_grinch_adcMult_elem_SHcut);
  HList.Add(h_grinch_le_elem);
  HList.Add(h_grinch_le_elem_multcut);
  HList.Add(h_grinch_le_elem_SHcut);
  HList.Add(h_grinch_te_elem_multcut);
  HList.Add(h_grinch_te_elem);
  HList.Add(h_grinch_te_elem_SHcut);

  HList.Add(h_grinch_ghits);
  HList.Add(h_grinch_elem);
  HList.Add(h_grinch_tot_all);
  HList.Add(h_grinch_tot_all_lecut);
  HList.Add(h_grinch_tot_elem);
  HList.Add(h_grinch_tot_elem_multcut); 
  HList.Add(h_grinch_tot_elem_SHcut); 
  HList.Add(h_grinch_mult_elem);

  HList.Add(h_grinch_mult_elem_SHcut);
  HList.Add(h_grinch_adc_chan_vs_num);
  HList.Add(h_grinch_adc_chan_vs_num_SHcut);
  HList.Add(h_grinch_adc_chan_vs_num_multcut);
  HList.Add(h_EvtHdr_TrigBits);

  HList.Add(h_grinch_cluster_size);
  HList.Add(h_grinch_cluster_size_good);
  HList.Add(h_grinch_cluster_size_bad);
  HList.Add(h_grinch_cluster_center);
  HList.Add(h_grinch_cluster_center_display);
  HList.Add(test_coord_histo);
  HList.Add(h_grinch_pmt_of_interest);
  HList.Add(h_grinch_pmt_good_hit);
  HList.Add(h_grinch_cluster_spread);
  HList.Add(h_grinch_cluster_center_spreadcut);
  HList.Add(h_grinch_cluster_size_cuts);
  HList.Add(h_grinch_cluster_center_bad_sh_e);
  HList.Add(h_grinch_cluster_center_good_sh_e);
  HList.Add(h_grinch_cluster_adc_tdc_test);

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


  Long64_t nentries = fchain->GetEntries();
  cout<<"nentries: "<<nentries<<endl;

  Int_t temp_cnt=0;

  Int_t zero_count = 0; 
  Int_t clustercnt = 0;
  Int_t max = 20000; // make this lower if you don't want to analyze all of the entries

  
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
    if (i%1000==0) cout << " Entry = " << i << endl;
    
    Double_t tot_e = ps_e+sh_e;
    Double_t rat = tot_e/pmom;
    h_ps_e->Fill(ps_e);
    h_TBB_e->Fill(tot_e);
    h_ratio_e->Fill(rat);
    h_EvtHdr_TrigBits ->Fill(fTrigBits);
    Bool_t sh_flag = kFALSE; //flag for when an event passes the shower cut
    Bool_t g_cluster_flag =kFALSE;//flag for when an event has a grinch cluster

    //// initialize some arrays
    for(Int_t n = 0; n< 64; n++) //initilalize root index array
      {
	adc_root_index_array[n]=0;
      }

    
    if (ps_e > 0.2 && abs(rat - 1.15) < 0.2 && abs(tot_e -3) < 0.5)// shower cut. abs(tot_e - 3.25)< 0.5 && ps_e > 0.2 && abs(rat - 1.2) < 0.2
       {
	////
	sh_flag = kTRUE; //throw up a flag that this is a good electron hit according to the shower detectors
	sh_ps_hit_cnt++;
	////
       }

    //// loop over the PMTs with ADC's that had a signal for this event
    for (Int_t ig=0;ig<GrinchADCNum;ig++) { //since PMTs with no signal are not written to the root file, this "ig" is NOT the ADC channel
      Int_t gindex_adc=adcGID[ig]; // we need to specifically ask which ADC channel this signal is for
      h_grinch_atime_elem->Fill(gindex_adc,adcGAtime[ig]); //the root tree needs this "ig" index and NOT the adc index to get values like the amplitude. 
      h_grinch_amp_elem->Fill(gindex_adc,adcGAmp[ig]);
      h_grinch_adc_chan_vs_num ->Fill(gindex_adc,grinch_adc[ig]);//
      h_grinch_adcMult_elem ->Fill(gindex_adc,adcGMult[ig]);//
      h_grinch_adc[gindex_adc]->Fill(grinch_adc[ig]);// ADC arrays for individual channels 
      h_grinch_amp_p[gindex_adc] ->Fill(adcGAmp[ig]);//
      
      
      //// Map the "gindex_adc" (which is the ADC Chan number) to the "ig index" which is what the branch needs.
      adc_root_index_array[gindex_adc]=ig;
      //// need this in order to go back after finding clusters to look at ADC values and whatnot.
      
      
      if (sh_flag) {//if an event passed the electron cuts
	h_grinch_amp_elem_SHcut ->Fill(adcGID[ig],adcGAmp[ig]);
	h_grinch_atime_elem_SHcut->Fill(adcGID[ig],adcGAtime[ig]);
	h_grinch_adc_chan_vs_num_SHcut ->Fill(adcGID[ig],grinch_adc[ig]);
	h_grinch_adcMult_elem_SHcut ->Fill(adcGID[ig],adcGMult[ig]);
      }
      
      
      if (adcGMult[ig]==1){ //checking ADC multiplicity 
	h_grinch_amp_elem_multcut ->Fill(adcGID[ig],adcGAmp[ig]);
	h_grinch_atime_elem_multcut->Fill(adcGID[ig],adcGAtime[ig]);
	h_grinch_adc_chan_vs_num_multcut ->Fill(adcGID[ig],grinch_adc[ig]);

	if(adcGAtime[ig] > 1230 && adcGAtime[ig]<1300){//might be the wrong atime range. Should check that. 
	  h_grinch_amp_elem_multcut_timecut->Fill(adcGID[ig],adcGAmp[ig]);
	}
      }
    }// end loop over ADC 

    for (int n = 0; n<510;n++) //initialize some more arrays
      {
	hit_flag_array[n]=0;
	mult_array[n]=0;
	root_index_array[n]=0;
      }
    
    
    Double_t goodhit=0;
    //// Loop over the PMTs that had a signal for this event 
    for (Int_t ig=0;ig<GrinchNum;ig++) { //since PMTs with no signal are not written to the root file, this "ig" is NOT the PMT number
      Int_t gindex=tdcGID[ig]; // we need to specifically ask which PMT this signal is for
      if (gindex <511 && gindex >-1) {
	h_grinch_elem->Fill(tdcGID[ig]);//the root tree needs this "ig" index and NOT the PMT index to get values like LE. 
	h_grinch_le_elem->Fill(tdcGID[ig],tdcGLe[ig]);
	h_grinch_te_elem->Fill(tdcGID[ig],tdcGTe[ig]);
	h_grinch_tot_elem->Fill(tdcGID[ig],tdcGTot[ig]);
	h_grinch_mult_elem->Fill(tdcGID[ig],tdcGMult[ig]);
	h_grinch_tot_all->Fill(tdcGTot[ig]);
	h_grinch_le[gindex]->Fill(tdcGLe[ig]);
      
	if (sh_flag) {// if the event passed the shower cuts
	  h_grinch_le_elem_SHcut->Fill(tdcGID[ig],tdcGLe[ig]);
	  h_grinch_te_elem_SHcut->Fill(tdcGID[ig],tdcGTe[ig]);
	  h_grinch_tot_elem_SHcut->Fill(tdcGID[ig],tdcGTot[ig]);
	  h_grinch_mult_elem_SHcut->Fill(tdcGID[ig],tdcGMult[ig]);
	 
	}
      
	if (tdcGMult[ig]==1){//checking multiplicity 
	  h_grinch_le_elem_multcut->Fill(tdcGID[ig],tdcGLe[ig]);
	  h_grinch_te_elem_multcut->Fill(tdcGID[ig],tdcGTe[ig]);
	  h_grinch_tot_elem_multcut->Fill(tdcGID[ig],tdcGTot[ig]);
	}
     
	if (abs(tdcGLe[ig]-900)<20) { //checking if event is within the correct timing for an event 
	  h_grinch_tot_lecut[gindex]->Fill(tdcGTot[ig]);
	}
       
	h_grinch_tot[gindex]->Fill(tdcGTot[ig]);


	// note for eric: I belive this cut is similar to the one in SBSGRINCH. So it would be like we are starting here. 
	if (abs(tdcGLe[ig]-900)<20){ // if the event is within the good timing cut for the grinch 
	  goodhit++;
	  ////////////////////////////
	  hit_flag_array[gindex] = 1; // mark that this PMT has a good hit
	  mult_array[gindex] = tdcGMult[ig]; // save it's multiplicity to the array. 
	  root_index_array[gindex]=ig; // Map the "gindex" (which is the PMT number) to the "ig index" which is what the branch needs.
	                               // I want to be able to go back after finding clusters to look at ADC and TDC and whatnot on those PMTs.   
	  h_grinch_pmt_good_hit ->Fill(gindex);
	  ////////////////////////////
	}

	if (abs(tdcGLe[ig]-900)<20) {
	  h_grinch_tot_all_lecut->Fill(tdcGTot[ig]);
	}	
	     
      }
      h_grinch_ghits->Fill(goodhit);
      
     
    }//end loop over TDC PMTs
  

    clusterstack = Clear_Stack(); //clearing the stack we are using for cluster finding for the new event
    
    
   
    if (temp_cnt <max){// temporary loop here to only do cluster finding on a certain number of events for debugging. 
      for(int m = 0; m<510;m++) //fill with 0 
	{
	  sum_array[m]=0;
	}

      //// Start of Cluster Finding
      for (int pmtcntr = 9; pmtcntr <=500 ; pmtcntr ++) // algorithm doesn't need the top row and bottom row.
	{
	  sum_array[pmtcntr] = Sum_Adjacent_Hits(pmtcntr); // retuns the number PMTs with hits adjacent to the PMT we are looking at
	
	  if (sum_array[pmtcntr] != 0)
	    {
	      //cout<<"sum_array [PMT] "<<pmtcntr<< " = "<<sum_array[pmtcntr]<<endl;
	    }
	  if (sum_array[pmtcntr] >= 2) // if the PMT has 2 or more neighboring PMTs with hits
	    {
	      //cout<< "Possible cluster around PMT "<<  pmtcntr << ": sum = " <<sum_array[pmtcntr]<<endl;	      
	      clusterstack.push(pmtcntr); // put that PMT in the stack to be analyzed to get cluster size
	      h_grinch_pmt_of_interest ->Fill(pmtcntr); 
	    }
	}
      temp_cnt ++; 

     
	  
      //cout<< "EVENT "<<i<<endl;
      //Print_Stack(clusterstack);
      //cout<<"clusterstack"<<endl;

      Int_t ClusterSize = Calculate_Cluster_Size(clusterstack); //send the stack of PMTs that may have a cluster around them
      //// cluster_members_stack is a global variable that is filled in the Calculate_Cluster_Size function. 
      //// I might want to make it a reference parameter instead to make this less confusting. 
      //// cluster_members_stack is now every PMT that is a part of this event with no duplicates. 
      h_grinch_cluster_size ->Fill(ClusterSize);
      //cout<< "cluster size = "<< ClusterSize<<endl;

      g_cluster_flag =kFALSE;
      if(ClusterSize >= 3)
	{
	  //cout<<"EVENT "<<i<<endl;
	  //cout<< "cluster size = "<< ClusterSize<<endl;
	  clustercnt++;
	  g_cluster_flag = kTRUE;
	
	  // Print_Stack(cluster_members_stack);
	  //cout<<"cluster_members_stack in Fill_Event"<<endl;

	  //cout<<"clustercnt = "<<clustercnt<<endl;
	  
      	 	  
	  /*
	    stack <Int_t> copystack = cluster_members_stack;
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
	  
	  Calculate_Cluster_Center(cluster_members_stack,clusterx,clustery); // clusterx and clustery are reference parameters so the function updates them.
	  //// cluster_memebers_stack is a global variable that was filled in Calculate_Cluster_Size. I should prob change it to a referene param so it's not so confusing 
	  //cout<<"cluster center is at "<<clusterx<<", "<<clustery<<endl;
	  h_grinch_cluster_center ->Fill(clusterx,clustery);
	  h_grinch_cluster_center_display ->Fill(clusterx,-clustery);

	  Fill_Cluster_Histos(cluster_members_stack); // filling ADC and TDC histograms for only the PMTs in the cluster. This does not necessarily have to be done in a sperate function like this. 
	  
	  if(!sh_flag) //filling different histos based on if there was a good hit in the shower or not
	    {
	      h_grinch_cluster_center_bad_sh_e ->Fill(clusterx,-clustery);
	      //Fill_Bad_Cluster_Histos(cluster_members_stack);
	      h_grinch_cluster_size_bad ->Fill(ClusterSize);
	    }
	  if(sh_flag)
	    {
	       h_grinch_cluster_center_good_sh_e ->Fill(clusterx,-clustery);
	       // Fill_Cluster_Histos(cluster_members_stack);
	       h_grinch_cluster_size_good ->Fill(ClusterSize);
	    }

	  Int_t vertspread = 0;
	  Int_t horizspread = 0; 

	  Find_Cluster_Extrema(cluster_members_stack,horizspread,vertspread); // function measures the height and width of the cluster
	  h_grinch_cluster_spread->Fill(horizspread,vertspread);
	  
          if(vertspread <=6 && horizspread <=5) //if the cluster is wider or taller than this, it is most likely multiple clusters. 
	    {
	      h_grinch_cluster_center_spreadcut ->Fill(clusterx,-clustery);
	      h_grinch_cluster_size_cuts ->Fill(ClusterSize);
	    }
	  
	  //cout<<"cluster is "<<horizspread<<" wide, "<< vertspread<<" tall"<<endl;


	  //cout<< "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " <<endl;
	}	
	
      //// Asking if the cluster had a good shower hit or if there was a good shower hit and no cluster? 
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

  

  cout<< "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  cout<<"cluster analysis on "<<max<<" events"<<endl;
  cout<<clustercnt<<" clusters found"<<endl;

  cout<<"good shower energy on "<< sh_ps_hit_cnt<<" events" <<endl;
  cout<<"good shower energy and NO good grinch cluster on "<< sh_hit_no_grinch_hit_cnt<<" events" <<endl;
  cout<<"good shower energy and good grinch cluster on "<< grinch_hit_sh_hit_cnt <<" events" <<endl;
  cout<<"NO good shower energy and good grinch cluster on "<< grinch_hit_no_sh_hit_cnt <<" events" <<endl;
  cout<<"NO good shower energy and NO grinch cluster on "<< no_sh_no_grinch_cnt <<" events" <<endl;

  //h_grinch_sh_hit_no_grinch_hit_cnt ->Fill(sh_hit_no_grinch_hit_cnt);
  //h_grinch_grinch_hit_sh_hit_cnt ->Fill(grinch_hit_sh_hit_cnt);
  // h_grinch_grinch_hit_no_sh_hit_cnt ->Fill(grinch_hit_no_sh_hit_cnt);
  //h_grinch_no_sh_no_grinch_cnt ->Fill(no_sh_no_grinch_cnt);
  // h_grinch_sh_ps_hit_cnt ->Fill(sh_ps_hit_cnt);

  Double_t grinch_hit_double = grinch_hit_sh_hit_cnt;
  Double_t sh_hit_double = sh_ps_hit_cnt;
  Double_t percent = 100 * grinch_hit_double/sh_hit_double; 

  h_grinch_efficiency_percent->Fill(percent);

  cout<<"Effiency on Good Hits: "<< grinch_hit_sh_hit_cnt<< "/" << sh_ps_hit_cnt <<" = " << percent <<" %"<<endl; 
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  

  TFile hsimc(outputhist,"recreate");
  HList.Write();

}// end main



////
//// See how many adjacent PMTs also had a hit
////
Int_t Sum_Adjacent_Hits(Int_t PMT)
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
 
  while (!searchstack.empty())
    {
      neighbor = PMT-searchstack.top();
      hit_sum = hit_sum + hit_flag_array[neighbor];     
      searchstack.pop();
    }
  
  //cout<< "hit_sum =" <<hit_sum <<endl;
  return hit_sum;

}//end Sum_Adjacent_Hits


Int_t Calculate_Cluster_Size( stack <Int_t> inputstack)
{ 
  //// This function takes in a stack of PMTs that we already determined may be part of a cluster
  //// It goes through the stack and checks neighbors and makes sure they haven't already been counted
  //// (hmmmmmm now that I'm thinking about it again, there may be an easier way to do this. we already checked neighbors in the process before this, so this method of having to check them again may be a bit redundent. We could make it so the Sum function also fills a stack of the PMT numbers of it's neighbors, then this function would just have to check for duplicates instead of checking neighbors AND checking duplicates). 
  
  stack <Int_t> already_counted;
  Int_t hittube;
  Int_t counter = 0;
  Int_t tube;

  stack <Int_t> searchstack;  
  
  while (!inputstack.empty())
    {
      // cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;      
      // cout<< "hittube " << inputstack.top() <<endl;
      // Print_Stack(inputstack);
      // cout<< "input stack"<<endl;

      hittube = inputstack.top();
      inputstack.pop();      
      if(Is_NOT_In_Stack(already_counted, hittube))//check if the tube has already been counted
      	{
       	  counter ++;
       	  already_counted.push(hittube);
	  // cout<< "hittube has not been counted yet. Adding to already_counted stack and printing."<<endl;
	  // Print_Stack(already_counted);
	  // cout<<"already_counted stack"<<endl;
       	}
      else
	{
	  //cout<< "hittube has already been counted. Not adding it to the already_counted stack. Printing."<<endl;
	}

      if (hittube >= 0 && hittube <= 7) // we are on the top row of the detector
	{
	  //cout<<"top row of the detector"<<endl;
	  if (hittube == 0)//we are on PMT 0, the top left of the detector
	    {
	      //cout<< "upper left corner."<<endl;
	      searchstack = top_row_left_search;
	    }
	  else if (hittube == 7)// we are on PMT 7, the top right of the detector
	    {
	      //cout<< "upper right corner. "<<endl;
	      searchstack = top_row_right_search;
	    }
	  else // we are on the top row and not on the ends
	    {
	      //cout<< "middle of the top row."<<endl;
	      searchstack = top_row_middle_search;
	    }
	}//end top row of the detector
      else if (hittube >= 501 && hittube <= 509) // we are on the bottom row of the detector
	{
	  //cout<< "bottom row of the detector"<<endl;
	  if(hittube == 501)// we are on the PMT 501, the bottom left of the detector
	    {
	      //cout<<"bottom left corner"<<endl;
	      searchstack = bottom_row_left_search;
	    }
	  else if (hittube == 509)// we are on PMT 509, the bottom right of the detector
	    {
	      //cout<<"bottom right corner"<<endl;
	      searchstack = bottom_row_right_search;
	    }
	  else // we are in the middle of the bottom row of the detector
	    {
	      //cout<<"bottom row middle of the detector"<<endl;
	      searchstack = bottom_row_middle_search;
	    }
	}// end bottom row of the detctor 

      else  if(row_array[hittube]%2 == 0)// short row of 8
	{
	  // cout<<"row of 8"<<endl;
	  if (col_array[hittube] == 0)//on left edge short row
	    {
	      // // check neighbors except -1
	      //cout<<"left edge of a row of 8"<<endl;
	      searchstack = short_row_on_left_search;	      
	    }
	  if(col_array[hittube] == 7)//on right edge short row
	    {
	      ////check neighbors except +1
	      //cout<<"right edge of a row of 8"<<endl;
	      searchstack = short_row_on_right_search;	      
	    }
	  else
	    {
	      //cout<<"middle of a row of 8"<<endl;
	      searchstack = middle_search;
	    }
	}//end  short row of 8

      else  // long row of 9
	{
	  //cout<< "row of 9"<<endl;
	  if (col_array[hittube] == 0) // on left edge long row
	    {
	      //cout<<"left edge of a row of 9"<<endl;
	      searchstack = long_row_on_left_search;
	    }
	  if (col_array[hittube] == 8)
	    {
	      //cout<<"right edge of a row of 9"<<endl;
	      searchstack = long_row_on_right_search;
	    }
	  else
	    {
	      //cout<< "middle of a row of 9"<<endl;
	      searchstack = middle_search;
	    }
	} //end row of 9

      
      // cout<< "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
      // cout <<"starting while loop for !searchstack.empty()"<<endl;
      while (!searchstack.empty())
	{
	  //Print_Stack(searchstack);
	  //cout<<"searchstack in while loop"<<endl;
	  tube = hittube + searchstack.top();
	  searchstack.pop();
	  if (hit_flag_array[tube] == 1 && Is_NOT_In_Stack(already_counted,tube))
	    {
	      //cout <<"tube "<<tube<< " has a hit and is not in the already_counted stack"<<endl;
	      counter ++;
	      already_counted.push(tube);
	    } //end if tube checked was not on stack and hit 
	  /*
	  else if(hit_flag_array[tube] == 1 && !Is_NOT_In_Stack(already_counted,tube))
	    {
	      cout<<"tube "<<tube<< "has a hit but is already in the stack"<<endl;
	    }
	  else 
	    {
	      cout<<"tube "<<tube<< " was not hit"<<endl;
	    }
	  */
	  //cout<< "counter = "<<counter <<endl;
	  // Print_Stack(already_counted);
	  // cout<<"already_counted stack after if statment"<<endl;
	}//end while searchstack not empty
      //cout<<"Ended while loop for !searchstack.empty()"<<endl;
      //cout<< "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
     
    }//end of while big stack
  
  // may also want to return a pointer to the "already_counted" stack? 

  cluster_members_stack = Clear_Stack();
  cluster_members_stack = already_counted; // probably should return a pointer or have this as a reference param instead of a global variable like this. 
 
  return counter;
}//end Calculate_Cluster_Size


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
{
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
{
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
{
  Int_t PMT;
  Int_t root_id;
  Int_t ADC_chan;
  Int_t ADC_root_id;
  while (!inputstack.empty())
    {
      PMT = inputstack.top();
      inputstack.pop();
      root_id = root_index_array[PMT];//get the index that the branches need in order to look at the specific PMT for this event
      h_grinch_cluster_le_elem -> Fill(PMT,tdcGLe[root_id]);
      h_grinch_cluster_te_elem -> Fill(PMT,tdcGTe[root_id]);
      h_grinch_cluster_tot_elem -> Fill(PMT,tdcGTot[root_id]);
      h_grinch_cluster_mult_elem -> Fill(PMT,tdcGMult[root_id]);
      h_grinch_cluster_tot[PMT]->Fill(tdcGMult[root_id]);
      

      //see if the PMT has a corresponding ADC channel 
      ADC_chan = tdc_to_adc_map_array[PMT];
      if(ADC_chan !=-1)//the PMT has a corresponding ADC channel. 
	{
	  ADC_root_id = adc_root_index_array[ADC_chan];
	  h_grinch_cluster_atime_elem ->Fill(ADC_chan, adcGAtime[ADC_root_id]);
	  h_grinch_cluster_amp_elem ->Fill(ADC_chan, adcGAmp[ADC_root_id]);
	  h_grinch_cluster_amp_p[ADC_chan] ->Fill(adcGAmp[ADC_root_id]);
	  h_grinch_cluster_int_elem ->Fill(ADC_chan, grinch_adc[ADC_root_id]);
	  h_grinch_cluster_adcMult_elem ->Fill(ADC_chan,adcGMult[ADC_root_id]); 
	  h_grinch_cluster_adc_tdc ->Fill(tdcGTot[root_id],adcGAmp[ADC_root_id]); 
	  h_grinch_cluster_adcInt_tdc ->Fill(tdcGTot[root_id],grinch_adc[ADC_root_id]); 
	  h_grinch_cluster_int_amp ->Fill(adcGAmp[ADC_root_id],grinch_adc[ADC_root_id]);//just testing the amp and integral are proportional to eachother 
	  if(ADC_chan == 23)
	    {
	      h_grinch_cluster_adc_tdc_test ->Fill(tdcGTot[root_id],adcGAmp[ADC_root_id]);
	    }
	}
    }
}//end Fill_Cluster_Histos


void Fill_Bad_Cluster_Histos(stack <Int_t> inputstack)
{
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


void Print_Stack(stack <Int_t> inputstack)
{
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
{
  stack <Int_t> clearstack; //default constructor is an empty stack
  return clearstack;
}// end Clear_Stack
