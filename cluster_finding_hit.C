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
//// Feb 2023
//// msatnik@email.wm.edu
////
//// Multiple Clusters!
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


stack <Int_t> cluster_stack[32];
stack <Int_t> cluster_stack_ordered[32];
stack <Int_t> neighbors_stack[510];
stack <Int_t> check_stack;
Int_t check_status[510];
Int_t cluster_stack_cnt = 0;
Bool_t clusterflag=kFALSE;

Int_t ClusterSize[32];
Double_t cluster_vert[32];
Double_t cluster_horiz[32];
Double_t cluster_xpos[32];
Double_t cluster_ypos[32];
Int_t horizspread[32];
Int_t vertspread[32];
Double_t cluster_mean_time[32];
Double_t cluster_variance[32];

Double_t LE_offset[510] = {0};

Int_t ClusterSize_ordered[32];
Double_t cluster_vert_ordered[32];
Double_t cluster_horiz_ordered[32];
Double_t cluster_xpos_ordered[32];
Double_t cluster_ypos_ordered[32];
Int_t horizspread_ordered[32];
Int_t vertspread_ordered[32];
Double_t cluster_mean_time_ordered[32];
Double_t cluster_variance_ordered[32];
Double_t grinch_dx[32];
Double_t grinch_dy[32];
Double_t projx;
Double_t projy;

struct Cluster{
  Int_t size = 0;
  Int_t PMT[32] = {0};
  Int_t rootID[32] = {0};
  Double_t LE[32] = {0};
  Double_t ToT[32] = {0};
  Double_t sqrt_variance = 0;
  Double_t mean_time = 0;
  Double_t row = 0;
  Double_t col = 0;
  Double_t xpos = 0;
  Double_t ypos = 0;
  Double_t horizspread = 0;
  Double_t vertspread = 0;
  Double_t ToT_sum = 0;
};




Int_t row;
Int_t col;


const Int_t N_ROW=60;

const Double_t grinch_distance = 0.48; //distance from focal plane. (Stolen from Andrew's code. May be different for different kinematics). 

Int_t row_array[510];
Int_t col_array[510];
Int_t tdc_to_adc_map_array[510]; // Needs to be updated when ADCs are moved to different Nino cards. 
Double_t xpos_array[510];
Double_t ypos_array[510];

Double_t tr_x[10];
Double_t tr_y[10];
Double_t tr_n;
Double_t tr_vz[100];
Double_t tr_tg_th[100];
Double_t tr_tg_ph[100];
Double_t tr_th[100];
Double_t tr_ph[100];
Double_t gem_track_nhits;
Double_t tr_p[100];
Double_t hcal_e;

UInt_t fTrigBits;
UInt_t fRun;

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

Double_t kine_W; 

Double_t hodo_tmean;

Int_t testcnt = 0;



Int_t NumHits;
Double_t tdcGHitLe[1000];
Double_t tdcGHitPmtnum[1000];
Double_t tdcGHitTot[1000];
Double_t tdcGHitCol[1000];
Double_t tdcGHitRow[1000];
Double_t tdcGHitX[1000];
Double_t tdcGHitY[1000];


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
TH1F* h_sh_e = new TH1F("h_sh_e"," ; Total Shower E ",300,0.0,3.0);
TH1F* h_sh_tot_e = new TH1F("h_sh_tot_e"," ; Total PRe+Shower E ",300,0.0,5.0);
TH1F* h_ratio_e = new TH1F("h_ratio_e"," ; Total PRe+Shower E/ track norm ",300,0.0,2.5);
TH1F* h_gTrigBits = new TH1F("h_gTrigBits", "; gTrigBits;",35,0,35);
TH1F* h_fTrigBits = new TH1F("h_fTrigBits", "; fTrigBits;",35,0,35);
TH2F* h_ps_sh_e = new TH2F("h_ps_sh_e", "; sh_e ; ps_e ;",300,0,3,200,0,2);
TH2F* h_ps_e_tr_p = new TH2F("h_ps_e_tr_p",";tr_p;ps_e;", 200,0,5,200,0,2);



TH1F *h_hodo_tmean = new TH1F("h_hodo_tmean",";bb.hodotdc.clus.tmean;",40,-20,20); 
TH1F* h_tr_x = new TH1F("h_tr_x", "bb.tr.x[0]", 200,-1,1);
TH1F* h_tr_y = new TH1F("h_tr_y", "bb.tr.y[0]", 200, -1,1);
TH1F* h_tr_p = new TH1F("h_tr_p","bb.tr.p[0]",200,0,5);
TH1F* h_tr_vz = new TH1F("h_tr_vz","bb.tr.vz[0])",200,-0.2,0.2);
TH1F* h_hcal_e = new TH1F("h_hcal_e","sbs.hcal.e" ,200,0,1);
TH1F* h_gem_track_nhits = new TH1F("h_gem_track_nhits","bb.gem.track.nhits",10,0,10);
TH1F* h_tr_tg_ph = new TH1F("h_tr_tg_ph", "bb.tr.tg_ph[0]",200,-0.2,0.2);
TH1F* h_tr_tg_th = new TH1F("h_tr_tg_th", "bb.tr.tg_th[0]",200,-0.3,0.3);
TH1F* h_tr_ph = new TH1F("h_tr_ph", "bb.tr.ph[0]",200,-0.2,0.2);
TH1F* h_tr_th = new TH1F("h_tr_th", "bb.tr.th[0]",200,-0.3,0.3);
TH1F* h_tr_n = new TH1F("h_tr_n","bb.tr.n",10,0,10);

TH1F* h_grinch_tr_vert_proj =  new TH1F("h_grinch_tr_vert_proj","track x projected to grinch window", 200,-2,2);
TH1F* h_grinch_tr_horiz_proj = new TH1F("h_grinch_tr_horiz_proj","track y projected to grinch window",200,-2,2);

TH1F* h_grinch_cluster_dx = new TH1F("h_grinch_cluster_dx","grinch cluster center x -  track x projection", 50,-1,1);
TH1F* h_grinch_cluster_dy = new TH1F("h_grinch_cluster_dy","grinch cluster center y -  track y projection", 50,-1,1);
TH1F* h_grinch_cluster_dx_best = new TH1F("h_grinch_cluster_dx_best","grinch cluster center x -  track x projection for best cluster", 50,-1,1);
TH1F* h_grinch_cluster_dy_best = new TH1F("h_grinch_cluster_dy_best","grinch cluster center y -  track y projection for best cluster", 50,-1,1);
TH1F* h_grinch_best_cluster_index = new TH1F("h_grinch_best_cluster_index","index of best grinch cluster (ranked by size)",10,0,10);
TH1F* h_grinch_cluster_dx_rejected = new TH1F("h_grinch_cluster_dx_rejected","grinch cluster center x -  track x projection for rejected clusters", 50,-1,1);
TH1F* h_grinch_cluster_dy_rejected = new TH1F("h_grinch_cluster_dy_rejected","grinch cluster center y -  track y projection for rejected clusters", 50,-1,1);
TH1F* h_grinch_rejected_cluster_index = new TH1F("h_grinch_rejected_cluster_index","index of rejected grinch clusters (ranked by size)",10,0,10);

TH1F* h_grinch_cluster_dx_multi_best = new TH1F("h_grinch_cluster_dx_multi_best","grinch cluster center x -  track x projection for best cluster", 50,-1,1);
TH1F* h_grinch_cluster_dy_multi_best = new TH1F("h_grinch_cluster_dy_multi_best","grinch cluster center y -  track y projection for best cluster", 50,-1,1);
TH1F* h_grinch_cluster_dx_multi_rejected = new TH1F("h_grinch_cluster_dx_multi_rejected","grinch cluster center x -  track x projection for rejected clusters", 50,-1,1);
TH1F* h_grinch_cluster_dy_multi_rejected = new TH1F("h_grinch_cluster_dy_multi_rejected","grinch cluster center y -  track y projection for rejected clusters", 50,-1,1);

TH1F* h_grinch_cluster_xdiff = new TH1F(" h_grinch_cluster_xdiff ", "x seperation between clusters",50,0,1 );
TH1F* h_grinch_cluster_ydiff = new TH1F(" h_grinch_cluster_ydiff ", "y seperation between clusters", 50,0,1);
TH2F* h_grinch_cluster_xdiff_ydiff = new TH2F("h_grinch_cluster_xdiff_ydiff", ";y seperation; x seperation", 50, 0, 1, 50, 0, 1);
TH1F* h_grinch_cluster_rdiff = new TH1F(" h_grinch_cluster_rdiff ", "radial seperation between clusters",20,0,0.2 );
TH1F* h_grinch_cluster_radius_ana = new TH1F(" h_grinch_cluster_radius_ana ", "r1+r2 - seperation",50,-0.1,0.1 );

TH1F* h_grinch_clustercnt = new TH1F("h_grinch_clustercnt","number of clusters in event with more than 2 pmts",10,0,10);

TH1F* h_grinch_cluster_tr_vert_proj =  new TH1F("h_grinch_cluster_tr_vert_proj","track x projected to grinch window", 200,-2,2);
TH1F* h_grinch_cluster_tr_horiz_proj = new TH1F("h_grinch_cluster_tr_horiz_proj","track y projected to grinch window",200,-2,2);

TH1F* h_grinch_cluster_center_xpos = new TH1F("h_grinch_cluster_center_xpos","cluster x at grinch window",100,-1,1);
TH1F* h_grinch_cluster_center_ypos = new TH1F("h_grinch_cluster_center_ypos","cluster y at grinch window",100,-1,1);

TH1F* h_grinch_projx =  new TH1F("h_grinch_projx","track x projected to grinch window",100,-1,1);
TH1F* h_grinch_projy =  new TH1F("h_grinch_projy","track y projected to grinch window",100,-1,1);


TH2F* h_grinch_cluster_tr_vert_proj_vs_cluster_vert =  new TH2F("h_grinch_cluster_tr_vert_proj_vs_cluster_vert",";track x projected to grinch window ; cluster vertical", 100,-1,1,59,0,60);
TH2F* h_grinch_cluster_tr_horiz_proj_vs_cluster_horiz = new TH2F("h_grinch_cluster_tr_horiz_proj_vs_cluster_horiz",";track y projected to grinch window ; cluster horizontal",100,-1,1,17,0,9);

TH2F* h_grinch_cluster_ToT_sum_vs_size = new TH2F("h_grinch_cluster_ToT_sum_vs_size", ";ToT sum ; cluster size", 200,0,200, 20,0,20);

TH2F* h_grinch_mirror_2_projy_ypos = new TH2F ("h_grinch_mirror_2_projy_","mirror 2 ;track y projection; cluster y position;", 100,-0.25,0.25,100,-0.25,0.25);
TH2F* h_grinch_mirror_1_projy_ypos = new TH2F ("h_grinch_mirror_1_projy_","mirror 1 ;track y projection; cluster y position;", 100,-0.25,0.25,100,-0.25,0.25);
TH2F* h_grinch_mirror_3_projy_ypos = new TH2F ("h_grinch_mirror_3_projy_","mirror 3 ;track y projection; cluster y position;", 100,-0.25,0.25,100,-0.25,0.25);
TH2F* h_grinch_mirror_4_projy_ypos = new TH2F ("h_grinch_mirror_4_projy_","mirror 4 ;track y projection; cluster y position;", 100,-0.25,0.25,100,-0.25,0.25);

TH2F* h_grinch_cluster_projx_xpos_best = new TH2F("h_grinch_cluster_projx_xpos_best","; cluster x position best cluster ; projected x at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projy_ypos_best = new TH2F("h_grinch_cluster_projy_ypos_best","; projected y at grinch from track; cluster y position best cluster ;",100,-0.2,0.2,100,-0.2,0.2);
TH2F* h_grinch_cluster_projx_xpos = new TH2F("h_grinch_cluster_projx_xpos","; projected x at grinch window from track ;cluster x position",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projy_ypos = new TH2F("h_grinch_cluster_projy_ypos","; projected y at grinch window from track ;cluster y position",100,-0.2,0.2,100,-0.2,0.2);
TH2F* h_grinch_cluster_projx_xpos_rejected = new TH2F("h_grinch_cluster_projx_xpos_rejected","; cluster x position rejected clusters ; projected x at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projy_ypos_rejected = new TH2F("h_grinch_cluster_projy_ypos_rejected","; cluster y position rejected clusters ; projected y at grinch window from track",100,-1,1,100,-1,1);

TH2F* h_grinch_cluster_projx_xpos_multi_best = new TH2F("h_grinch_cluster_projx_xpos_multi_best","; cluster x position best cluster ; projected x at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projy_ypos_multi_best = new TH2F("h_grinch_cluster_projy_ypos_multi_best","; cluster y position best cluster ; projected y at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projx_xpos_multi_rejected = new TH2F("h_grinch_cluster_projx_xpos_multi_rejected","; cluster x position rejected clusters ; projected x at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projy_ypos_multi_rejected = new TH2F("h_grinch_cluster_projy_ypos_multi_rejected","; cluster y position rejected clusters ; projected y at grinch window from track",100,-1,1,100,-1,1);


TH2F* h_grinch_cluster_center_vert_tr = new TH2F("h_grinch_cluster_center_vert_tr",";GRINCH cluster center vertical ; track x;", 120,0,60,250,-1,1);
TH2F* h_grinch_cluster_center_horiz_tr = new TH2F("h_grinch_cluster_center_horiz_tr",";GRINCH cluster center horizontal; track y;", 18,0,9,75,-0.3,0.3);
TH2F* h_grinch_cluster_center_horiz_tr_display = new TH2F("h_grinch_cluster_center_horiz_tr_display",";cluster center horizontal; bb.tr.y[0];", 18,-9,0,75,-0.3,0.3);

TH2F* h_grinch_CUT_cluster_center_vert_tr = new TH2F("h_grinch_CUT_cluster_center_vert_tr",";cluster center vertical ; bb.tr.x[0];", 120,0,60,250,-1,1);
TH2F* h_grinch_CUT_cluster_center_horiz_tr = new TH2F("h_grinch_CUT_cluster_center_horiz_tr",";cluster center horizontal; bb.tr.y[0];", 18,0,9,75,-0.3,0.3);

TH2F* h_grinch_cluster_track_xy = new TH2F("h_grinch_cluster_track_xy",";bb.tr.y[0];bb.tr.x[0];", 60,-0.3,0.3,200,-1,1);

TH2F* h_good_event_track_xy = new TH2F("h_good_event_track_xy",";bb.tr.y[0];bb.tr.x[0];",  60,-0.3,0.3,200,-1,1);

TH2F* h_good_electron_no_grinch_hit_track_xy = new TH2F("h_good_electron_no_grinch_hit_track_xy",";bb.tr.y[0];bb.tr.x[0];",  60,-0.3,0.3,200,-1,1);

TH2F* h_grinch_eff_xy = new TH2F("h_grinch_eff_xy",";bb.tr.y[0];bb.tr.x[0];",  32,-0.2,0.2,112,-0.7,0.7);

TH2F* h_grinch_cluster_track_subtraction = new TH2F("h_grinch_cluster_track_subtraction",";clusterx - tr_y[0]; clustery -tr_x[0];", 500,-1, 9,500,-1,60);
TH1F* h_grinch_cluster_center_vert_tr_subtraction =  new TH1F("h_grinch_cluster_center_vert_tr_subtraction", "clustery - tr_x[0]", 500,-1,60);
TH1F* h_grinch_cluster_center_horiz_tr_subtraction = new TH1F("h_grinch_cluster_center_horiz_tr_subtraction ", "clusterx - tr_y[0]", 500, -1,9);


TH2F* h_grinch_CUT_cluster_center_display = new TH2F("h_grinch_CUT_cluster_center_display", "; horizontal ; vertical ;",16,0,8,60,-60,0);


TH1F* h_grinch_cluster_size = new TH1F("h_grinch_cluster_size" ,"; Cluster Size ",18,2,20);
TH2F* h_grinch_cluster_center = new TH2F("h_grinch_cluster_center", "; horizontal ; vertical ;",16,0,8,60,0,60);
TH2F* h_grinch_cluster_center_display = new TH2F("h_grinch_cluster_center_display", "; horizontal ; vertical ;",16,0,8,60,-600);
TH1F* h_grinch_cluster_size_best = new TH1F("h_grinch_cluster_size_best" ,"; Cluster Size ",18,2,20);

TH1F* h_grinch_cluster_size_cut = new TH1F("h_grinch_cluster_size_cut" ,"; Cluster Size, mirror 2",18,2,20);

TH1F* h_grinch_cluster_tmean= new TH1F("h_grinch_cluster_tmean",";Cluster Mean Time (LE)", 240,-20,220);
TH1F* h_grinch_cluster_ToT_sum = new TH1F(" h_grinch_cluster_ToT_sum","Cluster ToT sum", 200,0,200);
TH1F* h_grinch_cluster_variance = new TH1F("h_grinch_cluster_variance","; Cluster Variance (LE)",40,0,20);
TH2F* h_grinch_cluster_vert_tmean =new TH2F("h_grinch_cluster_vert_tmean", ";cluster vert; cluster mean time (LE)",60,0,60,40,180,220);
TH2F* h_grinch_cluster_horiz_tmean = new TH2F("h_grinch_cluster_horiz_tmean", ";cluster horiz; cluster mean time (LE)",16,0,8,40,180,220);
TH2F* h_grinch_cluster_vert_variance  =new TH2F("h_grinch_cluster_vert_variance", ";cluster vert; cluster variance (LE)",60,0,60,40,0,20);
TH2F* h_grinch_cluster_horiz_variance = new TH2F("h_grinch_cluster_horiz_variance", ";cluster horiz; cluster variance (LE)",16,0,8,40,0,20);

TH2F* h_grinch_cluster_tmean_size = new TH2F("h_grinch_cluster_size_tmean",";size; tmean", 20,0,20,240,-20,220);

TH1F* h_grinch_cluster_tmean_ordered= new TH1F("h_grinch_cluster_tmean_ordered",";Cluster Mean Time (LE)", 40,180,220);
TH1F* h_grinch_cluster_variance_ordered = new TH1F("h_grinch_cluster_variance_ordered","; Cluster Variance (LE)",40,0,20);
TH1F* h_grinch_cluster_variance_mulit_best = new TH1F("h_grinch_cluster_variance_multi_best","; Cluster Variance (LE) for best cluster",40,0,20);
TH1F* h_grinch_cluster_variance_mulit_rejected = new TH1F("h_grinch_cluster_variance_multi_rejected","; Cluster Variance (LE) for rejected clusters",40,0,20);
TH1F* h_grinch_cluster_size_ordered = new TH1F("h_grinch_cluster_size_ordered","cluster size of [0] cluster (excluding 2)",15,0,15);
TH2F* h_grinch_cluster_vert_tmean_ordered =new TH2F("h_grinch_cluster_vert_tmean_ordered", ";cluster vert; cluster mean time (LE)",60,0,60,40,180,220);
TH2F* h_grinch_cluster_horiz_tmean_ordered = new TH2F("h_grinch_cluster_horiz_tmean_ordered", ";cluster horiz; cluster mean time (LE)",16,0,8,40,180,220);
TH2F* h_grinch_cluster_vert_variance_ordered  =new TH2F("h_grinch_cluster_vert_variance_ordered", ";cluster vert; cluster variance (LE)",60,0,60,40,0,20);
TH2F* h_grinch_cluster_horiz_variance_ordered = new TH2F("h_grinch_cluster_horiz_variance_ordered", ";cluster horiz; cluster variance (LE)",16,0,8,40,0,20);


TH1F* h_grinch_cluster_cnt= new TH1F("h_grinch_cluster_cnt","number GRINCH clusters in cluster event (including clusters of 2)",16,0,16);


TH2F* h_grinch_cluster_center_spreadcut = new TH2F("h_grinch_cluster_center_spreadcut",";horizontal; vertical;",16,0,8,60,-60,0);
TH2F* test_coord_histo = new TH2F(" test_coord_histo", "; horizontal ; vertical ;",18,0,9,60,0,60);
TH1F* h_grinch_pmt_good_hit = new TH1F("h_grinch_pmt_good_hit",";PMT ;",510,0,510);
TH2F* h_grinch_cluster_spread = new TH2F("h_grinch_cluster_spread",";cluster width;  cluster height;",8,1,9,59,1,60);



TH1F* h_grinch_cluster_elem = new TH1F("h_grinch_cluster_elem",";GRINCH TDC elemID;", 510,0,510);
TH2F* h_grinch_cluster_le_elem = new TH2F("h_grinch_cluster_le_elem","; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,420,-20,400);

TH2F* h_grinch_cluster_le_elem_corr = new TH2F("h_grinch_cluster_le_elem_corr","; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,200,0,400);

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


TGraph* background_graph = new TGraph();
TGraph* rate_graph = new TGraph();
TGraphErrors* offset_graph = new TGraphErrors();
TGraphErrors* mean_graph = new TGraphErrors();
TGraphErrors* sigma_graph = new TGraphErrors();
TGraphErrors* sigma_squared_graph = new TGraphErrors();
TGraphErrors* offset_cluster_graph = new TGraphErrors();
TGraphErrors* mean_cluster_graph = new TGraphErrors();
TGraphErrors* sigma_cluster_graph = new TGraphErrors();
TGraphErrors* offset_led_graph = new TGraphErrors();
TGraphErrors* mean_led_graph = new TGraphErrors();
TGraphErrors* sigma_led_graph = new TGraphErrors();
TGraphErrors* mean_led_cluster_subtract_graph = new TGraphErrors();


TH2F* h_grinch_le_elem = new TH2F("h_grinch_le_elem"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,1220,-20,1200);
TH2F* h_grinch_tot_elem = new TH2F("h_grinch_tot_elem"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) ",510,0,510,101,0,100);

TH2F* h_grinch_le_elem_hodo = new TH2F("h_grinch_le_elem_hodo"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) - hodo tmean ",510,0,510,2600,0,2600);


TH2F* h_grinch_le_elem_LED = new TH2F("h_grinch_le_elem_LED"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,2000,0,2000);

TH2F* h_grinch_mult_elem = new TH2F("h_grinch_mult_elem"," ; GRINCH TDC elemID ; GRINCH TDC Multiplicty ",510,0,510,10,0,10);
TH1F* h_grinch_elem = new TH1F("h_grinch_elem",";GRINCH elem ID;",510,0,510);
TH2F* h_grinch_amp_elem = new TH2F("h_grinch_amp_elem"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) ",64,0,63,1000,-100,1900);  
TH1F* h_grinch_atime = new TH1F("h_grinch_atime","GRINCH ADC aTime",1000,0,999);
TH2F* h_grinch_atime_elem = new TH2F("h_grinch_atime_elem", "GRINCH ADC elemID; GRINCH ADC aTime",64,0,63,1000,0,999);

TH2F* h_grinch_amp_elem_LED = new TH2F("h_grinch_amp_elem_LED"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) ",64,0,63,1000,-100,1900); 
TH1F* h_grinch_atime_LED = new TH1F("h_grinch_atime_LED","GRINCH ADC aTime LED",1000,0,999);
TH2F* h_grinch_atime_elem_LED = new TH2F("h_grinch_atime_elem_LED", "GRINCH ADC elemID; GRINCH ADC aTime LED",64,0,63,1000,0,999);


TH1F* h_grinch_xhit =  new TH1F("h_grinch_xhit","bb.grinch_tdc.hit.xhit",59,-0.9145,0.9145);
TH1F* h_grinch_yhit = new TH1F("h_grinch_yhit","bb.grinch_tdc.hit.yhit",17,-0.124,0.124);

 

TH2F* h_grinch_background_elem = new TH2F("h_grinch_background_elem","; GRINCH TDC elemID; Background",510,0,510,600,0,600);

TH2F* h_grinch_rate_elem = new TH2F("h_grinch_rate_elem"," ; GRINCH TDC elemID; calc. bkgnd rate",510,0,510,10000,0,1000000);

TH1F* h_grinch_hitslistlength = new TH1F("h_grinch_hitslistlength",";total fires in the multi-hit tdc;", 500,0,500);

//TH2F* h_grinch_hits_elem = new TH2F("h_grinch_hits_elem", " ; ; ")

TH1F* h_grinch_le_all = new TH1F("h_grinch_le_all","; GRINCH LE ALL ;", 2500,0,2500);

TH1F* h_grinch_le_all_LED = new TH1F("h_grinch_le_all_LED","; GRINCH LE ALL ;", 2500,0,2500);


TH1F* h_grinch_tot_all =new TH1F("h_grinch_tot_all","; GRINCH TOT ALL ;", 101,0,100);

TH2F* h_grinch_hit_le_elem = new TH2F("h_grinch_hit_le_elem"," ; GRINCH TDC elemID ; GRINCH TDC HIT LE (ns) ",510,0,510,2700,-100,2600);
TH1F* h_grinch_hit_le_all = new TH1F("h_grinch_hit_le_all","; GRINCH LE HIT ALL ;", 2500,0,2500);

TH1F* h_grinch_cluster_le_all = new TH1F("h_grinch_cluster_le_all","; GRINCH CLUSTER LE ALL ;", 1520,-20,1500);
TH1F* h_grinch_cluster_le_all_corr = new TH1F("h_grinch_cluster_le_all_corr","; GRINCH CLUSTER LE ALL ;", 2500,0,2500);
TH1F* h_grinch_cluster_tot_all = new TH1F("h_grinch_cluster_tot_all","; GRINCH CLUSTER TOT ALL ;", 101,0,100);

TH1F* h_grinch_cluster_tot[511];
TH1F* h_grinch_cluster_le[511];
TH1F* h_grinch_cluster_le_corr[511];
TH1F* h_grinch_cluster_amp[64];//maria 
TH1F* h_grinch_amp[64];
TH2F* h_grinch_cluster_le_tot[511];

Int_t sh_ps_hit_cnt = 0; //counter for the number of events that pass the shower cut
Int_t sh_hit_no_grinch_hit_cnt = 0; //counter for the number of events that pass the shower cut but do not have a grinch cluster
Int_t grinch_hit_sh_hit_cnt = 0; // counter for the number of event that pass the shower cut and have a grinch cluster
Int_t grinch_hit_no_sh_hit_cnt = 0;// counter for the number of events that have a grinch cluster but do not pass the shower cut
Int_t no_sh_no_grinch_cnt = 0;// counter for the number of events that do not pass the shower cut and do not have a grinch cluster
Int_t cluster_in_tr_cut_cnt = 0;
Int_t tr_cut_cnt = 0; 

Int_t cluster_in_tr_cut_cnt_2D[500][100] = {0};
Int_t tr_cut_cnt_2D[500][100] = {0};
Double_t nsteps_trx = 112;// 56 28;
Double_t nsteps_try = 32;//16 8;
Double_t stepsize = 0.0125; //0.025 0.5
Double_t min_trx = -0.7;
Double_t min_try = -0.2;
Double_t trx_steps_array[500] = {0};
Double_t try_steps_array[500] = {0};


TH1F* h_kine_W =  new TH1F("h_kine_W", "; W ;", 500,-10,10);


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
void Calculate_Cluster_Center_Row_Col(stack <Int_t> inputstack, double& horiz, double& vert);
void Calculate_Cluster_Center_ypos_xpos(stack <Int_t> inputstack, double& ypos, double& xpos);
void Find_Cluster_Extrema(stack <Int_t> inputstack, Int_t& horizspread, Int_t& vertspread);
void Fill_Cluster_Histos(stack <Int_t> inputstack, Int_t index);
void Fill_Cluster_Histos_struct(struct Cluster cluster_struct_array[], Int_t entries);
void make_tdc_to_adc_map_array(); // Probably uncessary: should add to database instead
Int_t Stack_Size(stack <Int_t> inputstack);
Double_t offset_gaus(Double_t *x, Double_t *par);
void Palette1();
void set_color_env();
void Calculate_Cluster_Time(stack <Int_t> inputstack, Double_t& mean, Double_t& variance);
void SelectionSort(Int_t a[], Int_t tracker[], Int_t n);
void StackToArray(stack <Int_t> s, Int_t a[]);



void cluster_finding_hit(TString basename="",Int_t nrun=2043,TString configfilename="run_list_grinch.txt", Int_t entries = -1){ //MAIN
  if (basename=="") {
    cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
    cin >> basename;
  }
  TString fullname = Form("%s_%d",basename.Data(),nrun);
  // gStyle->SetPalette(kDarkBodyRadiator);
  //gStyle->SetPalette(1);
   // Palette1();
  set_color_env();
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
      cout<< " may take a moment to load..." <<endl;
      fchain->Add(currentline);
    }   
  } 


  // reading in LE offsets from a config file

  Bool_t CORRECTING = kFALSE; 

  if(CORRECTING)
    {
      ifstream configfile2("/w/halla-scshelf2102/sbs/msatnik/GRINCH_macros/textfiles/offsets_3038_corrections.txt");
      TString currentline2;
      Int_t tokencounter = 0;
      Int_t pmtcounter = 0;
      TString sval;
      while( currentline2.ReadLine( configfile2 ) && !currentline.BeginsWith("#") )
	{
	  //cout<< "current line: "<< currentline2<< endl;
	  TObjArray *tokens = currentline2.Tokenize(" ");
	  Int_t ntokens = tokens ->GetEntries();
	  //cout << "ntokens :" << ntokens <<endl;
	  if (ntokens < 16)
	    {
	      cout<<"something wrong with the le offset file"<<endl;
	    }
	  else for (Int_t i = 0 ; i <16 ; i++)
		 {
		   sval = ( (TObjString*)(*tokens)[i] )->GetString();
		   LE_offset[pmtcounter] = sval.Atof();
		   pmtcounter ++;
		 }
	  delete tokens;
	}

    }

 
  // for (Int_t i = 0; i< 510 ; i++)
  //   {
  //     if (i%16 == 0)
  // 	{
  // 	  cout<< endl;
  // 	}
  //      cout << LE_offset[i] << " ";
  //   }
  // cout<<endl;
  // cout <<  "LE_offset[1] "<< LE_offset[1] <<endl;
  

  //  Int_t hcalNum;
  //fchain->SetBranchAddress("Ndata.sbs.hcal.tdcelemID",&hcalNum) ;
  
  fchain ->SetBranchStatus("*",0);


  fchain ->SetBranchStatus("bb.ps.e",1);
  fchain->SetBranchAddress("bb.ps.e",&ps_e) ; // preshower energy

  fchain->SetBranchStatus("bb.sh.e",1);
  fchain->SetBranchAddress("bb.sh.e",&sh_e) ;  // shower energy

  // fchain->SetBranchStatus("BB.gold.p",1);
  // fchain->SetBranchAddress("BB.gold.p",&pmom) ; //

  // fchain->SetBranchStatus("BB.gold.th",1) ; //
  // fchain->SetBranchAddress("BB.gold.th",&xptar) ; //

  // fchain->SetBranchStatus("BB.gold.ph",1) ; //
  // fchain->SetBranchAddress("BB.gold.ph",&yptar) ; //

  // fchain->SetBranchStatus("BB.gold.y",1) ; //
  // fchain->SetBranchAddress("BB.gold.y",&ytar) ; //

  fchain ->SetBranchStatus("bb.tr.x",1);
  fchain ->SetBranchAddress("bb.tr.x",&tr_x);

  fchain ->SetBranchStatus("bb.tr.y",1);
  fchain ->SetBranchAddress("bb.tr.y",&tr_y);

  fchain ->SetBranchStatus("bb.tr.n", 1);
  fchain ->SetBranchAddress("bb.tr.n", &tr_n);

  fchain ->SetBranchStatus("bb.tr.vz",1);
  fchain ->SetBranchAddress("bb.tr.vz", &tr_vz);

  fchain ->SetBranchStatus("bb.tr.tg_th",1);
  fchain ->SetBranchAddress("bb.tr.tg_th",&tr_tg_th);

  fchain ->SetBranchStatus("bb.tr.tg_ph",1);
  fchain ->SetBranchAddress("bb.tr.tg_ph",&tr_tg_ph);

  fchain ->SetBranchStatus("bb.tr.th",1);
  fchain ->SetBranchAddress("bb.tr.th",&tr_th);

  fchain ->SetBranchStatus("bb.tr.ph",1);
  fchain ->SetBranchAddress("bb.tr.ph",&tr_ph);

  fchain ->SetBranchStatus("bb.gem.track.nhits",1);
  fchain ->SetBranchAddress("bb.gem.track.nhits",&gem_track_nhits);

  fchain ->SetBranchStatus("bb.tr.p",1);
  fchain ->SetBranchAddress("bb.tr.p",&tr_p);

  fchain ->SetBranchStatus("sbs.hcal.e", 1);
  fchain ->SetBranchAddress("sbs.hcal.e", &hcal_e);

  fchain->SetBranchStatus("e.kine.W2", 1);
  fchain->SetBranchAddress("e.kine.W2", &kine_W);

  fchain ->SetBranchStatus("bb.hodotdc.clus.tmean",1); 
  fchain ->SetBranchAddress("bb.hodotdc.clus.tmean",&hodo_tmean); 


 



  fchain->SetBranchStatus("Ndata.bb.grinch_tdc.hit.pmtnum",1) ;
  fchain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.pmtnum",&GrinchNum) ;// 
  //was Ndata.bb.grinch_tdc.tdcelemID
  fchain->SetBranchStatus("bb.grinch_tdc.hit.pmtnum",1) ;     
  fchain->SetBranchAddress("bb.grinch_tdc.hit.pmtnum",&tdcGID) ; // The PMT number 
  fchain->SetBranchStatus("bb.grinch_tdc.hit.time",1) ; 
  fchain->SetBranchAddress("bb.grinch_tdc.hit.time",&tdcGLe) ; // leading edge
  //fchain -> SetBranchStatus("bb.grinch_tdc.tdc_te",1) ; // trailing edge
  //fchain -> SetBranchAddress("bb.grinch_tdc.tdc_te",&tdcGTe) ; // trailing edge
  fchain->SetBranchStatus("bb.grinch_tdc.hit.amp",1) ; 
  fchain->SetBranchAddress("bb.grinch_tdc.hit.amp",&tdcGTot) ; // tdc time over threshold (te-le)
  fchain -> SetBranchStatus("bb.grinch_tdc.hit.xhit",1);
  fchain -> SetBranchAddress("bb.grinch_tdc.hit.xhit",&tdcGHitX);
  fchain -> SetBranchStatus("bb.grinch_tdc.hit.yhit",1);
  fchain -> SetBranchAddress("bb.grinch_tdc.hit.yhit",&tdcGHitY);
  fchain -> SetBranchStatus("bb.grinch_tdc.hit.col",1);
  fchain -> SetBranchAddress("bb.grinch_tdc.hit.col",&tdcGHitCol);
  fchain -> SetBranchStatus("bb.grinch_tdc.hit.row",1);
  fchain -> SetBranchAddress("bb.grinch_tdc.hit.row",&tdcGHitRow);

  // fchain -> SetBranchStatus("Ndata.bb.grinch_tdc.tdcelemID",1);
  // fchain->SetBranchAddress("Ndata.bb.grinch_tdc.tdcelemID",&GrinchNum) ;
  // fchain -> SetBranchStatus("bb.grinch_tdc.tdcelemID",1);
  // fchain->SetBranchAddress("bb.grinch_tdc.tdcelemID",&tdcGID) ; // The PMT number 
  // fchain -> SetBranchStatus("bb.grinch_tdc.tdc",1);
  // fchain->SetBranchAddress("bb.grinch_tdc.tdc",&tdcGLe) ; // leading edge
  // fchain -> SetBranchStatus("bb.grinch_tdc.tdc_te",1);
  // fchain -> SetBranchAddress("bb.grinch_tdc.tdc_te",&tdcGTe) ; // trailing edge
  // fchain -> SetBranchStatus("bb.grinch_tdc.tdc_mult",1);
  // fchain->SetBranchAddress("bb.grinch_tdc.tdc_mult",&tdcGMult) ; //tdc multiplicity
  // fchain -> SetBranchStatus("bb.grinch_tdc.tdc_tot",1);
  // fchain->SetBranchAddress("bb.grinch_tdc.tdc_tot",&tdcGTot) ; // tdc time over threshold (te-le)
  
  
  //fchain->SetBranchStatus("bb.grinch_tdc.tdc_mult",1) ;
  //fchain->SetBranchAddress("bb.grinch_tdc.tdc_mult",&tdcGMult) ; //tdc multiplicity
  
  ////I need to go back to a-onl and see how I was looking at all the hits and get that set up in my ifarm   
  ////for when we are looking at all the TDC hits in the window.  grinch_tdc->SetStoreRawHits(kTRUE) in replay_GRINCH
  // fchain->SetBranchStatus("bb.grinch_tdc.hits.TDCelemID", 1); 
  // fchain->SetBranchAddress("bb.grinch_tdc.hits.TDCelemID", &tdcGID); //pmt number
  //
  // fchain->SetBranchStatus("Ndata.bb.grinch_tdc.hits.TDCelemID", 1); 
  // fchain->SetBranchAddress("Ndata.bb.grinch_tdc.hits.TDCelemID", &GrinchNum); 
  //
  // fchain->SetBranchStatus("bb.grinch_tdc.hits.t",1); 
  // fchain->SetBranchAddress("bb.grinch_tdc.hits.t",&tdcGLe); //leading edge for all hits in the window
  //
  // fchain->SetBranchStatus("bb.grinch_tdc.hits.t_te",1);
  // fchain->SetBranchAddress("bb.grinch_tdc.hits.t_te",&tdcGTe);
  //
  // fchain->SetBranchStatus("bb.grinch_tdc.hits.t_tot",1);
  // fchain->SetBranchAddress("bb.grinch_tdc.hits.t_tot",&tdcGTot);


  fchain->SetBranchStatus("Ndata.bb.grinch_adc.adcelemID",1) ;
  fchain->SetBranchAddress("Ndata.bb.grinch_adc.adcelemID",&GrinchADCNum) ; //number of PMTs with ADC channels with a signal in an event. 

  fchain->SetBranchStatus("bb.grinch_adc.adcelemID",1) ; 
  fchain->SetBranchAddress("bb.grinch_adc.adcelemID",&adcGID) ; // The ADC channel 

  fchain->SetBranchStatus("bb.grinch_adc.a_time",1) ; 
  fchain->SetBranchAddress("bb.grinch_adc.a_time",&adcGAtime) ; // the time the first(?)ADC signal went over threshold in the window  

  fchain->SetBranchStatus("bb.grinch_adc.a_amp_p",1) ; 
  fchain->SetBranchAddress("bb.grinch_adc.a_amp_p",&adcGAmp) ; // pedestal-subtracted amplitude

  fchain->SetBranchStatus("bb.grinch_adc.a_mult",1) ;
  fchain->SetBranchAddress("bb.grinch_adc.a_mult",&adcGMult) ;// ADC multiplicity 

  fchain->SetBranchStatus("bb.grinch_adc.a_p",1) ;  
  fchain->SetBranchAddress("bb.grinch_adc.a_p",&grinch_adc) ; // ADC Integral 

  fchain->SetBranchStatus("g.trigbits",1); 
  fchain->SetBranchAddress("g.trigbits",&gTrigBits); // The type of trigger the event was. (1 is the bbcal trig, 16 is the grinch LED) (May not be working)

  fchain->SetBranchStatus("fEvtHdr.fTrigBits",1);
  fchain->SetBranchAddress("fEvtHdr.fTrigBits",&fTrigBits);

  //fchain->SetBranchStatus("fEvtHdr.fRun",1);
  //fchain->SetBranchAddress("fEvtHdr.fRun",&fRun);

  
  
  
 
  //Making histos for the individual channels 
  for (Int_t ig=0;ig<511;ig++) {    
    h_grinch_cluster_tot[ig] = new TH1F(Form("h_grinch_cluster_tot_%d",ig),Form(" ; GRINCH TDC ToT CLUSTER PMT %d ; ",ig),50,0,50);
    //HList.Add(h_grinch_cluster_tot[ig]);
      h_grinch_cluster_le[ig] = new TH1F(Form("h_grinch_cluster_le_%d",ig),Form(" ; GRINCH TDC LE CLUSTER PMT %d ; ",ig),2000,0,2000);
      h_grinch_cluster_le_corr[ig] = new TH1F(Form("h_grinch_cluster_le_corr_%d",ig),Form(" ; GRINCH TDC LE CLUSTER PMT %d ; ",ig),2000,0,2000);
      // HList.Add(h_grinch_cluster_le[ig]);
    h_grinch_cluster_le_tot[ig] = new TH2F(Form("h_grinch_cluster_le_tot_%d",ig),Form("  GRINCH TDC LE vs ToT CLUSTER PMT %d; cluster LE ; cluster ToT; ",ig),40,180,220,50,0,50);
    // HList.Add(h_grinch_cluster_le_tot[ig]);
  }
 
  for (Int_t i=0;i<64;i++){   
    h_grinch_cluster_amp[i] = new TH1F(Form("h_grinch_cluster_amp_%d",i),Form(" ; ADC Amp CLUSTER %d ; ",i), 501, -0.5, 500.5); 
    h_grinch_amp[i] = new TH1F(Form("h_grinch_amp_%d",i),Form(" ; ADC Amp %d ; ",i), 501, -0.5, 500.5);
    //HList.Add(h_grinch_cluster_amp[i]);
  }
    
 
  //// functions that initialize a few things 
  Fill_Search_Stacks(); //sets up cluster finding stacks
  make_row_col_array(); //sets up an array to quickly get the row/col from PMT number
  make_tdc_to_adc_map_array();//sets up arrays to quickly convert between TDC number and ADC number 

  TH2F* h_color_test = new TH2F("h_color_test"," ",16,0,15,16,0,15);
  Double_t color_cnt = 1;
  Double_t normalized  =0;
  for (Int_t i = 1; i<=16; i++){
    for (Int_t j = 1; j<=16; j++){
      normalized = color_cnt * (0.390625); //normalizing 256 to 100. 100/256 = 0.390625
      h_color_test ->SetBinContent(j,i,normalized);
      color_cnt = color_cnt +1;
    }      
  }

  
  Double_t value_x = 0;
  Double_t value_y = 0;

  for (Int_t stepcnt_x = 0; stepcnt_x < nsteps_trx+1; stepcnt_x ++)
    {
      value_x = min_trx + stepsize*stepcnt_x;
      if (value_x >-0.0001 && value_x <0.0001){value_x = 0;}//// it was saving zero as 1.22 E-16, so I hardcoded it to zero. 
      trx_steps_array[stepcnt_x] =value_x;
      //cout<<"trx_steps_array["<<stepcnt_x<<"] = "<<value_x<<endl;
    }
  
  for (Int_t stepcnt_y = 0; stepcnt_y <nsteps_try+1; stepcnt_y ++)
    {
      value_y = min_try + stepsize*stepcnt_y;
      try_steps_array[stepcnt_y] = value_y;
      //cout<< "try_steps_array["<<stepcnt_y <<"] = "<<value_y<<endl;
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
  
  HList.Add(h_grinch_cluster_variance); 
  HList.Add(h_grinch_cluster_tmean);
  
  HList.Add(h_ps_e);  
  HList.Add(h_sh_tot_e); 
  HList.Add(h_ratio_e); 
 
  HList.Add(h_gTrigBits);

  HList.Add(h_grinch_cluster_size);
  HList.Add(h_grinch_cluster_size_cut);
  HList.Add(h_grinch_cluster_size_best);

  HList.Add(h_grinch_cluster_center);
  HList.Add(h_grinch_cluster_center_display);
  HList.Add(test_coord_histo); 
  HList.Add(h_grinch_pmt_good_hit);
  HList.Add(h_grinch_cluster_spread);
  HList.Add(h_grinch_cluster_center_spreadcut);
  
  HList.Add(h_grinch_cluster_adc_tdc_test);

  HList.Add(h_grinch_cluster_ToT_sum);


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
    //HList.Add(h_grinch_hit_le[ig]);  
  }

  TH1F* h_grinch_le[511];
  TH1F* h_grinch_le_LED[511];
  TH1F* h_grinch_tot[511];
  TH1F* h_grinch_hodo_le[511];
  TH2F* h_grinch_le_tot[511];
   for (Int_t ig=0;ig<511;ig++) {
    h_grinch_le[ig] = new TH1F(Form("h_grinch_le_%d",ig),Form(" ; GRINCH TDC LE (ns) PMT %d  ; ",ig),2600,0,2600);
    h_grinch_tot[ig] = new TH1F(Form("h_grinch_tot_%d",ig),Form(" ; GRINCH TDC ToT (ns) PMT %d  ; ",ig),50,0,50);
    // HList.Add(h_grinch_le[ig]);   
    h_grinch_le_tot[ig] = new TH2F(Form("h_grinch_le_tot_%d",ig),Form("  GRINCH TDC LE vs ToT PMT %d; LE ; ToT; ",ig),200,100,300,50,0,50);

    

    h_grinch_hodo_le[ig] = new TH1F(Form("h_grinch_hodo_le_%d",ig),Form(" ; GRINCH TDC LE (ns) -hodo_tmean  PMT %d  ; ",ig),2600,0,2600);
     h_grinch_le_LED[ig] = new TH1F(Form("h_grinch_le_LED_%d",ig),Form(" ; GRINCH TDC LE (ns) LED PMT %d  ; ",ig),2600,0,2600);
  }



  Long64_t nentries = fchain->GetEntries();
  cout<<"nentries: "<<nentries<<endl;

  Int_t temp_cnt=0;

  Int_t zero_count = 0; 
  Int_t clustercnt = 0;
  Int_t max = 0;
  if (entries == -1){
      max = nentries;
    }
  else{
    max = entries;
  }
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
  for (int entry = 0; entry < max ; entry++) {//nentries
    fchain->GetEntry(entry);
    if (entry%10000==0) cout << " Entry = " << entry << endl;

    if (entry ==2){
      //cout< "fRun "<<fRun <<endl;
    }
    
    Double_t tot_e = ps_e+sh_e;
    Double_t rat = tot_e/pmom;
    h_ps_e->Fill(ps_e);
    h_sh_e->Fill(sh_e);
    h_sh_tot_e->Fill(tot_e);
    h_ps_sh_e ->Fill(sh_e, ps_e);    
    h_ratio_e->Fill(rat);

    h_hodo_tmean ->Fill(hodo_tmean);

    h_tr_x ->Fill(tr_x[0]);
    h_tr_y ->Fill(tr_y[0]);

    h_tr_n ->Fill(tr_n);
    h_tr_vz ->Fill(tr_vz[0]);
    h_tr_tg_th ->Fill(tr_tg_th[0]);
    h_tr_tg_ph ->Fill(tr_tg_ph[0]);
    h_tr_th ->Fill(tr_th[0]);
    h_tr_ph ->Fill(tr_ph[0]);
    h_tr_p ->Fill(tr_p[0]);
    h_hcal_e ->Fill(hcal_e);
    h_gem_track_nhits ->Fill(gem_track_nhits);
    h_ps_e_tr_p ->Fill(tr_p[0],ps_e);
    
    h_gTrigBits ->Fill(gTrigBits);
    h_fTrigBits ->Fill(fTrigBits);


    h_grinch_tr_vert_proj ->Fill(tr_x[0] + tr_th[0]*grinch_distance); //small angle approx. tan(th) = sin(th) = th
    h_grinch_tr_horiz_proj ->Fill(tr_y[0] + tr_ph[0]*grinch_distance);
    

    // commented out trigger cut for GEn 
    // if(gTrigBits!=1)// gTrigBits == 1 is the bbcal trigger, == 16 is the grinch LED. 
    //   ///Only want to process bbcal trigger events.
    //   {
    // 	continue; //breaks out and goes to next entry 
    //   }

    Bool_t sh_flag = kFALSE; //flag for when an event passes the shower cut
    Bool_t g_cluster_flag =kFALSE;//flag for when an event has a grinch cluster
    Bool_t good_event_flag = kFALSE;


    //// initialize some arrays
    for(Int_t n = 0; n< 64; n++) //initilalize root index array
      {
	adc_root_index_array[n]=0;
      }
  
       h_kine_W ->Fill(kine_W);
    
    //if (ps_e > 0.2 && abs(rat - 1.15) < 0.1 && abs(tot_e -3.15) < 0.5)// shower cut. abs(tot_e - 3.25)< 0.5 && ps_e > 0.2 && abs(rat - 1.2) < 0.2
    // if(tr_n ==1 && abs(tr_vz[0])<0.08 &&  abs(tr_tg_th[0])<0.15 && abs(tr_tg_ph[0])<0.3 && gem_track_nhits > 3 && tr_p[0] >3.0 && tr_p[0]<4 && hcal_e > 0.025 && ps_e >0.22) //SBS 8 cuts abs(tr_tg_th[0])<0.15 && abs(tr_tg_ph[0])<0.3 
       // if(tr_n ==1 && abs(tr_vz[0])<0.05 && gem_track_nhits > 3 && tr_p[0] >1.4 && tr_p[0]<2.0 && hcal_e > 0.025 && ps_e >0.22  && abs(tr_tg_th[0])<0.1 && abs(tr_tg_ph[0])<0.03 && abs( kine_W- 0.9) < 0.2 ) //SBS 9 cuts && abs(tr_tg_th[0])<0.15 && abs(tr_tg_ph[0])<0.3 && abs(tr_vz[0])<0.08 && abs( kine_W- 0.9) < 0.2 
       //if(tr_n ==1 && abs(tr_vz[0])<0.05 && gem_track_nhits > 3 &&  ps_e < 0.1 &&  abs(tr_tg_th[0])<0.1 && abs(tr_tg_ph[0])<0.03 && sh_e < 1.1 && tr_p[0]<1.4 ) //SBS 9 PIONS could maybe make these tighter.  
       //if (tr_n ==1 && ps_e > 0.22 && gem_track_nhits > 3 && gTrigBits == 4 && abs(tr_vz[0])<0.05 )// gen test cut
       if(tr_n ==1 && ps_e < 0.2  &&  abs(tr_vz[0])<0.05)//
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
      h_grinch_atime ->Fill(adcGAtime[ig]);
      h_grinch_atime_elem ->Fill(gindex_adc,adcGAtime[ig]);

      if (gTrigBits == 16)
	  {
	    h_grinch_amp_elem_LED->Fill(gindex_adc,adcGAmp[ig]);
	    h_grinch_atime_LED ->Fill(adcGAtime[ig]);
	    h_grinch_atime_elem_LED ->Fill(gindex_adc,adcGAtime[ig]);
	  }

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
	//tdcGLe[ig] = tdcGLe[ig] - LE_offset[gindex];// rewriting 
	tdcGLe[ig] = tdcGLe[ig] +200 ;// rewriting 
	h_grinch_le[gindex]->Fill(tdcGLe[ig] - hodo_tmean); // 
	h_grinch_hodo_le[gindex]->Fill(tdcGLe[ig] -hodo_tmean); //
	
	h_grinch_le_tot[gindex] ->Fill(tdcGLe[ig], tdcGTot[ig]);
	h_grinch_le_all ->Fill(tdcGLe[ig]-hodo_tmean);
	h_grinch_tot_elem ->Fill(tdcGID[ig],tdcGTot[ig]);
	h_grinch_tot_all ->Fill(tdcGTot[ig]);
	h_grinch_tot[gindex] ->Fill(tdcGTot[ig]);
	
	grinch_minus_hodo_time = tdcGLe[ig] - hodo_tmean;
	h_grinch_le_elem_hodo ->Fill(tdcGID[ig], grinch_minus_hodo_time);

	if (gTrigBits == 16)
	  {
	    h_grinch_le_elem_LED->Fill(tdcGID[ig],tdcGLe[ig]);
	    h_grinch_le_all_LED ->Fill(tdcGLe[ig]);
	    h_grinch_le_LED[gindex]->Fill(tdcGLe[ig]);
	  }


	// note for eric: I belive this cut is similar to the one in SBSGRINCH. So it would be like we are starting here. 
	if (abs(tdcGLe[ig]-200)<15 && sh_flag && tdcGTot[ig]>0){ // if the event is within the good timing cut for the grinch // was 900 <25 for GMn // was -200
	  goodhit++;
	  ////////////////////////////
	  hit_flag_array[gindex] = 1; // mark that this PMT has a good hit
	  mult_array[gindex] = tdcGMult[ig]; // save it's multiplicity to the array.
	  h_grinch_mult_elem ->Fill(gindex, tdcGMult[ig]);
	  root_index_array[gindex]=ig; // Map the "gindex" (which is the PMT number) to the "ig index" which is what the branch needs.
	                               // I want to be able to go back after finding clusters to look at ADC and TDC and whatnot on those PMTs.   
	  h_grinch_le_elem->Fill(tdcGID[ig],tdcGLe[ig] - hodo_tmean);
	  h_grinch_xhit->Fill(tdcGHitX[ig]);
	  h_grinch_yhit->Fill(tdcGHitY[ig]);

	  // cout<<"xhit " << tdcGHitX[ig]<<endl;
	  // cout<<"yhit "<< tdcGHitY[ig]<<endl;
	  // cout<<endl;

	  //cout<<"row, col from array: "<<row_array[gindex]<<", "<<col_array[gindex]<<endl;
	  //cout<<"row, col from tree: "<<tdcGHitRow[ig]<<" , "<<tdcGHitCol[ig]<<endl;
	  //cout<<endl;

	  h_grinch_pmt_good_hit ->Fill(gindex);
	  
	  good_event_flag = kTRUE;
	  ////////////////////////////
	}
	
     
      }         
     
    }//end loop over TDC PMTs
  

    clusterstack = Clear_Stack(); //clearing the stack we are using for cluster finding for the new event
    
    
   
    if ( good_event_flag){// 
      for(int m = 0; m<510;m++) //fill with 0 
	{
	  sum_array[m]=0;
          check_status[m]=0;
	}

      //// Start of Cluster Finding
      
      cluster_stack[0] = Clear_Stack();
      cluster_stack_cnt=0;
      Int_t top;
      Int_t check_top = 0;
      stack <Int_t> neighbors_stack_temp;
      Int_t neighbors_top = 0;

     for (int pmtcntr = 0; pmtcntr<=509 ; pmtcntr++)//was 9 to 500
      {
	//row = hit_row_array[pmtcntr];
	//col = hit_col_array[pmtcntr];
	neighbors_stack[pmtcntr] = Clear_Stack();
	sum_array[pmtcntr] =  Sum_Adjacent_Hits(pmtcntr, neighbors_stack[pmtcntr]);
	check_status[pmtcntr] = sum_array[pmtcntr];
      }

     
     for (int i = 0;i< 16;i++)
       {
	 cluster_stack[i]=Clear_Stack();
       }
     check_stack = Clear_Stack();

     clusterflag=kFALSE;
     cluster_stack_cnt = 0;
     for (int pmtcntr = 0; pmtcntr<=509 ; pmtcntr++)// was 9 to 500
       {
	 if(sum_array[pmtcntr]>=1 && check_status[pmtcntr]!=-1 )
	   {
	     clusterflag = kTRUE;
	     check_stack.push(pmtcntr);
	     //cout<<"pmtcntr = "<<pmtcntr<<", check status = "<<check_status[pmtcntr]<<" , cluster_stack_cnt = "<<cluster_stack_cnt<<endl;
	     // Print_Stack(cluster_stack[cluster_stack_cnt]);
	     // cout<<"cluster_stack"<<endl;
	     // Print_Stack(check_stack);
	     // cout<<"check_stack"<<endl;

	     while(!check_stack.empty())
	       {
		 // cout<<"11111111111111111111111111111111"<<endl;
		 // Print_Stack(check_stack);
		 // cout<<"check_stack 1"<<endl;

		 check_top = check_stack.top();
		 cluster_stack[cluster_stack_cnt].push(check_top);
		 check_status[check_top] = -1;
		 check_stack.pop();
		 neighbors_stack_temp = neighbors_stack[check_top];
		
		 //cout<< "check_top = "<<check_top<<endl;
		 // Print_Stack(cluster_stack[cluster_stack_cnt]);
		 // cout<<"cluster_stack 1"<<endl;

		 // Print_Stack(neighbors_stack_temp);
		 // cout<<"neighbors_stack_temp 1"<<endl;
		 // cout<<"111111111111111111111111111111"<<endl;
 
		 while(!neighbors_stack_temp.empty())
		   {
		     neighbors_top = neighbors_stack_temp.top();
		     neighbors_stack_temp.pop();
		     if ( check_status[neighbors_top] != -1 && neighbors_top != check_top)
		       {
			 check_stack.push(neighbors_top);
		       }
		     check_status[neighbors_top] = -1;
		   }	
	       }
	     cluster_stack_cnt ++;
	   }

       }

    
 

 g_cluster_flag = kFALSE;

 if(cluster_stack_cnt >0)
   {
     h_grinch_cluster_cnt ->Fill(cluster_stack_cnt);
   }


 for (int i = 0; i < cluster_stack_cnt; i++)
   {
     ClusterSize[i]=Stack_Size(cluster_stack[i]);
     Calculate_Cluster_Center_Row_Col(cluster_stack[i],cluster_horiz[i],cluster_vert[i]);
     Calculate_Cluster_Center_ypos_xpos(cluster_stack[i],cluster_ypos[i],cluster_xpos[i]);
     Find_Cluster_Extrema(cluster_stack[i],horizspread[i],vertspread[i]); // function measures the height and width of the cluster
     Calculate_Cluster_Time(cluster_stack[i],cluster_mean_time[i],cluster_variance[i]);

     if(ClusterSize[i] > 1) 
       {
	 g_cluster_flag = kTRUE; 
	 clustercnt++;
	 Fill_Cluster_Histos(cluster_stack[i], i);

	 for (Int_t stepcnt_x = 0; stepcnt_x < nsteps_trx; stepcnt_x ++)
	   {
	     for (Int_t stepcnt_y = 0; stepcnt_y < nsteps_try+1; stepcnt_y ++)
	       {
		 if( tr_x[0] >= trx_steps_array[stepcnt_x] && tr_x[0] <  trx_steps_array[stepcnt_x+1] && tr_y[0] >= try_steps_array[stepcnt_y] && tr_y[0] <try_steps_array[stepcnt_y+1]){
		   cluster_in_tr_cut_cnt_2D[stepcnt_x][stepcnt_y] ++;
		 }		
	       }
	   }
     
       }
   }

 if(cluster_stack_cnt > 0)
   { 
     //cout<<"SelectionSort"<<endl;
     Int_t tracker[32]= {-1};
     for (Int_t i = 0 ; i<cluster_stack_cnt ; i++)
       {
	 ClusterSize_ordered[i] = ClusterSize[i]; //making a copy of the array so we can leave the unsorted version unchanged 
       }
     // cout<<"Before: "<<endl;
     // for (Int_t i = 0 ; i <cluster_stack_cnt; i++)
     //   {
     // 	 cout<< ClusterSize[i] <<", ";
     //   }
     // cout<<endl;

     SelectionSort(ClusterSize_ordered,tracker,cluster_stack_cnt); //arrays act like reference params

     // cout<< "After: "<<endl;
     // for (Int_t i = 0 ; i <cluster_stack_cnt; i++)
     //   {
     // 	 cout<< ClusterSize_ordered[i] <<", ";
     //   }
     // cout<<endl;
     // cout<<"tracker: "<<endl;
     // for (Int_t i = 0 ; i <cluster_stack_cnt; i++)
     //   {
     // 	 cout<< tracker[i] <<", ";
     //   }
     // cout<<endl;
     // cout<<endl;


     for (Int_t i = 0 ; i < cluster_stack_cnt ; i++)
       {
	 cluster_stack_ordered[i] = cluster_stack[ tracker[i] ];
	 cluster_vert_ordered[i] = cluster_vert[ tracker[i] ];
	 cluster_horiz_ordered[i] = cluster_horiz[ tracker[i] ];
	 cluster_xpos_ordered[i] = cluster_xpos[ tracker[i] ];
	 cluster_ypos_ordered[i] = cluster_ypos[ tracker[i] ];
	 vertspread_ordered[i] = vertspread[ tracker[i] ];
	 horizspread_ordered[i] = horizspread[ tracker[i]];
	 cluster_mean_time_ordered[i] = cluster_mean_time[ tracker[i] ];
	 cluster_variance_ordered[i] = cluster_variance[ tracker[i] ];
       }

     Int_t cluster_array[32];

     struct Cluster cluster_struct_array[32];

      for (Int_t i = 0 ; i < cluster_stack_cnt ; i++)
	{
	  cluster_struct_array[i].size = ClusterSize_ordered[i];
	  cluster_struct_array[i].mean_time = cluster_mean_time_ordered[i];
	  cluster_struct_array[i].sqrt_variance = cluster_variance_ordered[i];
	  cluster_struct_array[i].row = cluster_vert_ordered[i];
	  cluster_struct_array[i].col = cluster_horiz_ordered[i];
	  cluster_struct_array[i].xpos = cluster_xpos_ordered[i];
	  cluster_struct_array[i].ypos = cluster_ypos_ordered[i];
	  cluster_struct_array[i].horizspread = horizspread_ordered[i];
	  cluster_struct_array[i].vertspread = vertspread_ordered[i];
	  
	  StackToArray(cluster_stack_ordered[i],cluster_array);

	  Int_t sum = 0;
	  for ( Int_t j = 0 ; j < ClusterSize_ordered[i] ; j++)
	    {
	      Int_t rootID = root_index_array[cluster_array[j]];
	      cluster_struct_array[i].PMT[j] = cluster_array[j];
	      cluster_struct_array[i].rootID[j] = rootID;
	      cluster_struct_array[i].LE[j] =  tdcGLe[rootID];
	      cluster_struct_array[i].ToT[j] = tdcGTot[rootID];
	      sum = sum + tdcGTot[rootID];
	    }
	  cluster_struct_array[i].ToT_sum = sum;
	}


      Fill_Cluster_Histos_struct(cluster_struct_array, cluster_stack_cnt);


      for (Int_t i = 0 ; i< cluster_stack_cnt ; i++)
	{
	  if( cluster_struct_array[i].xpos > -0.4 && cluster_struct_array[i].xpos < - 0.05) // small slice to stay away from mirror boundaries for a test
	    {
	      h_grinch_cluster_size_cut -> Fill(cluster_struct_array[i].size);
	    }
	}

      
      Double_t x_diff = 0;
      Double_t y_diff = 0;
      Double_t r_diff =0;
      Double_t ri = 0;
      Double_t rj = 0;
      if (cluster_stack_cnt >= 2)
	{
	  for (Int_t i = 0 ; i < cluster_stack_cnt -1 ; i++)
	    {
	      for (Int_t j = i+1 ; j< cluster_stack_cnt ; j++)
		{
		  if( cluster_struct_array[i].xpos > -0.4 && cluster_struct_array[i].xpos < - 0.05 && cluster_struct_array[j].xpos > -0.4 && cluster_struct_array[j].xpos < - 0.05) // small slice to stay away from mirror boundaries for now
		    {

		      ri =  0.5 * 0.0310 *sqrt( pow(cluster_struct_array[i].horizspread,2) + pow(cluster_struct_array[i].vertspread,2));
		      rj =  0.5 * 0.0310 *sqrt( pow(cluster_struct_array[j].horizspread,2) + pow(cluster_struct_array[j].vertspread,2));
		      x_diff = abs(cluster_struct_array[i].xpos - cluster_struct_array[j].xpos);
		      y_diff = abs(cluster_struct_array[i].ypos - cluster_struct_array[j].ypos);
		      r_diff = sqrt( x_diff *x_diff + y_diff*y_diff);
		      h_grinch_cluster_radius_ana ->Fill( ri + rj - r_diff);
		      // if (ri +rj -r_diff < -0.03)//should be only uncorrelated with <-0.03 cut. May have worked? need to think about it more 
		      // 	{
		      // 	  h_grinch_cluster_xdiff ->Fill(x_diff);
		      // 	  h_grinch_cluster_ydiff ->Fill(y_diff);
		      // 	  h_grinch_cluster_xdiff_ydiff ->Fill(y_diff,x_diff);
		      // 	  h_grinch_cluster_rdiff ->Fill(r_diff);
		      // 	}

		      // ri = 0.5* cluster_struct_array[i].vertspread * 0.0310;
		      //rj = 0.5* cluster_struct_array[j].vertspread * 0.0310;

		      if(cluster_struct_array[i].size == 2 && cluster_struct_array[j].size ==2 )// what about only clusters sized 2 to keep things simple
			{
			  h_grinch_cluster_xdiff ->Fill(x_diff);
			  h_grinch_cluster_ydiff ->Fill(y_diff);
			  h_grinch_cluster_xdiff_ydiff ->Fill(y_diff,x_diff);
			  h_grinch_cluster_rdiff ->Fill(r_diff);
			}
		    }
		  
		}
	    }
	}
      

      //for (Int_t i = 0 ; i < cluster_stack_cnt ; i++)
	// {
	//   cout<< "size = "<<cluster_struct_array[i].size <<endl;
	//   cout<<"xpos, ypos = "<<cluster_struct_array[i].xpos<<", "<< cluster_struct_array[i].ypos << endl;
	//   cout<<"row, col = "<< cluster_struct_array[i].row <<", " <<cluster_struct_array[i].col<<endl;
	//   cout<<"horizspread, vertspread = "<<cluster_struct_array[i].horizspread <<", "<<cluster_struct_array[i].vertspread << endl;
	//   cout<<"mean time, variance = "<< cluster_struct_array[i].mean_time <<", "<< cluster_struct_array[i].sqrt_variance << endl; 

	//   cout<< "PMT number: ";
	//   for (Int_t j = 0; j < cluster_struct_array[i].size ; j++)
	//     {
	//       cout<< cluster_struct_array[i].PMT[j]<< ", ";
	//     }
	//   cout<<endl;
	//   cout<< "rootID: ";
	//   for (Int_t j = 0; j < cluster_struct_array[i].size ; j++)
	//     {
	//       cout<< cluster_struct_array[i].rootID[j]<< ", ";
	//     }
	//   cout<<endl;
	//   cout<< "LE: ";
	//   for (Int_t j = 0; j < cluster_struct_array[i].size ; j++)
	//     {
	//       cout<< cluster_struct_array[i].LE[j]<< ", ";
	//     }
	//   cout<<endl;
	//   cout<< "ToT: ";
	//   for (Int_t j = 0; j < cluster_struct_array[i].size ; j++)
	//     {
	//       cout<< cluster_struct_array[i].ToT[j]<< ", ";
	//     }
	//   cout<<endl;
	// }
	//cout<<endl;

     
     
     //projx = tr_x[0] + tr_th[0]*grinch_distance;//small angle approx. tan(th) = sin(th) = th
     // projy = tr_y[0] + tr_ph[0]*grinch_distance;

     projx = tr_x[0] + tan(tr_th[0])*grinch_distance;
     projy = tr_y[0] + tan(tr_ph[0])*grinch_distance;


      for (Int_t i = 0 ; i < cluster_stack_cnt  ; i++)
	{
	  if(cluster_struct_array[i].size >2)
	    {
	      if (projx <= -0.55)
		{
		  h_grinch_mirror_1_projy_ypos ->Fill(projy, cluster_struct_array[i].ypos);
		}
	      if (projx >= -0.4 && projx < -0.1)
		{
		  h_grinch_mirror_2_projy_ypos ->Fill(projy, cluster_struct_array[i].ypos);
		}
	      if (projx > 0.1 && projx < 0.4)
		{
		  h_grinch_mirror_3_projy_ypos ->Fill(projy, cluster_struct_array[i].ypos);
		}
	      if (projx > 0.6)
		{
		  h_grinch_mirror_4_projy_ypos ->Fill(projy, cluster_struct_array[i].ypos);
		}  
	    }
	}

     h_grinch_projx -> Fill(projx);
     h_grinch_projy -> Fill(projy);

     Double_t min_dx = 500;
     Int_t min_tracker = 0;
     Int_t cluster_of_three_or_more_pmts_cnt = 0;
     
     for (Int_t i = 0 ; i < cluster_stack_cnt ; i++)
       {
	 grinch_dx[i] = cluster_xpos_ordered[i] - projx;
	 grinch_dy[i] = cluster_ypos_ordered[i] - projy;

	 if(ClusterSize_ordered[i] >= 2)
	   {
	     h_grinch_cluster_dx ->Fill(grinch_dx[i]);
	     h_grinch_cluster_dy ->Fill(grinch_dy[i]);

	    
	     h_grinch_cluster_projx_xpos ->Fill(projx, cluster_xpos_ordered[i]);
	     h_grinch_cluster_projy_ypos ->Fill(cluster_ypos_ordered[i], projy);
	       
	     
	     cluster_of_three_or_more_pmts_cnt ++;
	     
	     if(abs(grinch_dx[i]) < abs(min_dx))
	       {
		 min_dx = abs(grinch_dx[i]);
		 min_tracker = i;
	       }
	   

	     // cout<< "x pos from cluster: " <<cluster_xpos_ordered[i]<<endl;
	     // cout<< "y pos from cluster: " <<cluster_ypos_ordered[i]<<endl;
	     // cout<<"grinch dx = "<<grinch_dx[i] <<endl;
	     // cout<<"grinch dy = "<<grinch_dy[i] <<endl;
	     // cout<< "cluster size " <<ClusterSize_ordered[i]<<endl;
	   }
       }
    
     for (Int_t i = 0 ; i < cluster_stack_cnt ; i++)
       {
	 if (i != min_tracker && min_dx != 500 && ClusterSize_ordered[i] >2)
	   {
	     h_grinch_cluster_dx_rejected ->Fill(grinch_dx[i]);
	     h_grinch_cluster_dy_rejected ->Fill(grinch_dy[i]);
	     h_grinch_rejected_cluster_index ->Fill(i);
	     h_grinch_cluster_projx_xpos_rejected ->Fill(cluster_xpos_ordered[i], projx);
	     h_grinch_cluster_projy_ypos_rejected ->Fill(cluster_ypos_ordered[i], projy);
	     if(cluster_of_three_or_more_pmts_cnt > 1)
	       {
		 h_grinch_cluster_dx_multi_rejected ->Fill(grinch_dx[i]);
		 h_grinch_cluster_dy_multi_rejected ->Fill(grinch_dy[i]);
		 h_grinch_cluster_projx_xpos_multi_rejected ->Fill(cluster_xpos_ordered[i], projx);
		 h_grinch_cluster_projy_ypos_multi_rejected ->Fill(cluster_ypos_ordered[i], projy);
		 h_grinch_cluster_variance_mulit_rejected ->Fill(cluster_variance_ordered[i]);
	       }
	   }
       }

     if(ClusterSize_ordered[min_tracker] > 2 && min_dx != 500)
       {
	 h_grinch_cluster_dx_best ->Fill(grinch_dx[min_tracker]);
	 h_grinch_cluster_dy_best ->Fill(grinch_dy[min_tracker]);
	 if(cluster_of_three_or_more_pmts_cnt > 1)
	   {
	     h_grinch_best_cluster_index ->Fill(min_tracker);
	     h_grinch_cluster_dx_multi_best ->Fill(grinch_dx[min_tracker]);
	     h_grinch_cluster_dy_multi_best ->Fill(grinch_dy[min_tracker]);
	     h_grinch_cluster_variance_mulit_best ->Fill(cluster_variance_ordered[min_tracker]);
	   }
	 
	 // cout<<"min tracker "<< min_tracker <<endl;
	 // cout<< "min_dx " << min_dx <<endl;
	 // cout<< "cluster size of best cluster " << ClusterSize_ordered[min_tracker]<<endl;
	 // cout<<endl;
       }
	
     if (ClusterSize_ordered[min_tracker]>2 && min_dx != 500) // looking at element with smallest dx  
       {
     	 h_grinch_cluster_size_ordered ->Fill(ClusterSize_ordered[min_tracker]);
     	 h_grinch_cluster_tmean_ordered ->Fill(cluster_mean_time_ordered[min_tracker]);
     	 h_grinch_cluster_variance_ordered ->Fill(cluster_variance_ordered[min_tracker]);
     	 h_grinch_cluster_vert_tmean_ordered ->Fill(cluster_vert_ordered[min_tracker],cluster_mean_time_ordered[min_tracker]);
     	 h_grinch_cluster_horiz_tmean_ordered ->Fill(cluster_horiz_ordered[min_tracker],cluster_mean_time_ordered[min_tracker]);
     	 h_grinch_cluster_vert_variance_ordered ->Fill(cluster_vert_ordered[min_tracker],cluster_variance_ordered[min_tracker]);
     	 h_grinch_cluster_horiz_variance_ordered ->Fill(cluster_horiz_ordered[min_tracker],cluster_variance_ordered[min_tracker]);

     	 h_grinch_cluster_tr_vert_proj ->Fill(projx); 
     	 h_grinch_cluster_tr_horiz_proj ->Fill(projy);
     	 h_grinch_cluster_tr_vert_proj_vs_cluster_vert ->Fill(projx,cluster_vert_ordered[min_tracker]); 
     	 h_grinch_cluster_tr_horiz_proj_vs_cluster_horiz ->Fill(projy,cluster_horiz_ordered[min_tracker]);
	 h_grinch_cluster_projx_xpos_best ->Fill(cluster_xpos_ordered[min_tracker], projx);
	 h_grinch_cluster_projy_ypos_best ->Fill(projy, cluster_ypos_ordered[min_tracker]);
	 h_grinch_cluster_size_best ->Fill(ClusterSize_ordered[min_tracker]);

	  if(cluster_of_three_or_more_pmts_cnt > 1)
	    {
	      h_grinch_cluster_projx_xpos_multi_best ->Fill(cluster_xpos_ordered[min_tracker], projx);
	      h_grinch_cluster_projy_ypos_multi_best ->Fill(cluster_ypos_ordered[min_tracker], projy);

	    }
       }
     
     // if (ClusterSize_ordered[0]>2) // looking at 0th element 
     //   {
     // 	 h_grinch_cluster_size_ordered ->Fill(ClusterSize_ordered[0]);
     // 	 h_grinch_cluster_tmean_ordered ->Fill(cluster_mean_time_ordered[0]);
     // 	 h_grinch_cluster_variance_ordered ->Fill(cluster_variance_ordered[0]);
     // 	 h_grinch_cluster_vert_tmean_ordered ->Fill(cluster_vert_ordered[0],cluster_mean_time_ordered[0]);
     // 	 h_grinch_cluster_horiz_tmean_ordered ->Fill(cluster_horiz_ordered[0],cluster_mean_time_ordered[0]);
     // 	 h_grinch_cluster_vert_variance_ordered ->Fill(cluster_vert_ordered[0],cluster_variance_ordered[0]);
     // 	 h_grinch_cluster_horiz_variance_ordered ->Fill(cluster_horiz_ordered[0],cluster_variance_ordered[0]);

     // 	 h_grinch_cluster_tr_vert_proj ->Fill(tr_x[0] + tr_th[0]*grinch_distance); //small angle approx. tan(th) = sin(th) = th
     // 	 h_grinch_cluster_tr_horiz_proj ->Fill(tr_y[0] + tr_ph[0]*grinch_distance);
     // 	 h_grinch_cluster_tr_vert_proj_vs_cluster_vert ->Fill(tr_x[0] + tr_th[0]*grinch_distance,cluster_vert_ordered[0]); //small angle approx. tan(th) = sin(th) = th
     // 	 h_grinch_cluster_tr_horiz_proj_vs_cluster_horiz ->Fill(tr_y[0] + tr_ph[0]*grinch_distance,cluster_horiz_ordered[0]);

     // 	 Int_t cluster_array[32] ={0};
     // 	 StackToArray(cluster_stack_ordered[0],cluster_array);   
     //   }

     

   }
 

 // if(cluster_stack_cnt > 0){
 //   cout<<"-------------------------"<<endl;
 //   cout <<"Event #"<< entry<< endl;
 //   cout<<"clusters found: "<<cluster_stack_cnt <<endl;
 //   cout<<"-------------------------"<<endl;
 //   cout<< "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " <<endl;
 // }

     
      temp_cnt ++; 

     	
      //// Asking if the cluster had a good shower hit or if there was a good shower hit and no cluster
      if(sh_flag && g_cluster_flag)
	{
	  grinch_hit_sh_hit_cnt++;
	}
      if(sh_flag && !g_cluster_flag)
	{
	  sh_hit_no_grinch_hit_cnt++;
          h_good_electron_no_grinch_hit_track_xy ->Fill(tr_y[0], tr_x[0]);
	  //cout<< "electron cuts passed but no grinch cluster "<<  sh_hit_no_grinch_hit_cnt<<endl;
	  // cout<< "track for that event: tr_y = " <<tr_y[0] <<" tr_x = "<<tr_x[0]<<endl;
	}
      if(!sh_flag && g_cluster_flag)
	{
	  grinch_hit_no_sh_hit_cnt++;
	}
      if(!sh_flag && !g_cluster_flag)
	{
	  no_sh_no_grinch_cnt++;
	}
      
      
    }// end "if good event"
       

  }//end event loop

  Bool_t FITTING = kFALSE;

  if (FITTING)
    {

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
 

      //// NEED to Uncomment 
      //This only works with all the TDC hits. Although it does work for giving the offset gaussian fit a place to start maybe? 

      for (Int_t i =0;i<510;i++)
	{
	  for (Int_t bin = 50; bin<=100 ;bin++)// was 500, 700
	    { 
	      bincontent =  h_grinch_le[i]->GetBinContent(bin);
	      background_sum[i]= background_sum[i] + bincontent;
	      histo_entries = h_grinch_le[i]->GetEntries();
	    }
	  background_fit[i] = new TF1(Form("background_fit_%d",i),"pol0",50,100);// was 500, 700
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
      rate_graph ->GetYaxis() -> SetTitle("Calculated Background rate (kHz)");
 

  

  

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
      Double_t mean_cluster_error[510];
      Double_t sigma_cluster[510];
      Double_t sigma_cluster_error[510];
      Double_t offset_cluster[510];
      Double_t offset_cluster_error[510];
      Double_t sigma_cluster_squared[510];
      Double_t sigma_cluster_squared_error[510];

      TF1* fit_led[510];
      Double_t mean_led[510];
      Double_t mean_led_error[510];
      Double_t sigma_led[510];
      Double_t sigma_led_error[510];
      Double_t offset_led[510];
      Double_t offset_led_error[510];
      Double_t sigma_squared_led[510];
      Double_t sigma_squared_led_error[510];

      Double_t led_cluster_subtract[510];

   
  

      //// FITTING. Need to uncomment after testing 
      for (Int_t i = 0;i<510; i++)
	{
	  //fit[i] = new TF1(Form("fit_%d",i),"gaus",800,1000);     
	  //fit[i] = new TF1(Form("fit_%d",i),offset_gaus,800,1000,4);
	  fit[i] = new TF1(Form("fit_%d",i),offset_gaus,150,250,4);
	  fit[i] ->SetParameters(2*background_per_bin[i],200,5,background_per_bin[i]);//need to do some stats on the histo to give these params?
	  fit[i] ->SetParNames("Amp.","Mean","Sigma","Offset");
	  TFitResultPtr r = h_grinch_le[i] ->Fit(fit[i],"RQ0");
	  Int_t fitStatus = r;
	  if (fitStatus != 0)
	    {
	      //cout<<"fit error for PMT "<< i <<" ,"<< " Status = "<<r<<endl;

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

	  if(mean[i] <1 || mean[i]>400 || mean_error[i] > 10 || sigma[i]>20 || sigma[i] < 1 || sigma_error[i] >10 || offset_error[i]>50)
	    {
	      //cout<<"bad fit in PMT "<< i <<endl;
	      mean[i] = -10;
	      mean_error[i] = 0;
	      sigma[i]= -10;
	      sigma_error[i]= 0;
	      offset[i]= -10;
	      offset_error[i]=0;
	    }
	}


      ////FITTING. Need to uncomment after testing. 
      for (Int_t i = 0; i<510;i++){
	fit_cluster[i] = new TF1(Form("fit_cluster_%d",i),offset_gaus,150,250,4);
	fit_cluster[i] ->SetParameters(background_per_bin[i],200,5,background_per_bin[i]/2);//need to do some stats on the histo to give these params?
	fit_cluster[i] ->SetParNames("Amp.","Mean","Sigma","Offset");
	TFitResultPtr r_cluster = h_grinch_cluster_le[i] ->Fit(fit_cluster[i],"RQ0");
	Int_t fitStatus_cluster = r_cluster;
	if (fitStatus_cluster != 0)
	  {
	    //cout<<"fit error for PMT "<< i <<" ,"<< " Status = "<<r<<endl;

	    mean_cluster[i] = -10;
	    mean_cluster_error[i] = 0;
	    sigma_cluster[i]= -10;
	    sigma_cluster_error[i]= 0;
	    offset_cluster[i]= -10;
	    offset_cluster_error[i]=0;
	    // sigma_squared[i] = -10;
	    //sigma_squared_error[i] = 0; 
	    continue;	  
	  }

	mean_cluster[i] = fit_cluster[i] ->GetParameter(1);
	mean_cluster_error[i] = fit_cluster[i] ->GetParError(1);      
	sigma_cluster[i] = abs(fit_cluster[i] ->GetParameter(2));
	sigma_cluster_error[i] = fit_cluster[i] ->GetParError(2);
	offset_cluster[i] = fit_cluster[i] ->GetParameter(3);
	offset_cluster_error[i] = fit_cluster[i] ->GetParError(3);  
      

	if( abs(mean_cluster[i]-200)>30 || mean_cluster_error[i]>10 || sigma_cluster[i]>20 || sigma_cluster[i] < 1 || sigma_cluster_error[i] >10 || offset_cluster_error[i]>50)
	  {
	    //cout<<"bad fit in PMT "<< i <<endl;
	    mean_cluster[i] = -10;
	    mean_cluster_error[i] = 0;
	    sigma_cluster[i]= -10;
	    sigma_cluster_error[i]= 0;
	    offset_cluster[i]= -10;
	    offset_cluster_error[i]=0;
	  }     
      }


      ////FITTING Need to uncomment after testing.  
      for (Int_t i = 0; i<510;i++)
	{
	  fit_led[i] = new TF1(Form("fit_led_%d",i),offset_gaus,1100,1150,4);
	  fit_led[i] ->SetParameters(100,1120,5,5);//need to do some stats on the histo to give these params?
	  fit_led[i] ->SetParNames("Amp.","Mean","Sigma","Offset");
	  TFitResultPtr r_led = h_grinch_le_LED[i] ->Fit(fit_led[i],"RQ0");
	  Int_t fitStatus_led = r_led;
	  if (fitStatus_led != 0)
	    {
	      //cout<<"fit error for PMT "<< i <<" ,"<< " Status = "<<r<<endl;

	      mean_led[i] = -10;
	      mean_led_error[i] = 0;
	      sigma_led[i]= -10;
	      sigma_led_error[i]= 0;
	      offset_led[i]= -10;
	      offset_led_error[i]=0;
	      // sigma_squared[i] = -10;
	      //sigma_squared_error[i] = 0; 
	      continue;	  
	    }

	  mean_led[i] = fit_led[i] ->GetParameter(1);
	  mean_led_error[i] = fit_led[i] ->GetParError(1);      
	  sigma_led[i] = abs(fit_led[i] ->GetParameter(2));
	  sigma_led_error[i] = fit_led[i] ->GetParError(2);
	  offset_led[i] = fit_led[i] ->GetParameter(3);
	  offset_led_error[i] = fit_led[i] ->GetParError(3);  
      

	  if( abs(mean_led[i]-1120)>50 || mean_led_error[i]>10 || sigma_led[i]>20 || sigma_led[i] < 1 || sigma_led_error[i] >10 || offset_led_error[i]>50)
	    {
	      //cout<<"bad fit in PMT "<< i <<endl;
	      mean_led[i] = -10;
	      mean_led_error[i] = 0;
	      sigma_led[i]= -10;
	      sigma_led_error[i]= 0;
	      offset_led[i]= -10;
	      offset_led_error[i]=0;
	    }
	}


      //// NEED TO UNCOMMENT AFTER TESTING 
      for (Int_t i = 0; i<510; i++)
	{
	  if (mean_led[i] != -10 && mean_cluster[i] != -10)
	    {
	      led_cluster_subtract[i] = mean_led[i] - mean_cluster[i];
	    }
	  else {
	    led_cluster_subtract[i] = -10;
	  }
	}
  

    
 

 



      //// NEED TO UNCOMMENT 
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

    
      // //writing to a file NEED TO UNCOMMENT 
      ofstream fw(Form("/w/halla-scshelf2102/sbs/msatnik/GRINCH_macros/textfiles/offsets_%d_%d.txt",nrun,max), std::ofstream::out);
      if (fw.is_open())
	{
	  fw<<"PMT number"<<"\n";
	  for (int i = 0;i<510;i++)
	    {    
	      if (i%16 ==0 && i!=0) {fw<<"\n";}
	      fw<<i<<" ";
	    }
	  fw << "0 0"<<"\n"; //last 2 channels are for triggers 
	  fw<<"\n";


	  fw<<"LE mean for cluster events"<<"\n";
	  for (int i = 0;i<510;i++)
	    {    
	      if (i%16 ==0 && i!=0) {fw<<"\n";}
	      fw<<round(mean_cluster[i])<<" ";
	    }
	  fw << "0 0"<<"\n"; //last 2 channels are for triggers
	  fw<<"\n";

	  fw<<"LE mean for clusters - 200"<<"\n";
	  for (int i = 0;i<510;i++)
	    {    
	      if (i%16 ==0 && i!=0) {fw<<"\n";}
	      fw<<round(mean_cluster[i]-200)<<" ";
	    }
	  fw << "0 0"<<"\n"; //last 2 channels are for triggers
	  fw<<"\n";

	  fw<<"LE mean for LED events"<<"\n";
	  for (int i = 0;i<510;i++)
	    {    
	      if (i%16 ==0 && i!=0) {fw<<"\n";}
	      fw<<round(mean_led[i])<<" ";
	    }
	  fw << "0 0"<<"\n"; //last 2 channels are for triggers
	  fw<<"\n";

	  fw<<"led mean - 1120  "<<"\n";
	  for (int i = 0;i<510;i++)
	    {    
	      if (i%16 ==0 && i!=0) {fw<<"\n";} 
	      if(mean_led[i]!= -10){
		fw<<round(mean_led[i] - 1120)<<" ";
	      }
	      else  {
		fw<<-100 << " ";
	      }
	    }
	  fw << "0 0"<<"\n"; //last 2 channels are for triggers
	  fw<<"\n";
       
	  fw<<"led mean - cluster mean"<<"\n";
	  for (int i = 0;i<510;i++)
	    {    
	      if (i%16 ==0 && i!=0) {fw<<"\n";} 
	      fw<<round(led_cluster_subtract[i])<<" ";
	    }
	  fw << "0 0"<<"\n"; //last 2 channels are for triggers
	  fw<<"\n";

	  fw.close();
	}
      else cout<<"problem with opening file"<<endl;

  



      Double_t zero[510] = {0}; 
  
      mean_graph = new TGraphErrors(510,pmt,mean,zero,mean_error);
      sigma_graph = new TGraphErrors(510,pmt,sigma,zero,sigma_error);
      offset_graph = new TGraphErrors(510,pmt,offset,zero,offset_error);
 

      mean_graph->SetTitle(Form("GRINCH LE Mean run %d",nrun)); 
      mean_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
      mean_graph ->GetYaxis() -> SetTitle("Mean from Gausian fit");
      sigma_graph->SetTitle(Form("GRINCH LE Sigma run %d",nrun)); 
      sigma_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
      sigma_graph ->GetYaxis() -> SetTitle("Sigma from Gausian fit");
      offset_graph->SetTitle(Form("GRINCH LE Offset run %d",nrun)); 
      offset_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
      offset_graph ->GetYaxis() -> SetTitle("Offset from Gausian fit");

  

      mean_cluster_graph = new TGraphErrors(510,pmt,mean_cluster,zero,mean_cluster_error);
      sigma_cluster_graph = new TGraphErrors(510,pmt,sigma_cluster,zero,sigma_cluster_error);
      offset_cluster_graph = new TGraphErrors(510,pmt,offset_cluster,zero,offset_cluster_error);

      mean_cluster_graph->SetTitle(Form("GRINCH LE Mean CLUSTERS run %d",nrun)); 
      mean_cluster_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
      mean_cluster_graph ->GetYaxis() -> SetTitle("Mean from Gausian fit CLUSTERS");
      sigma_cluster_graph->SetTitle(Form("GRINCH LE Sigma CLUSTERS run %d",nrun)); 
      sigma_cluster_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
      sigma_cluster_graph ->GetYaxis() -> SetTitle("Sigma from Gausian fit CLUSTERS");
      offset_cluster_graph->SetTitle(Form("GRINCH LE Offset CLUSTERS run %d",nrun)); 
      offset_cluster_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
      offset_cluster_graph ->GetYaxis() -> SetTitle("Offset from Gausian fit CLUSTERS");

      mean_led_graph = new TGraphErrors(510,pmt,mean_led,zero,mean_led_error);
      sigma_led_graph = new TGraphErrors(510,pmt,sigma_led,zero,sigma_led_error);
      offset_led_graph = new TGraphErrors(510,pmt,offset_led,zero,offset_led_error);

      mean_led_graph->SetTitle(Form("GRINCH LE Mean LED run %d",nrun)); 
      mean_led_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
      mean_led_graph ->GetYaxis() -> SetTitle("Mean from Gausian fit LED");
      sigma_led_graph->SetTitle(Form("GRINCH LE Sigma LED run %d",nrun)); 
      sigma_led_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
      sigma_led_graph ->GetYaxis() -> SetTitle("Sigma from Gausian fit LED");
      offset_led_graph->SetTitle(Form("GRINCH LE Offset LED run %d",nrun)); 
      offset_led_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
      offset_led_graph ->GetYaxis() -> SetTitle("Offset from Gausian fit LED");
   

      mean_led_cluster_subtract_graph = new TGraphErrors(510,pmt,led_cluster_subtract,zero,zero);
      mean_led_cluster_subtract_graph -> SetTitle(Form("GRINCH LE Mean LED - LE Mean Cluster run %d",nrun)); 
      mean_led_cluster_subtract_graph ->GetXaxis() ->SetTitle("GRINCH PMT");
      mean_led_cluster_subtract_graph ->GetYaxis() -> SetTitle("Mean from LED Gausian fit - Mean from Cluster");
  

    } //END "if (FITTING)"


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
  
  
    
  




  //h_grinch_sh_hit_no_grinch_hit_cnt ->Fill(sh_hit_no_grinch_hit_cnt);
  //h_grinch_grinch_hit_sh_hit_cnt ->Fill(grinch_hit_sh_hit_cnt);
  // h_grinch_grinch_hit_no_sh_hit_cnt ->Fill(grinch_hit_no_sh_hit_cnt);
  //h_grinch_no_sh_no_grinch_cnt ->Fill(no_sh_no_grinch_cnt);
  // h_grinch_sh_ps_hit_cnt ->Fill(sh_ps_hit_cnt);

  Double_t grinch_hit_double = grinch_hit_sh_hit_cnt;
  Double_t sh_hit_double = sh_ps_hit_cnt;
  Double_t percent_overall = 100 * grinch_hit_double/sh_hit_double; 

  Double_t cluster_in_tr_cut_cnt_double = cluster_in_tr_cut_cnt;
  Double_t tr_cut_cnt_double = tr_cut_cnt ; 
  Double_t tr_cut_percent = 100 * cluster_in_tr_cut_cnt_double/tr_cut_cnt_double;


  h_grinch_efficiency_percent->Fill(percent_overall);

 

  Double_t cluster_in_tr_cut_cnt_double_2D[500][100] = {0};
  Double_t tr_cut_cnt_double_2D[500][100] = {0};
  Double_t tr_cut_percent_2D[500][100] = {0};
 
  for (Int_t stepcnt_x = 0; stepcnt_x < nsteps_trx; stepcnt_x ++)
    {
      for (Int_t stepcnt_y = 0; stepcnt_y < nsteps_try; stepcnt_y ++)
	{
	  Double_t percent = 0;
	  cluster_in_tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] = cluster_in_tr_cut_cnt_2D[stepcnt_x][stepcnt_y];
	  tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] = tr_cut_cnt_2D[stepcnt_x][stepcnt_y];
	  if(tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] == 0 && cluster_in_tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] == 0){
	    percent = 0;
	     h_grinch_eff_xy ->SetBinContent(stepcnt_y+1,stepcnt_x+1, 0);
	  }
	  else  if( tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] != 0 && cluster_in_tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y] ==0)
	    {
	     percent = 0;
	      h_grinch_eff_xy ->SetBinContent(stepcnt_y+1,stepcnt_x+1, 0);
	    }
	  else{
	    percent = 100 * cluster_in_tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y]/tr_cut_cnt_double_2D[stepcnt_x][stepcnt_y];
	     h_grinch_eff_xy ->SetBinContent(stepcnt_y+1,stepcnt_x+1, percent);
	  }	  
	  tr_cut_percent_2D[stepcnt_x][stepcnt_y] = percent;	
	  if ( tr_cut_percent_2D[stepcnt_x][stepcnt_y]  !=0){
	    // cout<<"["<<trx_steps_array[stepcnt_x]<<"] ["<<try_steps_array[stepcnt_y]<<"]       "<< cluster_in_tr_cut_cnt_2D[stepcnt_x][stepcnt_y] <<" / "<< tr_cut_cnt_2D[stepcnt_x][stepcnt_y]  << " = "<< tr_cut_percent_2D[stepcnt_x][stepcnt_y]  << " %"<<endl; //uncomment for efficiences 
	  }
	}		
    }
 
  h_grinch_eff_xy ->SetContour(255);

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

  cout<< "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;

  cout<<"Effiency on Good Hits: "<< grinch_hit_sh_hit_cnt<< "/" << sh_ps_hit_cnt <<" = " << percent_overall <<" %"<<endl; 

  //cout<<"Effiency on the small track region we cut on for a test: " << cluster_in_tr_cut_cnt<<" / "<< tr_cut_cnt << " = "<< tr_cut_percent << " %"<<endl;

  cout<<" written to file "<<outputhist<<endl;
  cout<<"cluster analysis on "<<max<<" events"<<endl;
  cout<<clustercnt<<" clusters found"<<endl;
  cout<<"hit_sum = "<< hit_sum << endl;

  cout<<"good shower energy on "<< sh_ps_hit_cnt<<" events" <<endl;
  cout<<"good shower energy and NO good grinch cluster on "<< sh_hit_no_grinch_hit_cnt<<" events" <<endl;
  cout<<"good shower energy and good grinch cluster on "<< grinch_hit_sh_hit_cnt <<" events" <<endl;
  //cout<<"NO good shower energy and good grinch cluster on "<< grinch_hit_no_sh_hit_cnt <<" events" <<endl;
  //cout<<"NO good shower energy and NO grinch cluster on "<< no_sh_no_grinch_cnt <<" events" <<endl;

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



Int_t Sum_Adjacent_Hits(Int_t PMT, stack <Int_t>& neighbors_stack)
{ //// Finds how many adjacent PMTs also had a hit and also returns the list via reference parameter as a stack 
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

  if (row == 0)//top row of the detector
    {
      //cout<<"top row"<<endl;
      if (PMT == 0)
	{
	  searchstack = top_row_left_search;
	  // cout<<"PMT 0"<<endl;
	}
      else if (PMT == 7)
	{
	  searchstack = top_row_right_search;
	  //cout<<"PMT 7"<<endl;
	}
      else 
	{
	  searchstack = top_row_middle_search;
	  //cout<<"PMT "<<PMT<<endl;
	}
      //Print_Stack(searchstack);
    }

  else if (row == 59)// bottom row of the detector
    {
      if (PMT == 501)
	{
	  searchstack = bottom_row_left_search;
	}
      else if (PMT == 509)
	{
	  searchstack = bottom_row_right_search;
	}
      else 
	{
	  searchstack = bottom_row_middle_search;
	}
    }

  else if (row%2 == 0)//short row of 8
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
      if (col == 0 )
	{
	  searchstack = long_row_on_left_search;
	}
      else if(col == 8)
	{
	  searchstack = long_row_on_right_search;
	}
      else //not on an edge
	{
	  //cout<<"long row in the middle"<<endl;
	  searchstack = middle_search;
	}
    }

  //Print_Stack(searchstack);
  //cout<<"searchstack"<<endl;

  stack <Int_t> neighborstack = Add_to_Stack(searchstack,PMT); //add "PMT" integer to every element of the stack 
  //Print_Stack(neighborstack);
  //cout<<"neighborstack"<<endl;
  
  while (!neighborstack.empty())
    {
      neighbor = neighborstack.top();
      hit_sum = hit_sum + hit_flag_array[neighbor];   
      if(hit_flag_array[neighbor] == 1)
	{
	  neighbors_stack.push(neighbor);
	}
      neighborstack.pop();
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
  //// If we are on the one of the edges of the PMT array, we need to adjust slighty. 
  ////
  //// short_row_on_left : exclude [-1]
  //// short_row_on_right: exclude [+1]
  //// top_row_middle : exclue [-8], [-9]
  //// top_row_left : exclude [-1], [-8], [-9] this is only for PMT 0
  //// top_row_right: exclude [-9] ,[-8], [+1] this is only for PMT 7
  //// bottom_row_middle : exclude [+8], [+9] 
  //// long_row_on_right : exclude [-8], [+1], [+9]
  //// long_row_on_left : exclude [-9], [-1], [+8]
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


void Calculate_PMT_Coord(double& horiz, double& vert, Int_t PMT) 
// this is probably uncessary as these should already be in the database.
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


void Calculate_Cluster_Center_Row_Col(stack <Int_t> inputstack, double& horiz, double& vert)
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


void Calculate_Cluster_Center_ypos_xpos(stack <Int_t> inputstack, double& ypos, double& xpos)
{//// Takes the average of the horizontal and vertical positions of the PMTs in the cluster. 
 //// The result is returned via reference parameter. 
  Double_t xpos_sum=0;
  Double_t ypos_sum=0;
  Double_t counter=0;
  Double_t temp_xpos = 0;
  Double_t temp_ypos = 0;
  Int_t PMT=0;

  while (!inputstack.empty())
    {
      PMT = inputstack.top();
      temp_xpos = xpos_array[PMT];
      temp_ypos = ypos_array[PMT];
      xpos_sum = xpos_sum + temp_xpos;
      ypos_sum = ypos_sum + temp_ypos;
      counter = counter +1;
      inputstack.pop();
    }
  
  xpos = xpos_sum/counter;
  ypos = ypos_sum/counter;  
} //end calculate cluster center


void Calculate_Cluster_Time(stack <Int_t> inputstack, Double_t& mean, Double_t& variance)
{ //// Calculates the mean and variance of the Leading Edge (LE) of the PMTs in a cluster
  //// (we may want to consider using something other than mean and variance since it's not necessarily a normal distribution)
  stack <Int_t> inputstackcopy = inputstack; 
  Double_t time_sum=0;
  Double_t N = Stack_Size(inputstack);
  Double_t variance_squared_sum = 0;
  Double_t variance_squared = 0;
  Int_t PMT;
  Int_t root_id;
  Double_t var_sq=0;
  Double_t diff = 0;

  while(!inputstack.empty())
    {
      PMT = inputstack.top();
      inputstack.pop();
      root_id = root_index_array[PMT];//get the index that the branches need in order to look at the specific PMT for this event
      time_sum = time_sum + tdcGLe[root_id];
      //cout<<"[PMT "<<PMT<<"] LE = "<< grinch_tdc_le[root_id] <<endl;
    }
  mean = time_sum/N;
  //cout<<"mean = "<<mean<<endl;

  while(!inputstackcopy.empty())
    {
      PMT = inputstackcopy.top();
      inputstackcopy.pop();
      root_id = root_index_array[PMT];//get the index that the branches need in order to look at the specific PMT for this event
      diff = tdcGLe[root_id] -mean;
      var_sq = TMath::Power(diff,2);
      variance_squared_sum = variance_squared_sum + var_sq;
      //cout<<"[PMT "<<PMT<< "] diff = "<<diff<< ", var sq = " << var_sq <<endl;
    }
  variance_squared = variance_squared_sum / N;
  // cout<<"variance_squared = " <<variance_squared<<endl;
  variance = TMath::Sqrt(variance_squared); 
  //cout<<"variance= "<<variance<<endl;

}//end Calculate_Cluster_Time

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


void Fill_Cluster_Histos_struct(struct Cluster cluster_struct_array[], Int_t entries)
{ //// Filling some histograms for the cluster using the struct info. 
  //// This is just to keep things cleaner in the main program. 
  Int_t pmt;
  Int_t LE;
  Int_t ToT;

  Double_t projx = tr_x[0] + tan(tr_th[0])*grinch_distance;
  Double_t projy = tr_y[0] + tan(tr_ph[0])*grinch_distance;
  // && abs(cluster_struct_array[i].xpos - projx)<0.15
 
  for (Int_t i = 0; i < entries ; i++)
    {
      if (cluster_struct_array[i].size > 1 && abs(cluster_struct_array[i].xpos - projx)<0.15){
	h_grinch_cluster_tmean -> Fill(cluster_struct_array[i].mean_time);
	h_grinch_cluster_variance ->Fill(cluster_struct_array[i].sqrt_variance);
	h_grinch_cluster_size -> Fill(cluster_struct_array[i].size);
	h_grinch_cluster_center_xpos ->Fill(cluster_struct_array[i].xpos);
	h_grinch_cluster_center_ypos -> Fill(cluster_struct_array[i].ypos);
	h_grinch_cluster_center ->Fill(cluster_struct_array[i].col,cluster_struct_array[i].col);
	h_grinch_cluster_spread->Fill(cluster_struct_array[i].horizspread,cluster_struct_array[i].vertspread);
	h_grinch_cluster_ToT_sum ->Fill(cluster_struct_array[i].ToT_sum);
	h_grinch_cluster_ToT_sum_vs_size -> Fill(cluster_struct_array[i].ToT_sum, cluster_struct_array[i].size);
	h_grinch_cluster_tmean_size ->Fill(cluster_struct_array[i].size, cluster_struct_array[i].mean_time);
	//plot vs prshower energy?

	for (Int_t j = 0 ; j < cluster_struct_array[i].size ; j++)
	  {
	    pmt = cluster_struct_array[i].PMT[j];
	    LE = cluster_struct_array[i].LE[j];
	    ToT = cluster_struct_array[i].ToT[j]; 

	    h_grinch_cluster_elem -> Fill(pmt);
	    h_grinch_cluster_le_elem ->Fill(pmt, LE);
	    h_grinch_cluster_le_elem_corr ->Fill(pmt, LE - LE_offset[pmt]);
	    h_grinch_cluster_tot_elem ->Fill(pmt, ToT);
	    h_grinch_cluster_le[pmt]->Fill(LE);
	    h_grinch_cluster_le_corr[pmt]->Fill(LE - LE_offset[pmt]);
	    h_grinch_cluster_tot[pmt]->Fill(ToT);
	    h_grinch_cluster_le_all -> Fill(LE);
	    h_grinch_cluster_le_all_corr -> Fill(LE -LE_offset[pmt]);
	    h_grinch_cluster_tot_all ->Fill(ToT);
	    h_grinch_cluster_le_tot[pmt] ->Fill(LE, ToT);
	  }
      }
    }
}// end Fill_Cluster_Histos_struct



void Fill_Cluster_Histos(stack <Int_t> inputstack, Int_t index)
{ //// Fills the hitograms for the cluser in this function to keep the main loop clearer. 
  //// Several of these are now in the "Fill_Cluster_Histos_struct" funcition instead. 
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
      // h_grinch_cluster_le_elem -> Fill(PMT,tdcGLe[root_id] - hodo_tmean);
      // h_grinch_cluster_te_elem -> Fill(PMT,tdcGTe[root_id]);
      //h_grinch_cluster_tot_elem -> Fill(PMT,tdcGTot[root_id]);
      h_grinch_cluster_mult_elem -> Fill(PMT,tdcGMult[root_id]);
      //h_grinch_cluster_tot[PMT]->Fill(tdcGTot[root_id]);
      // h_grinch_cluster_le[PMT]->Fill(tdcGLe[root_id]-hodo_tmean);
      // h_grinch_cluster_le_all -> Fill(tdcGLe[root_id]-hodo_tmean);
      // h_grinch_cluster_tot_all ->Fill(tdcGTot[root_id]);
      //h_grinch_cluster_elem ->Fill(PMT);
      //h_grinch_cluster_le_tot[PMT] ->Fill(tdcGLe[root_id], tdcGTot[root_id]);

      grinch_minus_hodo_time = tdcGLe[root_id] - hodo_tmean;
      h_grinch_cluster_le_elem_hodo ->Fill(PMT, grinch_minus_hodo_time);
      
	 // Print_Stack(cluster_stack[i]);
	 // cout<<"cluster size = "<< ClusterSize[i] <<endl;
	 // cout<<"cluster center is at "<<cluster_horiz[i]<<", "<<cluster_vert[i]<<endl;
	 // cout<<"cluster center is at "<<cluster_ypos[i]<<", "<<cluster_xpos[i]<<endl;
	 // cout<<"cluster is "<<horizspread[i]<<" wide, "<< vertspread[i]<<" tall"<<endl;


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

  //h_grinch_cluster_center ->Fill(cluster_horiz[index],cluster_vert[index]);
  h_grinch_cluster_center_display ->Fill(cluster_horiz[index],-cluster_vert[index]);

  //h_grinch_cluster_center_xpos ->Fill(cluster_xpos[index]);
  // h_grinch_cluster_center_ypos ->Fill(cluster_ypos[index]);

  h_grinch_cluster_center_horiz_tr->Fill(cluster_horiz[index], tr_y[0]); //horizontal
  h_grinch_cluster_center_horiz_tr_display->Fill(-cluster_horiz[index], tr_y[0]); //horizontal with neg clusterx
  h_grinch_cluster_center_vert_tr ->Fill(cluster_vert[index], tr_x[0]); //vertical 
      
  //h_grinch_cluster_spread->Fill(horizspread[index],vertspread[index]);
  h_grinch_cluster_track_xy ->Fill(tr_y[0], tr_x[0]);  
	
  //h_grinch_cluster_tmean ->Fill(cluster_mean_time[index]);
  // h_grinch_cluster_variance ->Fill(cluster_variance[index]);
  h_grinch_cluster_vert_tmean ->Fill(cluster_vert[index],cluster_mean_time[index]);
  h_grinch_cluster_horiz_tmean ->Fill(cluster_horiz[index],cluster_mean_time[index]);
  h_grinch_cluster_vert_variance ->Fill(cluster_vert[index],cluster_variance[index]);
  h_grinch_cluster_horiz_variance ->Fill(cluster_horiz[index],cluster_variance[index]);

  //h_grinch_cluster_size ->Fill(ClusterSize[index]);
}////end Fill_Cluster_Histos




void SelectionSort(Int_t a[], Int_t tracker[], Int_t n)
{ //// Uses a classic selection sort algorithm to sort from largest to smallest 
  //// tracker[] shows the final positions of the elements: 
  //// Example: Array={20,23,22} becomes {23,22,20} with tracker={1,2,0}
  for (Int_t m =0 ; m < n ; m++)
    {
      tracker[m] = m;
    }

  for (Int_t i = 0; i < n; i++)
    {
      Int_t max_index = i;
      Int_t max_element = a[i];
      
      for (Int_t j = i+1; j < n ;j++)
	{
	  if (a[j] > max_element)
	    {
	      max_element = a[j];
	      max_index = j;
	    }
	}
      Int_t temp = a[i];
      a[i] = a[max_index];
      a[max_index] = temp;
      temp = tracker[i];
      tracker[i] = tracker[max_index];
      tracker[max_index] = temp;
    }
}//end SelectionSort


void StackToArray(stack <Int_t> s, Int_t a[])
{
  stack <Int_t> scopy= s;
  Int_t i = 0;
  while (!scopy.empty())
    {
      a[i] = scopy.top();
      scopy.pop();
      i++;
    }
}

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
	  xpos_array[pmt] = -0.9145 + 0.0310*i; 
	  if(N_COL == 8)//short row
	    {
	      ypos_array[pmt]= 0.1085 - 0.0310*j;
	    }
	  else
	    {//long row
	      ypos_array[pmt]= 0.1085 + 0.0155 - 0.0310*j; // offset by half of a PMT width
	    }
	  pmt ++;
	}
    }
  
    // for (Int_t r = 0; r< 510; r++)
    // {
    //   cout<< "PMT "<<r<<": ROW "<< row_array[r] <<", COL"<< col_array[r]<<endl;
    //   cout<< "PMT "<<r<<": xpos "<<xpos_array[r]<<", ypos "<< ypos_array[r]<<endl;
    // }
  
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



void Palette1()
{
  Int_t n =256;
  static Int_t colors[256];
  static Bool_t initialized = kFALSE;
  
  Double_t Red[5] =  {0.00, 0.00, 0.00, 1.00, 1.00};
  Double_t Green[5] =  {0.00, 0.00, 1.00, 1.00, 1.00};
  Double_t Blue[5] =  {0.00,1.00, 0.00, 0.00, 1.00};
  Double_t Length[5] = {0.00,0.8,0.9, 0.95, 1.00};

  if (!initialized)
    {
      Int_t FI = TColor::CreateGradientColorTable(5,Length,Red,Green,Blue,256);
      for (int i = 0;i<n;i++) colors [i] = FI+i;
      initialized = kTRUE;
      return;
    }
  gStyle ->SetPalette(n,colors);
}


void set_color_env()
{ 	
  //Creates a nice color gradient.
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t stops[NRGBs] = { 0.00, 0.45, 0.90, 0.95, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);   //Chooses custom gradient NCont.
}
