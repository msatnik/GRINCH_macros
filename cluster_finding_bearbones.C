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
//// cluster_finding_barebones.C
//// Maria Satik
//// Feb 2023
//// msatnik@wm.edu
////
//// This program analyzes GRINCH data and finds clusters of 2 or more PMTs in trigger events. 
//// As of this moment (March 3, 2023), it is using electron cuts for GEn data taken January 2023. 
//// This program finds clusters by utilizing stacks to simulate recursion.
//// The data for each event is orgainzed into an array of "Cluster" structs. Where the largest cluster is the [0]th element.
////        cluster_struct_array[0].size would be the number of PMTs fired in the largest cluster for that event.  
////        cluster_struct_array[0].PMT[0] would be the PMT number of the [0]th PMT to fire in the largest cluster (no particualr PMT order). 
////        Scroll down to just below this header to see all the memebers. 
////        These get re-written for every event and need to be written into the tree instead. 
//// 
//// 
//// Running the program: Easiest way is to run from an a-onl machine. ssh into there and in the terminal type "gogrinch". 
//// The cofig file to put the rootfile paths in is located in /adaqfs/home/a-onl/sbs/Grinch_replay/grinch/macros/RunLists.
//// Here is an example where I have the rootfiles listed in the config file "run_list_grinch_3038.txt" in that location for 50000 events
////
////        gogrinch 
////        cd grinch/macros/
////        analyzer
////        .x cluster_finding_barebones.C("grinch",3038,"run_list_grinch_3038.txt",50000) 
////        .q (when you are all done)
//// 
//// This runs the program for that many events and fills lots of histograms. A few that are possibly useful: 
////        h_grinch_cluster_le_elem -> Draw("colz")
////        h_grinch_cluster_le_all ->Draw()
////        h_grinch_cluster_tot_all ->Draw()
////        h_grinch_projx_xpos ->Draw("colz")
////        h_grinch_cluster_cnt ->Draw()
//// The program cluster_finding_hit.C does everything this program does but also has options for fitting. 
/////////////////////////////////////////////////////

//// struct to hold the final cluster information for each cluster in each event. 
struct Cluster{
  Int_t size = 0; // number of PMTs in the cluster
  Int_t PMT[32] = {0}; // PMT number for each PMT in cluster
  Int_t rootID[32] = {0}; // index the tree needs for each PMT
  Double_t LE[32] = {0}; // Leading Edge: time the singal crosses the threshold (relative to trigger) for each PMT
  Double_t ToT[32] = {0}; // Time Over Threshold = Leading Edge - Trailing Edge. "Size" of the signal for each PMT
  Double_t mean_time = 0; // the average of the LE of all the PMTs in the cluster. 
  Double_t sqrt_variance = 0; // the square root of the calculated variance of the PMTs in the cluster. (This is the standard deviation if you assume a normal distribution). 
  Double_t row = 0; // Calcuated row, col of the cluster. (It is just an average of the location of each PMT in the cluster)
  Double_t col = 0; 
  Double_t xpos = 0; // Same as row, col calculation but are in terms of the trasport coordinate system given in the database. 
  Double_t ypos = 0;
  Double_t horizspread = 0; // how many columns the cluster occupies. (How wide the cluster is)
  Double_t vertspread = 0; // how many rows the cluster occupies. (How tall the cluster is). 
};

////global variables
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

Int_t row;
Int_t col;


const Int_t N_ROW=60;

const Double_t grinch_distance = 0.48; //distance from focal plane. (Stolen from Andrew's code. May be different for different kinematics). 

Int_t row_array[510];
Int_t col_array[510];
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
Double_t tdcGMult[1000];
Double_t tdcGHitX[1000];
Double_t tdcGHitY[1000];
Double_t tdcGHitRow[1000];
Double_t tdcGHitCol[1000];

Double_t kine_W; 

Double_t hodo_tmean;



Double_t gTrigBits; // Should be "1" if it is a bbcal trigger, 16 if it is a GRINCH LED event (I think)

 
Int_t hit_flag_array[510]; // array to help keep track of which PMTs were hit during an event
Int_t mult_array[510]; // array for the multiplicities of each PMT during an event
Int_t sum_array[510]; // array to help with cluster finding: records how many adjacent PMTs have a hit for each PMT
Int_t root_index_array[510]; // An array to help us convert betweeen the PMT number and the index that root uses for an event

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

TH1F* h_grinch_hitslistlength = new TH1F("h_grinch_hitslistlength",";total fires in the multi-hit tdc;", 500,0,500);

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



TH1F* h_grinch_clustercnt = new TH1F("h_grinch_clustercnt","number of clusters in event with more than 2 pmts",10,0,10);

TH1F* h_grinch_cluster_tr_vert_proj =  new TH1F("h_grinch_cluster_tr_vert_proj","track x projected to grinch window", 200,-2,2);
TH1F* h_grinch_cluster_tr_horiz_proj = new TH1F("h_grinch_cluster_tr_horiz_proj","track y projected to grinch window",200,-2,2);

TH1F* h_grinch_cluster_center_xpos = new TH1F("h_grinch_cluster_center_xpos","cluster x at grinch window",100,-1,1);
TH1F* h_grinch_cluster_center_ypos = new TH1F("h_grinch_cluster_center_ypos","cluster y at grinch window",100,-1,1);

TH1F* h_grinch_projx =  new TH1F("h_grinch_projx","track x projected to grinch window",100,-1,1);
TH1F* h_grinch_projy =  new TH1F("h_grinch_projy","track y projected to grinch window",100,-1,1);


TH2F* h_grinch_cluster_tr_vert_proj_vs_cluster_vert =  new TH2F("h_grinch_cluster_tr_vert_proj_vs_cluster_vert",";track x projected to grinch window ; cluster vertical", 100,-1,1,59,0,60);
TH2F* h_grinch_cluster_tr_horiz_proj_vs_cluster_horiz = new TH2F("h_grinch_cluster_tr_horiz_proj_vs_cluster_horiz",";track y projected to grinch window ; cluster horizontal",100,-1,1,17,0,9);

TH2F* h_grinch_cluster_projx_xpos_best = new TH2F("h_grinch_cluster_projx_xpos_best","; cluster x position best cluster ; projected x at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projy_ypos_best = new TH2F("h_grinch_cluster_projy_ypos_best","; cluster y position best cluster ; projected y at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projx_xpos = new TH2F("h_grinch_cluster_projx_xpos","; cluster x position ; projected x at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projy_ypos = new TH2F("h_grinch_cluster_projy_ypos","; cluster y position ; projected y at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projx_xpos_rejected = new TH2F("h_grinch_cluster_projx_xpos_rejected","; cluster x position rejected clusters ; projected x at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projy_ypos_rejected = new TH2F("h_grinch_cluster_projy_ypos_rejected","; cluster y position rejected clusters ; projected y at grinch window from track",100,-1,1,100,-1,1);

TH2F* h_grinch_cluster_projx_xpos_multi_best = new TH2F("h_grinch_cluster_projx_xpos_multi_best","; cluster x position best cluster ; projected x at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projy_ypos_multi_best = new TH2F("h_grinch_cluster_projy_ypos_multi_best","; cluster y position best cluster ; projected y at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projx_xpos_multi_rejected = new TH2F("h_grinch_cluster_projx_xpos_multi_rejected","; cluster x position rejected clusters ; projected x at grinch window from track",100,-1,1,100,-1,1);
TH2F* h_grinch_cluster_projy_ypos_multi_rejected = new TH2F("h_grinch_cluster_projy_ypos_multi_rejected","; cluster y position rejected clusters ; projected y at grinch window from track",100,-1,1,100,-1,1);


TH2F* h_grinch_cluster_center_vert_tr = new TH2F("h_grinch_cluster_center_vert_tr",";GRINCH cluster center vertical ; track x;", 120,0,60,250,-1,1);
TH2F* h_grinch_cluster_center_horiz_tr = new TH2F("h_grinch_cluster_center_horiz_tr",";GRINCH cluster center horizontal; track y;", 18,0,9,75,-0.3,0.3);
TH2F* h_grinch_cluster_center_horiz_tr_display = new TH2F("h_grinch_cluster_center_horiz_tr_display",";cluster center horizontal; bb.tr.y[0];", 18,-9,0,75,-0.3,0.3);


TH2F* h_grinch_cluster_track_xy = new TH2F("h_grinch_cluster_track_xy",";bb.tr.y[0];bb.tr.x[0];", 60,-0.3,0.3,200,-1,1);

TH2F* h_good_event_track_xy = new TH2F("h_good_event_track_xy",";bb.tr.y[0];bb.tr.x[0];",  60,-0.3,0.3,200,-1,1);


TH2F* h_grinch_eff_xy = new TH2F("h_grinch_eff_xy",";bb.tr.y[0];bb.tr.x[0];",  32,-0.2,0.2,112,-0.7,0.7);


TH1F* h_grinch_cluster_size = new TH1F("h_grinch_cluster_size" ,"; Cluster Size ",18,2,20);
TH2F* h_grinch_cluster_center = new TH2F("h_grinch_cluster_center", "; horizontal ; vertical ;",16,0,8,60,0,60);
TH2F* h_grinch_cluster_center_display = new TH2F("h_grinch_cluster_center_display", "; horizontal ; vertical ;",16,0,8,60,-600);

TH1F* h_grinch_cluster_tmean= new TH1F("h_grinch_cluster_tmean",";Cluster Mean Time (LE)", 40,180,220);
TH1F* h_grinch_cluster_variance = new TH1F("h_grinch_cluster_variance","; Cluster Variance (LE)",40,0,20);
TH2F* h_grinch_cluster_vert_tmean =new TH2F("h_grinch_cluster_vert_tmean", ";cluster vert; cluster mean time (LE)",60,0,60,40,180,220);
TH2F* h_grinch_cluster_horiz_tmean = new TH2F("h_grinch_cluster_horiz_tmean", ";cluster horiz; cluster mean time (LE)",16,0,8,40,180,220);
TH2F* h_grinch_cluster_vert_variance  =new TH2F("h_grinch_cluster_vert_variance", ";cluster vert; cluster variance (LE)",60,0,60,40,0,20);
TH2F* h_grinch_cluster_horiz_variance = new TH2F("h_grinch_cluster_horiz_variance", ";cluster horiz; cluster variance (LE)",16,0,8,40,0,20);

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
TH2F* h_grinch_cluster_spread = new TH2F("h_grinch_cluster_spread",";cluster width;  cluster height;",8,1,9,59,1,60);


TH1F* h_grinch_cluster_elem = new TH1F("h_grinch_cluster_elem",";GRINCH TDC elemID;", 510,0,510);
TH2F* h_grinch_cluster_le_elem = new TH2F("h_grinch_cluster_le_elem","; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,200,0,400);

TH2F* h_grinch_cluster_le_elem_hodo = new TH2F("h_grinch_cluster_le_elem_hodo","; GRINCH TDC elemID ; GRINCH TDC LE (ns) - hodo tmean ",510,0,510,200,800,1000);

TH2F* h_grinch_cluster_tot_elem = new TH2F("h_grinch_cluster_tot_elem"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) ",510,0,510,110,-10,100);
TH2F* h_grinch_cluster_mult_elem = new TH2F("h_grinch_cluster_mult_elem"," ; GRINCH TDC elemID ; GRINCH TDC Mult ",510,0,510,10,0,10);


TH2F* h_grinch_le_elem = new TH2F("h_grinch_le_elem"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,2000,0,2000);
TH2F* h_grinch_tot_elem = new TH2F("h_grinch_tot_elem"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) ",510,0,510,101,0,100);

TH2F* h_grinch_le_elem_hodo = new TH2F("h_grinch_le_elem_hodo"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) - hodo tmean ",510,0,510,2600,0,2600);


TH2F* h_grinch_mult_elem = new TH2F("h_grinch_mult_elem"," ; GRINCH TDC elemID ; GRINCH TDC Multiplicty ",510,0,510,10,0,10);
TH1F* h_grinch_elem = new TH1F("h_grinch_elem",";GRINCH elem ID;",510,0,510);



TH1F* h_grinch_xhit =  new TH1F("h_grinch_xhit","bb.grinch_tdc.hit.xhit",59,-0.9145,0.9145);
TH1F* h_grinch_yhit = new TH1F("h_grinch_yhit","bb.grinch_tdc.hit.yhit",17,-0.124,0.124);


TH1F* h_grinch_le_all = new TH1F("h_grinch_le_all","; GRINCH LE ALL ;", 2500,0,2500);


TH1F* h_grinch_tot_all =new TH1F("h_grinch_tot_all","; GRINCH TOT ALL ;", 101,0,100);


TH1F* h_grinch_cluster_le_all = new TH1F("h_grinch_cluster_le_all","; GRINCH CLUSTER LE ALL ;", 2500,0,2500);
TH1F* h_grinch_cluster_tot_all = new TH1F("h_grinch_cluster_tot_all","; GRINCH CLUSTER TOT ALL ;", 101,0,100);

TH1F* h_grinch_cluster_tot[511];
TH1F* h_grinch_cluster_le[511];
TH2F* h_grinch_cluster_le_tot[511];


TH1F* h_kine_W =  new TH1F("h_kine_W", "; W ;", 500,-10,10);





//function declarations
void make_row_col_array(); 
void Fill_Search_Stacks();
Bool_t Is_NOT_In_Stack(stack <Int_t> checkstack, Int_t PMT);
void Print_Stack(stack <Int_t> inputstack);
stack <Int_t> Clear_Stack();
Int_t Sum_Adjacent_Hits(Int_t PMT, stack <Int_t>& neighbors_stack);
stack <Int_t> Add_to_Stack(stack <Int_t> inputstack, Int_t PMT);
void Calculate_PMT_Coord(double& horiz, double& vert, Int_t PMT); 
void Calculate_Cluster_Center_Row_Col(stack <Int_t> inputstack, double& horiz, double& vert);
void Calculate_Cluster_Center_ypos_xpos(stack <Int_t> inputstack, double& ypos, double& xpos);
void Find_Cluster_Extrema(stack <Int_t> inputstack, Int_t& horizspread, Int_t& vertspread);
void Fill_Cluster_Histos(stack <Int_t> inputstack, Int_t index);
void Fill_Cluster_Histos_struct(struct Cluster cluster_struct_array[], Int_t entries);
Int_t Stack_Size(stack <Int_t> inputstack);
void Calculate_Cluster_Time(stack <Int_t> inputstack, Double_t& mean, Double_t& variance);
void SelectionSort(Int_t a[], Int_t tracker[], Int_t n);
void StackToArray(stack <Int_t> s, Int_t a[]);



void cluster_finding_barebones(TString basename="",Int_t nrun=2043,TString configfilename="run_list_grinch.txt", Int_t entries = -1){ //MAIN
  if (basename=="") {
    cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
    cin >> basename;
  }
  TString fullname = Form("%s_%d",basename.Data(),nrun);
  // gStyle->SetPalette(kDarkBodyRadiator);
  //gStyle->SetPalette(1);
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
  TString inputroot;
  inputroot="rootfiles/"+fullname+".root";
  
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
  
  fchain ->SetBranchStatus("*",0);


  fchain ->SetBranchStatus("bb.ps.e",1);
  fchain->SetBranchAddress("bb.ps.e",&ps_e) ; // preshower energy

  fchain->SetBranchStatus("bb.sh.e",1);
  fchain->SetBranchAddress("bb.sh.e",&sh_e) ;  // shower energy

  fchain->SetBranchStatus("BB.gold.p",1);
  fchain->SetBranchAddress("BB.gold.p",&pmom) ; //

  fchain->SetBranchStatus("BB.gold.th",1) ; //
  fchain->SetBranchAddress("BB.gold.th",&xptar) ; //

  fchain->SetBranchStatus("BB.gold.ph",1) ; //
  fchain->SetBranchAddress("BB.gold.ph",&yptar) ; //

  fchain->SetBranchStatus("BB.gold.y",1) ; //
  fchain->SetBranchAddress("BB.gold.y",&ytar) ; //

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
  fchain->SetBranchAddress("Ndata.bb.grinch_tdc.hit.pmtnum",&GrinchNum) ;//number of PMTs with a signal in an event (I think)  
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
  

  fchain->SetBranchStatus("g.trigbits",1); 
  fchain->SetBranchAddress("g.trigbits",&gTrigBits); // The type of trigger the event was. (1 is the bbcal trig, 16 is the grinch LED) (May not be working)

  fchain->SetBranchStatus("fEvtHdr.fTrigBits",1);
  fchain->SetBranchAddress("fEvtHdr.fTrigBits",&fTrigBits);

  //fchain->SetBranchStatus("fEvtHdr.fRun",1);
  //fchain->SetBranchAddress("fEvtHdr.fRun",&fRun);

  
  
  
 
  //Making histos for the individual channels 
  for (Int_t ig=0;ig<511;ig++) {    
    h_grinch_cluster_tot[ig] = new TH1F(Form("h_grinch_cluster_tot_%d",ig),Form(" ; GRINCH TDC ToT CLUSTER PMT %d ; ",ig),50,0,50);
    h_grinch_cluster_le[ig] = new TH1F(Form("h_grinch_cluster_le_%d",ig),Form(" ; GRINCH TDC LE CLUSTER PMT %d ; ",ig),2000,0,2000); 
    h_grinch_cluster_le_tot[ig] = new TH2F(Form("h_grinch_cluster_le_tot_%d",ig),Form("  GRINCH TDC LE vs ToT CLUSTER PMT %d; cluster LE ; cluster ToT; ",ig),40,180,220,50,0,50);
  }
 
    
 
  //// functions that initialize a few things 
  Fill_Search_Stacks(); //sets up cluster finding stacks
  make_row_col_array(); //sets up an array to quickly get the row/col from PMT number

  
 

 

  TH1F* h_grinch_le[511];
  TH1F* h_grinch_tot[511];
  TH1F* h_grinch_hodo_le[511];
  TH2F* h_grinch_le_tot[511];
  for (Int_t ig=0;ig<511;ig++) {
    h_grinch_le[ig] = new TH1F(Form("h_grinch_le_%d",ig),Form(" ; GRINCH TDC LE (ns) PMT %d  ; ",ig),2600,0,2600);
    h_grinch_tot[ig] = new TH1F(Form("h_grinch_tot_%d",ig),Form(" ; GRINCH TDC ToT (ns) PMT %d  ; ",ig),50,0,50);
    h_grinch_le_tot[ig] = new TH2F(Form("h_grinch_le_tot_%d",ig),Form("  GRINCH TDC LE vs ToT PMT %d; LE ; ToT; ",ig),200,100,300,50,0,50);
    h_grinch_hodo_le[ig] = new TH1F(Form("h_grinch_hodo_le_%d",ig),Form(" ; GRINCH TDC LE (ns) -hodo_tmean  PMT %d  ; ",ig),2600,0,2600);
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
  double x=0;
  double y=0;  
  for (Int_t pmtnum = 0; pmtnum<510;pmtnum++)
    {
     Calculate_PMT_Coord(x, y, pmtnum);
    }
 
 

  //// Loop over the number of entries /////////
  for (int entry = 0; entry < max ; entry++) {//nentries
    fchain->GetEntry(entry);
    if (entry%5000==0) cout << " Entry = " << entry << endl;

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
    

   

    Bool_t electron_flag = kFALSE; //flag for when an event passes the shower cut
    Bool_t g_cluster_flag =kFALSE;//flag for when an event has a grinch cluster
    Bool_t good_event_flag = kFALSE;

    h_kine_W ->Fill(kine_W);
    
    if (tr_n ==1 && ps_e > 0.22 && gem_track_nhits > 3 && gTrigBits == 4 && abs(tr_vz[0])<0.05 && abs(tr_tg_th[0])<0.1 && abs(tr_tg_ph[0])<0.03)// gen test cut
      {
	////
	electron_flag = kTRUE; //throw up a flag that this is a good electron hit
	////	
	h_good_event_track_xy->Fill(tr_y[0], tr_x[0]);  
	
      } // END energy and track cut

    for (int n = 0; n<510;n++) //initialize some more arrays
      {
	hit_flag_array[n]=0;
	mult_array[n]=0;
	root_index_array[n]=0;
      }

   
    Double_t grinch_minus_hodo_time;
    Double_t goodhit=0;
    hit_sum = hit_sum + GrinchNum; //adding up all of the hits from the multi-hit TDC to normalize later
    h_grinch_hitslistlength ->Fill(GrinchNum);
    //// Loop over the PMTs that had a signal for this event 
    for (Int_t ig=0 ; ig < GrinchNum ; ig++) { //since PMTs with no signal are not written to the root file, this "ig" is NOT the PMT number
      Int_t gindex=tdcGID[ig]; // we need to specifically ask which PMT this signal is for
      h_grinch_elem->Fill(gindex);
      if (gindex <511 && gindex >-1) {
	//	h_grinch_le_elem->Fill(tdcGID[ig],tdcGLe[ig] - hodo_tmean);
	h_grinch_le[gindex]->Fill(tdcGLe[ig] - hodo_tmean); // 
	h_grinch_hodo_le[gindex]->Fill(tdcGLe[ig] -hodo_tmean); //
	h_grinch_le_tot[gindex] ->Fill(tdcGLe[ig], tdcGTot[ig]);
	h_grinch_le_all ->Fill(tdcGLe[ig]-hodo_tmean);
	h_grinch_tot_elem ->Fill(tdcGID[ig],tdcGTot[ig]);
	h_grinch_tot_all ->Fill(tdcGTot[ig]);
	h_grinch_tot[gindex] ->Fill(tdcGTot[ig]);
	
	grinch_minus_hodo_time = tdcGLe[ig] - hodo_tmean;
	h_grinch_le_elem_hodo ->Fill(tdcGID[ig], grinch_minus_hodo_time);

	if (abs(tdcGLe[ig]-200)<20 && tdcGTot[ig]>0 && electron_flag){ // if the signal in the PMT is within the good timing cut for the grinch and it was an electron event // was 900 <25 for GMn
	  good_event_flag = kTRUE;
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
	  
	  ////////////////////////////
	}     
      }         
     
    }//end loop over TDC PMTs
  

    clusterstack = Clear_Stack(); //clearing the stack we are using for cluster finding for the new event. Important to do. 
    
   
    if ( good_event_flag)
      {

      for(int m = 0; m<510;m++) //fill with 0 
	{
	  sum_array[m]=0;
          check_status[m]=0;
	}

      //// Start of Cluster Finding
      
      cluster_stack_cnt=0;
      Int_t top;
      Int_t check_top = 0;
      stack <Int_t> neighbors_stack_temp;
      Int_t neighbors_top = 0;

      for (int pmtcntr = 0; pmtcntr<=509 ; pmtcntr++)
	{
	  neighbors_stack[pmtcntr] = Clear_Stack();
	  sum_array[pmtcntr] =  Sum_Adjacent_Hits(pmtcntr, neighbors_stack[pmtcntr]);
	  check_status[pmtcntr] = sum_array[pmtcntr];
	}

      for (int i = 0;i < 16 ; i++)
	{
	  cluster_stack[i]=Clear_Stack();
	}
      check_stack = Clear_Stack();

      clusterflag=kFALSE;
      cluster_stack_cnt = 0;
      for (int pmtcntr = 0; pmtcntr<=509 ; pmtcntr++)
	{
	  if(sum_array[pmtcntr]>=1 && check_status[pmtcntr]!=-1 ) // if a PMT has 1 or more neighboring PMTs with hits and it hasn't already been counted
	    {
	      clusterflag = kTRUE;
	      check_stack.push(pmtcntr); // put that PMT in the stack to go check it. 

	      while(!check_stack.empty())
		{
		  check_top = check_stack.top(); // PMT on the top of the stack has neighbors and hasn't been counted yet
		  check_stack.pop();
		  cluster_stack[cluster_stack_cnt].push(check_top); // put it in the current cluster stack
		  neighbors_stack_temp = neighbors_stack[check_top]; // get the neighbors stack for the PMT we are looking at. 
		  check_status[check_top] = -1; // mark that this PMT has been checked 
		
		  while(!neighbors_stack_temp.empty()) //looping through the neighbors of the PMT we are looking at
		    {
		      neighbors_top = neighbors_stack_temp.top();
		      neighbors_stack_temp.pop();
		      if ( check_status[neighbors_top] != -1 && neighbors_top != check_top) // if we haven't checked the neighbor yet 
			{
			  check_stack.push(neighbors_top); // put it in the stack to go check it
			}
		      check_status[neighbors_top] = -1; // mark that we put it in the stack to check (possibly redundent to have this here)
		    }	
		}
	      cluster_stack_cnt ++; // if we are here, that means that we emptied out the check_stack and spirdered through all the PMTs connecting to the first PMT we put in the stack. 
	      ////The next PMT that we come across that hasn't been checked is thus physcially seperated from the last cluster. 
	    }
	}//// end loop over all PMTs

   
      g_cluster_flag = kFALSE;

      if(cluster_stack_cnt > 0)
	{
	  h_grinch_cluster_cnt ->Fill(cluster_stack_cnt);
	}


      for (int i = 0; i < cluster_stack_cnt; i++)
	{
	  //// Calling functions to analyze the clusters we just found
	  ClusterSize[i]=Stack_Size(cluster_stack[i]);
	  Calculate_Cluster_Center_Row_Col(cluster_stack[i],cluster_horiz[i],cluster_vert[i]);
	  Calculate_Cluster_Center_ypos_xpos(cluster_stack[i],cluster_ypos[i],cluster_xpos[i]);
	  Find_Cluster_Extrema(cluster_stack[i],horizspread[i],vertspread[i]); // function measures the height and width of the cluster
	  Calculate_Cluster_Time(cluster_stack[i],cluster_mean_time[i],cluster_variance[i]); // calculates the mean time and variance 

	  if(ClusterSize[i] > 2) // clusters of size 2 seem to mostly be noise, so I am filling a few histograms without them 
	    {
	      g_cluster_flag = kTRUE; 
	      clustercnt++;
	      Fill_Cluster_Histos(cluster_stack[i], i); //function that fills some histograms 
	    }
	}

      if(cluster_stack_cnt > 0)
	{ //// Ordering the clusters from largest to smallest 
	  //cout<<"SelectionSort"<<endl;
	  Int_t tracker[32]= {-1}; // keeps track of how the array is re-ordered 
	  for (Int_t i = 0 ; i<cluster_stack_cnt ; i++)
	    {
	      ClusterSize_ordered[i] = ClusterSize[i]; //making a copy of the array so we can leave the unsorted version unchanged 
	    }

	  SelectionSort(ClusterSize_ordered,tracker,cluster_stack_cnt); // Sorting from Largest cluster to smallest cluster. Arrays act like reference params

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

	  struct Cluster cluster_struct_array[32]; // Array of "Cluster" structs that hold all the info 

	  for (Int_t i = 0 ; i < cluster_stack_cnt ; i++)
	    {
	      //// Filling the struct with the info that has only one value per cluster 
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

	      for ( Int_t j = 0 ; j < ClusterSize_ordered[i] ; j++) //looping over all the PMTs in the cluster
		{
		  //// Filling the struct with the info for each pmt in the cluster 
		  Int_t rootID = root_index_array[cluster_array[j]];
		  cluster_struct_array[i].PMT[j] = cluster_array[j];
		  cluster_struct_array[i].rootID[j] = rootID;
		  cluster_struct_array[i].LE[j] =  tdcGLe[rootID];
		  cluster_struct_array[i].ToT[j] = tdcGTot[rootID];
		}
	    }


	  Fill_Cluster_Histos_struct(cluster_struct_array, cluster_stack_cnt); //function to fill some histograms using the data in the struct. 

	  //// Everything below here is trying to determine what the "best" cluster might be.

	  ////****************************************************************************


	  //// print statements for debugging 
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

	  projx = tr_x[0] + tan(tr_th[0])*grinch_distance;
	  projy = tr_y[0] + tan(tr_ph[0])*grinch_distance;

	  h_grinch_projx -> Fill(projx);
	  h_grinch_projy -> Fill(projy);

	  Double_t min_dx = 500;
	  Int_t min_tracker = 0;
	  Int_t cluster_of_three_or_more_pmts_cnt = 0;
     
	  for (Int_t i = 0 ; i < cluster_stack_cnt ; i++)
	    {
	      grinch_dx[i] = cluster_xpos_ordered[i] - projx;
	      grinch_dy[i] = cluster_ypos_ordered[i] - projy;

	      if(ClusterSize_ordered[i] > 2)
		{
		  h_grinch_cluster_dx ->Fill(grinch_dx[i]);
		  h_grinch_cluster_dy ->Fill(grinch_dy[i]);

		  h_grinch_cluster_projx_xpos ->Fill(cluster_xpos_ordered[i], projx);
		  h_grinch_cluster_projy_ypos ->Fill(cluster_ypos_ordered[i], projy);
	     
		  cluster_of_three_or_more_pmts_cnt ++;
	     
		  if(abs(grinch_dx[i]) < abs(min_dx))
		    {
		      min_dx = abs(grinch_dx[i]);
		      min_tracker = i;
		    }
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
	      h_grinch_cluster_projy_ypos_best ->Fill(cluster_ypos_ordered[min_tracker], projy);

	      if(cluster_of_three_or_more_pmts_cnt > 1)
		{
		  h_grinch_cluster_projx_xpos_multi_best ->Fill(cluster_xpos_ordered[min_tracker], projx);
		  h_grinch_cluster_projy_ypos_multi_best ->Fill(cluster_ypos_ordered[min_tracker], projy);

		}
	    }
	}
 

      // if(cluster_stack_cnt > 0){
      //   cout<<"-------------------------"<<endl;
      //   cout <<"Event #"<< entry<< endl;
      //   cout<<"clusters found: "<<cluster_stack_cnt <<endl;
      //   cout<<"-------------------------"<<endl;
      //   cout<< "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " <<endl;
      // }

     
      temp_cnt ++;
      
      
    }// end "if good event"
       

  }//end event loop

    
 
  cout<<"cluster analysis on "<<max<<" events"<<endl;
  cout<<clustercnt<<" clusters found"<<endl;
 
  cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
  

}// end main



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
 
  for (Int_t i = 0; i < entries ; i++)
    {
      if (cluster_struct_array[i].size > 2){
	h_grinch_cluster_tmean -> Fill(cluster_struct_array[i].mean_time);
	h_grinch_cluster_variance ->Fill(cluster_struct_array[i].sqrt_variance);
	h_grinch_cluster_size -> Fill(cluster_struct_array[i].size);
	h_grinch_cluster_center_xpos ->Fill(cluster_struct_array[i].xpos);
	h_grinch_cluster_center_ypos -> Fill(cluster_struct_array[i].ypos);
	h_grinch_cluster_center ->Fill(cluster_struct_array[i].col,cluster_struct_array[i].col);
	h_grinch_cluster_spread->Fill(cluster_struct_array[i].horizspread,cluster_struct_array[i].vertspread);

	for (Int_t j = 0 ; j < cluster_struct_array[i].size ; j++)
	  {
	    pmt = cluster_struct_array[i].PMT[j];
	    LE = cluster_struct_array[i].LE[j];
	    ToT = cluster_struct_array[i].ToT[j]; 

	    h_grinch_cluster_elem -> Fill(pmt);
	    h_grinch_cluster_le_elem ->Fill(pmt, LE);
	    h_grinch_cluster_tot_elem ->Fill(pmt, ToT);
	    h_grinch_cluster_le[pmt]->Fill(LE);
	    h_grinch_cluster_tot[pmt]->Fill(ToT);
	    h_grinch_cluster_le_all -> Fill(LE);
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
  Double_t grinch_minus_hodo_time;

  while (!inputstack.empty())
    {
      PMT = inputstack.top();
      inputstack.pop();
      root_id = root_index_array[PMT];//get the index that the branches need in order to look at the specific PMT for this event
      h_grinch_cluster_le_elem -> Fill(PMT,tdcGLe[root_id] - hodo_tmean);
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
{ // Converts a stack into an array. 
  stack <Int_t> scopy= s;
  Int_t i = 0;
  while (!scopy.empty())
    {
      a[i] = scopy.top();
      scopy.pop();
      i++;
    }
}// end StackToArray

void make_row_col_array()// making it easy to convert from PMT num to row, col
//// this is also probably something that is in the database or should go in the database
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
//// Example: Adding 20 to stack |-9,-8,-1,1,8,9| becomes |11,12,19,21,28,29|
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


