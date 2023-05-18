#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
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
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>

#include "TGLHistPainter.h"
#include "TGLBoxPainter.h"
#include <stdlib.h>

#include <string>
#include "TH1.h"
#include "TRandom.h"
#include <vector>
#include <TGClient.h>
#include <TF1.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <stack>



//using namespace std;




TString inputroot;
TFile *file;
TTree *eventtree;
int N_evts;


#define N_ROW 60

//TTree *eventtree;

//int Event=-1;
int Event;
Int_t event_cnt = 0;


//UTEventMin *event;
const Int_t maxHitsPerEvent = 100;
double* ped=NULL;
int* Hit_Events=NULL;
TCanvas* c;
TCanvas* c_heat;
TCanvas* c2;
TCanvas* allC;
TCanvas* cADC;
TCanvas* allC_heat;

//int N_evts;

int eventID;

const int N_PMT=512;

const Int_t FORWARD = 1;
const Int_t BACKWARD = 2;

Int_t hit_flag_array[510];
const Int_t WAS_HIT = 1;
const Int_t NOT_HIT = 0;

//Array to count hits on each tube for entire run.
double TubeHits[N_PMT];
double multiTubeHits[N_PMT];
//Array to count hits on each tube for individual events.
double EvtTubeHits[N_PMT] = {};


// SBS-Offline variables
Int_t ntdc=510;
Double_t grinch_tdc_pmt[1000];
Double_t grinch_tdc_le[1000];
Double_t grinch_tdc_te[1000];
Double_t grinch_tdc_tot[1000];
Double_t grinch_tdc_mult[1000];

Double_t grinch_pmt_row[1000];
Double_t grinch_pmt_col[1000];

Double_t ps_e; 
Double_t sh_e; 
Double_t ntrack; 
Double_t pmom; 

Int_t GrinchADCNum;
Double_t adcGID[1000]; // raw adc
Double_t adcGAtime[1000]; // raw adc
Double_t adcGAmp[1000]; // raw adc
Double_t grinch_adc[1000];//maria
Double_t adcGMult[1000]; // raw adc 

Double_t tr_x[10];
Double_t tr_y[10];
Double_t tr_n;
Double_t tr_vz[100];
Double_t tr_tg_th[100];
Double_t tr_tg_ph[100];
Double_t gem_track_nhits;
Double_t tr_p[100];
Double_t hcal_e;


// Histograms
TH1F *h_ntdc = new TH1F("h_ntdc"," ; Nhit ; Counts",512,0,511); // 
TH2F *h_tdc_vs_pmt = new TH2F("h_tdc_vs_pmt"," ; PMT # ; TDC LE", 512,0,511,150,1,150);
TH2F *h_adc_vs_chan = new TH2F("h_adc_vs_chan"," ; chan # ; ADC counts", 64,0,63,200,1,200);
TH2F *GRINCH_mult_vs_pmt = new TH2F("GRINCH_mult_vs_pmt"," ; PMT # ; PMT hit multiplicity", 510,0,510,5,1,6);
TH1F *h_GRINCH_mult = new TH1F("h_GRINCH_mult","; multiplicity ; counts",510,0,510);

TH2F* h_grinch_cluster_le_elem = new TH2F("h_grinch_cluster_le_elem"," ; GRINCH TDC elemID ; GRINCH TDC LE (ns) ",510,0,510,100,850,950);
TH2F* h_grinch_cluster_te_elem = new TH2F("h_grinch_cluster_te_elem"," ; GRINCH TDC elemID ; GRINCH TDC TE (ns)",510,0,510,100,850,950);
TH2F* h_grinch_cluster_tot_elem = new TH2F("h_grinch_cluster_tot_elem"," ; GRINCH TDC elemID ; GRINCH TDC TOT (ns) ",510,0,510,110,-10,100);
TH2F* h_grinch_cluster_mult_elem = new TH2F("h_grinch_cluster_mult_elem"," ; GRINCH TDC elemID ; GRINCH TDC Mult ",510,0,510,15,0,15);
TH2F* h_grinch_cluster_adcMult_elem = new TH2F("h_grinch_cluster_adcMult_elem"," ; GRINCH ADC elem id ; GRINCH ADC Mult. ",63,0,63,10,0,10);
TH2F* h_grinch_cluster_atime_elem = new TH2F("h_grinch_cluster_atime_elem"," ; GRINCH ADC elemID ; GRINCH ADC time (ns) ",63,0,63,500,0,2000);
TH2F* h_grinch_cluster_amp_elem = new TH2F("h_grinch_cluster_amp_elem"," ; GRINCH ADC elem id ; GRINCH ADC amp (mV) ",63,0,63,200,0,200);    
TH2F* h_grinch_cluster_int_elem = new TH2F("h_grinch_cluster_int_elem","; GRINCH ADC elemID ; GRINCH ADC channel ",64,0,63,200,0,200);

TH2F* h_grinch_cluster_adc_tdc =new TH2F("h_grinch_cluster_adc_tdc","; GRINCH ADC AMP ; GRINCH TDC ToT ",64,0,63,110,-10,100);




Int_t i,j;
double cut = 2000;
double MaxHits = 100;

Int_t next_evt=0;
Int_t usr_evt=0;

Int_t tdc_event[N_PMT][10];
Int_t tdcTimeArray[N_PMT];
//Double_t tdcTime[N_PMT];
double tdcTime=0;
Int_t sum_array[510];
Int_t root_index_array[510];
Int_t adc_root_index_array[510];

Int_t row_array[510];
Int_t col_array[510];
Int_t tdc_to_adc_map_array[510]; // Needs to be updated when ADCs are moved to different Nino cards. 



TPaveLabel** Signal_Label=new TPaveLabel*[N_PMT];
TEllipse** PMT=new TEllipse*[N_PMT];
TPaveLabel** PMT_Label=new TPaveLabel*[N_PMT];
Int_t eventEntries;

Color_t hit_color[10]={kRed, kPink+1, kMagenta, kViolet+7, kBlue, kAzure-4, kCyan, kTeal+3, kGreen, kYellow};
//Color_t hit_color_excess = kBlack;

TControlBar* bar = new TControlBar("","");



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
stack <Int_t> cluster_members_stack;

//function headers
void make_row_col_array();
void Fill_Search_Stacks();
void Print_Stack(stack <Int_t> inputstack);
void set_color_env();
void draw_event_display();
void fill_event(Int_t Direction);
//Int_t Sum_Adjacent_Hits(Int_t row, Int_t col, Int_t PMT);
Int_t Sum_Adjacent_Hits(Int_t PMT,stack <Int_t>& neighbors_stack);
Bool_t Is_NOT_In_Stack(stack <Int_t> checkstack, Int_t PMT);
void show(int option=0);
stack <Int_t> Clear_Stack();
stack <Int_t> Add_to_Stack(stack <Int_t> inputstack, Int_t PMT);
void Calculate_PMT_Coord(double& horiz, double& vert, Int_t PMT);
void Calculate_Cluster_Center(stack <Int_t> inputstack, double& horiz, double& vert);
void Test_Reference(double& a, double& b, Int_t test);
void Find_Cluster_Extrema(stack <Int_t> inputstack, Int_t& horizspread, Int_t& vertspread);
void Make_Histos_for_the_Cluster(stack <Int_t> inputstack);
void make_tdc_to_adc_map_array();
Int_t Stack_Size(stack <Int_t> inputstack);




void GRINCH_eventDisplay_m_simple()//main
{

  //inputroot="../rootfiles/grinch_13369_fev_0_nev_-1_fseg_0_lseg_0.root";
  // inputroot = "../rootfiles/grinch_13620_fev_0_nev_-1_fseg_0_lseg_0.root";

  inputroot = "/volatile/halla/sbs/seeds/rootfiles/hcal_gmn_fullreplay_13559_stream0_seg1_1_1.root";


  file = new TFile(inputroot); 
  eventtree = (TTree*) file->Get("T");
  eventEntries = eventtree -> GetEntries();
  
  cout<< "eventEntries = "<< eventEntries<<endl;


 
  // Fill the arrays from ROOT branches (variables to be used anywhere)
  ///////////////////////////////////////////////////////////////////////////////
  eventtree->SetBranchAddress("Ndata.bb.grinch_tdc.tdc",&ntdc) ; // Number of TDC (LE) in event
  eventtree->SetBranchAddress("bb.grinch_tdc.tdcelemID",&grinch_tdc_pmt) ; // PMT number
  eventtree->SetBranchAddress("bb.grinch_tdc.tdc",&grinch_tdc_le) ; // TDC LE
  eventtree->SetBranchAddress("bb.grinch_tdc.tdc_te",&grinch_tdc_te) ; // TDC TE
  eventtree->SetBranchAddress("bb.grinch_tdc.tdc_tot",&grinch_tdc_tot);//ToT
  eventtree->SetBranchAddress("bb.grinch_tdc.tdc_mult",&grinch_tdc_mult) ; // TDC multiplicity
  eventtree->SetBranchAddress("bb.grinch_tdc.tdcrow",&grinch_pmt_row) ; // TDC row
  eventtree->SetBranchAddress("bb.grinch_tdc.tdccol",&grinch_pmt_col) ; // TDC column

  eventtree->SetBranchAddress("bb.ps.e",&ps_e) ; 
  eventtree->SetBranchAddress("bb.sh.e",&sh_e) ;
  eventtree->SetBranchAddress("bb.tr.n",&ntrack) ;
  eventtree->SetBranchAddress("BB.gold.p",&pmom) ;

  eventtree->SetBranchAddress("Ndata.bb.grinch_adc.adcelemID",&GrinchADCNum) ;
  eventtree->SetBranchAddress("bb.grinch_adc.adcelemID",&adcGID) ;
  eventtree->SetBranchAddress("bb.grinch_adc.a_time",&adcGAtime) ;
  eventtree->SetBranchAddress("bb.grinch_adc.a_amp_p",&adcGAmp) ; //amp_p
  eventtree->SetBranchAddress("bb.grinch_adc.a_mult",&adcGMult) ;//maria
  eventtree->SetBranchAddress("bb.grinch_adc.a_p",&grinch_adc) ; // ADC counts //maria a_p

  eventtree ->SetBranchAddress("bb.tr.x",&tr_x);
  eventtree ->SetBranchAddress("bb.tr.y",&tr_y);
  eventtree  ->SetBranchAddress("bb.tr.n", &tr_n);
  eventtree  ->SetBranchAddress("bb.tr.vz", &tr_vz);
  eventtree  ->SetBranchAddress("bb.tr.tg_th",&tr_tg_th);
  eventtree  ->SetBranchAddress("bb.tr.tg_ph",&tr_tg_ph);
  eventtree  ->SetBranchAddress("bb.gem.track.nhits",&gem_track_nhits);
  eventtree  ->SetBranchAddress("bb.tr.p",&tr_p);
  eventtree  ->SetBranchAddress("sbs.hcal.e", &hcal_e);
  ////////////////////////////////////////////////////////////////////////////////// 
 

  make_row_col_array();
  make_tdc_to_adc_map_array(); 

  Fill_Search_Stacks();
  
  

  //Create GUI
  bar = new TControlBar("vertical", "GRINCH Event Display");
  bar->AddButton("GRINCH Event Display", "");
  bar->AddButton("", "");
  bar->AddButton("First Event", "show(-2)");
  bar->AddButton("Next Event", "show(1)");
  bar->AddButton("Prev. Event", "show(6)");
  bar->AddButton("Last Event", "show(-3)");
  bar->AddButton("", "");
  bar->AddButton("Choose Event #", "show(0)");
  bar->AddButton("", "");
  //bar->AddButton("Heat Map", "show(4)");
  //bar->AddButton("PMT Counter", "show(5)");
  bar->AddButton("", "");
  bar->AddButton("Exit",".q");
  bar->Show();

}// end main



void fill_event(Int_t Direction)
{  
 
  if ( allC != NULL)
    {
      allC->Clear();
    }
  else
    {
      allC=new TCanvas("allC","Multiplicity vs. PMT",1000,800);
    }

  //allC -> SetCanvasSize(700, 400);
  //allC -> SetWindowSize(700, 450);

  allC -> Divide(2,2);

  if( cADC != NULL)
    {
      cADC ->Clear();
    }
  else
    {
      cADC = new TCanvas("cADC","cADC",1000,800);
    }

  cADC ->Divide(2,2);

  if ( c != NULL)
    {
      c->Clear();
    }
  else
    {
      c=new TCanvas("c");
    }

   c -> Range(0,0,1.5, 1.8);// This changes the range of the canvas
   c -> SetCanvasSize(875,875);
  // c -> SetWindowSize(400, 2500);

  if( c2 )
    {
      c->cd();
    }


  // draw_event_display();

  Int_t hit_row_array[510];
  Int_t hit_col_array[510];
   
   
  const Int_t nChanADC = 64;
  const Int_t nChanVETROC = 510;

  //const Int_t TDC_array = 1000;
  const Int_t TDC_array = 10;


  Bool_t eventflag = kFALSE; 
  Bool_t clusterflag= kFALSE;

  stack <Int_t> clusterstack;
  stack <Int_t> emptystack;

  Bool_t sh_energy_flag = kFALSE;
 
  stack <Int_t> cluster_stack;
  stack <Int_t> neighbors_stack;

  while (!clusterflag){
    eventflag=kFALSE;

    while (!eventflag){
      //cout<<"event  = "<<event_cnt<<endl;

      eventtree->GetEntry(event_cnt);

      Double_t tot_e = ps_e+sh_e;
      Double_t rat = tot_e/pmom;

      sh_energy_flag = kFALSE;
      // if(abs(tot_e - 3)< 0.5 && ps_e > 0.2 && abs(rat - 1.2) < 0.2 )
      if(tr_n ==1 && abs(tr_vz[0])<0.08 && gem_track_nhits > 3 && tr_p[0] >1.4 && tr_p[0]<2.0 && hcal_e > 0.025 && ps_e >0.22)//SBS9
	{
	  sh_energy_flag = kTRUE;
	
    
      
	  for (int n=0; n<64; n++)
	    {
	      adc_root_index_array[n]=0;
	    }


	  for(int r=0; r<N_PMT; r++)  //Reset individual events tube hits counter for each event.
	    {
	      EvtTubeHits[r] = 0;
	      tdcTimeArray[r] = 0;
	      hit_flag_array[r]=0;
	      hit_row_array[r]= -1;
	      hit_col_array[r] = -1;
	      root_index_array[r]=0;
	    }
   
	  GRINCH_mult_vs_pmt -> Reset();
	  h_grinch_cluster_le_elem -> Reset();
	  h_grinch_cluster_te_elem -> Reset();
	  h_grinch_cluster_tot_elem -> Reset();
	  h_grinch_cluster_mult_elem -> Reset();
	  h_grinch_cluster_adcMult_elem -> Reset();
	  h_grinch_cluster_atime_elem -> Reset();
	  h_grinch_cluster_amp_elem -> Reset();
	  h_grinch_cluster_int_elem -> Reset();
	  h_grinch_cluster_adc_tdc ->Reset();
      
      
	  // I want to add ADC info to look at the clusters. 
	  for(Int_t ig=0; ig<GrinchADCNum; ig++)
	    {
	      Int_t gindex_adc = adcGID[ig];
	      adc_root_index_array[gindex_adc] = ig;
	    }
      
 
	  for(int j=0; j<nChanVETROC; j++)
	    {       
	      Double_t temp_mult=0;
	      Int_t temp_hit_flag=0;//was it hit or not
	      Int_t temp_row = -1;
	      Int_t temp_col = -1;
	  
	
	      for (Int_t ig=0;ig<ntdc;ig++) {
		if (grinch_tdc_pmt[ig] == j && abs(grinch_tdc_le[ig]-900) < 20)
		  {
		    temp_mult=grinch_tdc_mult[ig];
		    temp_hit_flag = 1;
		    temp_row = grinch_pmt_row[ig];
		    temp_col = grinch_pmt_col[ig];
		    root_index_array[j] = ig; // Map the "gindex" (which is the PMT number) to the "ig index" which is what the branch needs.
		    // I want to be able to go back after finding clusters to look at ADC and TDC and whatnot on those PMTs.
		    // cout<<"temp_row= "<<temp_row <<", temp_col= "<<temp_col<<" for PMT "<< j <<endl;

		  }
	      }
	    
	      h_GRINCH_mult->SetBinContent(j,temp_mult);
	      GRINCH_mult_vs_pmt->Fill( j , temp_mult );
	
	      tdcTimeArray[j] = temp_mult;
	      hit_flag_array[j] = temp_hit_flag;
	      hit_row_array[j] = temp_row;
	      hit_col_array[j] = temp_col;

	      if(tdcTimeArray[j] > 0)
		{
		  EvtTubeHits[j] = EvtTubeHits[j] + 1;
		  if(tdcTimeArray[j] <=10)
		    {
		      eventflag = kTRUE;
		    }
		}	  	
	    }
	}
      if(!eventflag){ // check to see if we found a hit in the detector for this event
	if(Direction == FORWARD && event_cnt < (eventEntries -1)){
	  event_cnt ++;
	} 
	else if(Direction == BACKWARD && event_cnt > 0){
	  event_cnt --;
	}
	else {
	  cout<<"You have asked to go out bounds. Breaking loop."<<endl;
	  eventflag = kTRUE;//break out of the while loop. 
	}
	//cout<<"No signals within cuts found for this event. Going to next event."<<endl;
      }
   
    }//end of event while loop
    
    // The event we are on now has a hit in at least 1 PMT 
    

    //cout<<"EVENT # = "<<event_cnt<<endl;
   
    for(int m = 0; m<510;m++) //fill with 0 
      {
	sum_array[m]=0;
      }
    Int_t row=-1;
    Int_t col=-1;
    clusterstack = Clear_Stack(); // clear the stack 
    
    //clusterflag = kFALSE;
    //// Look for clusters
    Int_t top;
    
    cluster_stack = Clear_Stack();
    for (int pmtcntr = 9; pmtcntr<=500 ; pmtcntr++)
      {
	row = hit_row_array[pmtcntr];
	col = hit_col_array[pmtcntr];
	neighbors_stack = Clear_Stack();
	//sum_array[pmtcntr] =  Sum_Adjacent_Hits(row,col,pmtcntr);
	sum_array[pmtcntr] =  Sum_Adjacent_Hits(pmtcntr, neighbors_stack);
	if (sum_array[pmtcntr] != 0)
	  {
	    //cout<<"sum_array [PMT] "<<pmtcntr<< " = "<<sum_array[pmtcntr]<<endl;
	  }
	if (sum_array[pmtcntr] >= 2)
	  {
	    //cout<< "Possible cluster around PMT "<<  pmtcntr << ": sum = " <<sum_array[pmtcntr]<<endl;
	    clusterflag = kTRUE;
	    clusterstack.push(pmtcntr); // put that PMT in the stack for the cluster size to be counted 
	    //Print_Stack(neighbors_stack);
	    //cout<<"neighbors stack test"<<endl;
	    while(!neighbors_stack.empty())
	      {
		top = neighbors_stack.top();
		if( Is_NOT_In_Stack(cluster_stack,top))
		  {
		    cluster_stack.push(top);
		  }
		neighbors_stack.pop();		
	      }
	  }
      }
     
    
    
    if(!clusterflag){
      if(Direction == FORWARD && event_cnt < (eventEntries -1)){
	event_cnt ++;
      } 
      else if(Direction == BACKWARD && event_cnt > 0){
	event_cnt --;
      }
      else {
	cout<<"You have asked to go out bounds. Breaking loop."<<endl;
	clusterflag = kTRUE;//break out of the while loop. 
      }
      // cout<<"No clusters within cuts found for this event. Going to next event."<<endl;
    }
  
  }//end of cluster while loop

 

  Print_Stack(cluster_stack); 
  cout<<"cluster_stack"<<endl;
  
  Int_t cluster_size =  Stack_Size(cluster_stack);

  cout<<"stack size = "<< cluster_size <<endl;
 


 
  //want to send cluster_members_stack to a function that makes histos of ADC/TDC for each PMT in the cluster. 
  Make_Histos_for_the_Cluster(cluster_stack);

  if(sh_energy_flag)
    {
      cout << "event had good shower energy"<<endl;
    }
  else if (!sh_energy_flag)
    {
      cout << "event did not have good shower energy"<<endl;
    }


  stack <Int_t> copystack = cluster_stack;
  Int_t pmttest=0;
  double x=0;
  double y=0;
  Int_t test = 0;
  double clusterx = 0;
  double clustery = 0;
   
  
  
  Calculate_Cluster_Center(cluster_stack,clusterx,clustery);
  cout<<"cluster center is at "<<clusterx<<", "<<clustery<<endl;
  
  Int_t vertspread = 0;
  Int_t horizspread = 0; 

  Find_Cluster_Extrema(cluster_stack,horizspread,vertspread);
  cout<<"cluster is "<<horizspread<<" wide, "<< vertspread<<" tall"<<endl;

  cout<< "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " <<endl;
  
  draw_event_display();

  for (Int_t j = 0; j<510 ; j++)//filling display with appropriate colors
    {
      if(tdcTimeArray[j] > 0)	 
	{
	  Int_t d = j;
	  Int_t color_shift = tdcTimeArray[d]-1;


	  //      cout<<"PMT: "<<d<<" value: "<<tdcTime<<endl;
	     
	  Signal_Label[d]->SetFillColor(hit_color[color_shift]);

	  if(tdcTimeArray[j]<=10)
	    {
	      Signal_Label[d]->SetTextColor(kBlack);
	      PMT[d]->SetFillColor(hit_color[color_shift]);
	      PMT_Label[d]->SetFillColor(hit_color[color_shift]);
	      // if (j>9 && j<500 &&  sum_array[j]>=2)
	      // 	{
	      // 	  PMT_Label[d]->SetFillColor(kWhite);
	      // 	}

	    }else{
	    Signal_Label[d]->SetTextColor(kRed);
	    PMT[d]->SetFillColor(kBlack);
	    PMT_Label[d]->SetFillColor(kRed);
	  }

	  // Signal_Label[d]->SetTextSize(2);
	  // // cout<<"Multiplicity = "<<tdcTimeArray[d]<<" in PMT = "<<d<<endl;
	  // Signal_Label[d]->SetLabel(Form("%i",tdcTimeArray[d]));
	  // Signal_Label[d]->Draw("same");

	}
    }
  
  

  allC->cd(1);
  //style(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1,0);  
  h_grinch_cluster_le_elem ->Draw("COLZ"); 
  allC->cd(2);
  h_grinch_cluster_te_elem ->Draw("COLZ");
  allC->cd(3);
  h_grinch_cluster_tot_elem ->Draw("COLZ");
  allC->cd(4);
  h_grinch_cluster_mult_elem ->Draw("COLZ");

  c->SetTitle(Form("Event %d",event_cnt));
  allC->SetTitle(Form("Event %d",event_cnt));

  cADC ->cd(1);
  h_grinch_cluster_amp_elem -> Draw("COLZ");
  cADC ->cd(2);
  h_grinch_cluster_int_elem -> Draw("COLZ");
  cADC ->cd(3);
  //h_grinch_cluster_atime_elem ->Draw("COLZ");
  h_grinch_cluster_adc_tdc ->Draw("COLZ");
  cADC ->cd(4);
  h_grinch_cluster_adcMult_elem ->Draw("COLZ");

  c->Modified(); c->Update(); gSystem->ProcessEvents();
  allC->Modified(); allC->Update(); gSystem->ProcessEvents();
  cADC->Modified(); cADC->Update(); gSystem->ProcessEvents();
  //allC2->Modified(); allC2->Update(); gSystem->ProcessEvents();
  //allC3->Modified(); allC3->Update(); gSystem->ProcessEvents();

  //c->SetTitle(Form("Event %d",Event));
  //allC->SetTitle(Form("Event %d",Event));

} // END FILL EVENT





void make_row_col_array()// making it easy to convert from PMT num to row, col
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



void set_color_env(){
  //Create a nice color gradient.
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);   //Chooses custom gradient NCont.
}


void draw_event_display()
{
  //*************************************
  //DOUBLE LOOP FOR EVENT DISPLAY CREATION
  ///////////////////////////////////////

  if (c == NULL) c=new TCanvas("c");

  //load_data();

  double PMT_dia=2./(59.0); //the whole canvas has a max size of 1. We divide by the number of rows. 
  double x,y;//center of PMT. Draws PMT array.

  double labelsize=0.2;
  Int_t k=0;

  double y_off = -0.25;

  //cout<<"Creating Event Display..."<<endl;
  for ( i=0; i<N_ROW; i++ )
    {
      int N_COL = 8 + i%2;
      for ( j=0; j<N_COL; j++ )
	{
	  if ( i%2==0 )
	    {
	      x=(j+1.5)*PMT_dia +0.1;
	      y=N_ROW*PMT_dia-(i*sin(60*TMath::DegToRad())+0.5)*PMT_dia + y_off;
	    }
	  else
	    {
	      x=(j+1)*PMT_dia  +0.1;
	      y=N_ROW*PMT_dia-(i*sin(60*TMath::DegToRad())+0.5)*PMT_dia + y_off;
	    }
	  // cout<<"ROW: "<<i<<" PMT: "<<k<<" X: "<<x<<" Y: "<<y<<endl;
	  PMT_Label[k]=new TPaveLabel(x-PMT_dia*labelsize, y-PMT_dia*labelsize, x+PMT_dia*labelsize, y+PMT_dia*labelsize, Form("%d",k));
	  PMT[k]=new TEllipse(x, y, PMT_dia*0.5, PMT_dia*0.5);

	  if ( i==0 && j==0 )
	    PMT[k]->Draw();
	  else
	    PMT[k]->Draw("same");

	  PMT_Label[k]->SetFillColor(0);
	  //PMT_Label[k]->SetBorderSize(1);
	  //PMT_Label[k]->SetTextSize(0.5);
	  PMT_Label[k]->SetTextSize(1);
	  PMT_Label[k]->Draw("same");

	  Signal_Label[k]=new TPaveLabel(x-1.5*PMT_dia*labelsize,y+PMT_dia*labelsize,x+1.5*PMT_dia*labelsize,y+2*PMT_dia*labelsize,"");
	  Signal_Label[k]->SetTextColor(2);
	  Signal_Label[k]->SetFillColor(0);
	  Signal_Label[k]->SetLabel("");
	  k++;
	}
    }//End building event display geometry.
  //cout<<"...Event Display Created"<<endl;
  //*************************************
} // Event Display function



Int_t Sum_Adjacent_Hits(Int_t PMT, stack <Int_t>& neighbors_stack)
{
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

  stack <Int_t> neighborstack = Add_to_Stack(searchstack,PMT); //"PMT" to every element of the stack 
  //Print_Stack(neighborstack);
  //cout<<"neighborstack"<<endl;
  
  neighbors_stack.push(PMT);
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

}// end Sum_Adjacent_Hits


/*
Int_t Sum_Adjacent_Hits(Int_t row, Int_t col, Int_t PMT)
{   
  if (row == -1 || col == -1 || PMT < 9 || PMT >500)
    {
      // don't want to sum these. return zero
      
      return 0;
    }

  Int_t hit_sum = -1;//
 
  //cout << "Sum_Adjacent_Hit:  "<<row <<" " <<col<< " "<<PMT<<endl;
  if (row%2==0)//short row of 8
    {
      if (col == 0)//on the left edge
	{
	  hit_sum = 1 + hit_flag_array[PMT-9]+hit_flag_array[PMT-8]+hit_flag_array[PMT+1]+hit_flag_array[PMT+8]+hit_flag_array[PMT+9];
	  //sum up -9, -8, +1, +8, +9 (excluding -1)
	}
      else if (col == 7)//on the right edge
	{
	  // sum up -9, -8, -1, +8,+9 (excluding +1)
	  hit_sum = 1 + hit_flag_array[PMT-9]+hit_flag_array[PMT-8]+hit_flag_array[PMT-1]+hit_flag_array[PMT+8]+hit_flag_array[PMT+9];
	}
      else // not on an edge 
	{
	  // sum up all of them  -9, -8, -1, 1, 8, 9
	  hit_sum = 1 + hit_flag_array[PMT-9]+hit_flag_array[PMT-8]+hit_flag_array[PMT-1]+hit_flag_array[PMT+1]+hit_flag_array[PMT+8]+hit_flag_array[PMT+9];
	}
    }
  else // long row of 9
    {
      if (col == 0 || col == 8)//on an edge of a long row
	{
	  // don't sum anything. I don't want to cluster find on the long row edges. 
	  hit_sum=0;
	  //cout<< "Sum_Adjacent_Hits returning zero because it is on an outer col"<<endl;
	}
      else // not on an edge
	{
	  hit_sum = 1 + hit_flag_array[PMT-9]+hit_flag_array[PMT-8]+hit_flag_array[PMT-1]+hit_flag_array[PMT+1]+hit_flag_array[PMT+8]+hit_flag_array[PMT+9];
	}
    }
  // cout<<"Sum_Adjacent_Hits returning at end of function: sum= "<< hit_sum << endl;
  return hit_sum;
}
*/

stack <Int_t> Clear_Stack()
{
  stack <Int_t> clearstack; //default constructor is an empty stack
  return clearstack;
}


stack <Int_t> Add_to_Stack(stack <Int_t> inputstack, Int_t PMT)
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
} 


Bool_t Is_NOT_In_Stack(stack <Int_t> checkstack, Int_t PMT)
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
}


Int_t Stack_Size(stack <Int_t> inputstack)
{
  Int_t counter = 0;
  while(!inputstack.empty())
    {
      counter ++;
      inputstack.pop();
    }
  return counter;
}

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
}




void Calculate_PMT_Coord(double& horiz, double& vert, Int_t PMT)
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
}


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

void Test_Reference(double& a, double& b, Int_t test)
{
  a = a+1;
  b = b+2;
  cout<<"test: "<<test<<endl;

  return;
}



void Make_Histos_for_the_Cluster(stack <Int_t> inputstack)
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
      h_grinch_cluster_le_elem -> Fill(PMT,grinch_tdc_le[root_id]);
      h_grinch_cluster_te_elem -> Fill(PMT,grinch_tdc_te[root_id]);
      h_grinch_cluster_tot_elem -> Fill(PMT,grinch_tdc_tot[root_id]);
      h_grinch_cluster_mult_elem -> Fill(PMT,grinch_tdc_mult[root_id]);

      ADC_chan = tdc_to_adc_map_array[PMT];
      if(ADC_chan !=-1)//the PMT has a corresponding ADC channel. 
	{
	  cout<<"hit on PMT "<<PMT<<" with ADC chan "<< ADC_chan<<endl;
	  ADC_root_id = adc_root_index_array[ADC_chan];
	  h_grinch_cluster_atime_elem ->Fill(ADC_chan, adcGAtime[ADC_root_id]);
	  h_grinch_cluster_amp_elem ->Fill(ADC_chan, adcGAmp[ADC_root_id]);
	  h_grinch_cluster_int_elem ->Fill(ADC_chan, grinch_adc[ADC_root_id]);
	  h_grinch_cluster_adcMult_elem ->Fill(ADC_chan,adcGMult[ADC_root_id]); 

	  h_grinch_cluster_adc_tdc ->Fill(adcGAmp[ADC_root_id],grinch_tdc_tot[root_id]); 
	  //cout<<"PMT "<<PMT<<": ADC amp= "<<adcGAmp[ADC_root_id]<<" TDC ToT= "<< grinch_tdc_tot[root_id]<<endl;
	}

    }

  //ADCs? Need to maybe make another array for their id's? 
}



void make_tdc_to_adc_map_array() //get the ADC chan number from the PMT number (if there is one)
{
  //Fill the array so that it is -1 where the PMT index does not have a corresponding TDC channel, and the ADC channel number where it does. 
  for (Int_t i = 0; i<510; i ++)//initialize all to -1
    {           
      tdc_to_adc_map_array[i]= -1;
    }
  
  Int_t ADC_cnt = 0; //ADCs are on PMTs 192 to 255 as of 2/1/2022
  for (Int_t PMT = 192; PMT <= 255; PMT ++)
    {
      tdc_to_adc_map_array[PMT] = ADC_cnt;
      ADC_cnt++;
    } 
}//end make_tdc_to_adc_map_array



// LOAD THE GUI WITH BUTTONS
void show(int option=0)
{

  //cout<<" Option: "<<option<<endl;

  int event;
  int n;
  

  double* adc=new double[N_PMT];
  //  Double_t* tdc=new double[N_PMT][70];
  int* Nadc=new int[N_PMT];
  int* Ntdc=new int[N_PMT];

  
 
  //cout << "TOTAL EVENTS: " << eventEntries << endl;

  //cout << "BEFORE EVENT: " << event_cnt << endl;


  //MENU SELECTION

  if ( option==-2 )   //First Event
    {
      event_cnt = 0;
      fill_event(FORWARD);
    }
  else if ( option==-3 ) //Last Event 
    {
      event_cnt = eventEntries - 1;
      fill_event(BACKWARD);
    }
  else if ( option==0 ) //Ask Event
    {
      Printf("Please Enter Event Number:");
      cin>>usr_evt;
      cout<<"usr_evt "<< usr_evt <<endl;
      if(usr_evt >= 0 && usr_evt < (eventEntries -1)){    
	event_cnt = usr_evt;
	fill_event(FORWARD);     
      }
      else {cout<<"Invalid event number from user input"<<endl;}
    }

  else if ( option==1 ) //Next Event
    {
      //should add a bounds check
      if (event_cnt < eventEntries -1){
	event_cnt++;
	fill_event(FORWARD);
      } else {
	cout<<"You are at the last event in the tree. Can't go farther."<<endl;
      }
    }

  
  else if (option == 6)//Prev event
    {
      if (event_cnt > 0)
	{
	  event_cnt --;
	  fill_event(BACKWARD);	  
	}
      else{
	cout<<"You are at the first entry in the tree. Can't go back farther"<<endl;
      }
     
    }


 
  // cout<<"AFTER EVENT: "<<event_cnt<<endl;

}// END show function
