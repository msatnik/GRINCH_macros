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

// global params

// function declarations 
Double_t gausFit(Double_t *x, Double_t *par);
Double_t gausFit_offset(Double_t *x, Double_t *par);
void WriteToCSV(const string& filename, const vector<double>& data);
void ReadCSVFile(const string& filename, vector<double>& data);
void shiftTH2D(TH2D* hist, vector<double>& shifts);
void FitBinsGausTH2D_vectors(TH2D* hist,vector<double>& amp, vector<double>& amp_err, vector<double>& mean, vector<double>& mean_err,vector<double>& sigma, vector<double>& sigma_err,vector<double>& offset,  vector<double>& offset_err, vector<double>& binCenters);
TH2D* doubleTH2DRange(TH2D* originalHist);


//// **********************************************************************************************////
//// This program aligns the leading edge for the GRINCH tdc signals through fitting and adjusting a
//// TH2D of PMT vs LE that is created by the GRINCH data parsing program. Offsets are read in from 
//// a csv file and applied to each PMT bin in the TH2D. The TH2D with shifted bins is then fit to gaussians
//// and the resulting means are read out to a csv file. It's intended to be used iterativly, first with all the 
//// offsets being zero. Then read in the output as the input and repeate until you like the results. 
////
//// Use the accomanying shell script run_GRINCH_LE_align.sh to do this process automattically. 
//// ./run_GRINCH_LE_align.sh <number of iterations>
//// *********************************************************************************************////
void GRINCH_LE_align(string inputstring = "zeros", string outputstring = "try0"){ // main
  // cout<<"------------------------------------------------------------------------------------------------- "<<endl;
  // cout<<"Give the name of the csv file you want to read in, followed by the name of the file you want it to write to."<<endl;
  // cout<<"Something like:  .x GRINCH_LE_allign.C(``input'', ``output'')"<<endl;
  // cout<<"You will need to manually set the path the rootfile containing the TH2D you want to fit."<<endl;
  // cout<<"------------------------------------------------------------------------------------------------- "<<endl;

  /// set up some gStyle things to make the plots look nicer
  gStyle->SetNumberContours(255); //smooths out the color gradient
  gStyle->SetOptFit(111111);//// Instead of right-clicking on every stat box and typing a million 1's!!!!
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kCool);

  // Load rootfiles for histograms
   TFile *f1 = TFile::Open("../output/sbs8.root"); // Load rootfile
  //TFile *f1 = TFile::Open("../output/sbs9_mar14_size.root"); // Load rootfile
 
 
  //TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/cmjackso/GEn/out_files/GRINCH_GEn2_He3.root"); // this need to be manually set to the rootfile that you are reading the TH2D from. 

  // Load Histograms
  TH2D *hist1 = (TH2D*)f1->Get("h_BBgr_hit_time_elemID_bestcluster"); // Jack's parsing script has it named "pmtnum" and Maria's has it "elemID"

  //// make a clone
  TH2D *clonedHist = (TH2D*)(hist1->Clone("clonedHist"));
  double yMin = clonedHist ->GetYaxis()->GetXmin();
  double yMax = clonedHist ->GetYaxis()->GetXmax();
  
  /// make vectors to store the fit results
  vector<double> amp(510,0);
  vector<double> amp_err(510,0);
  vector<double> mean(510,0);
  vector<double> mean_err(510,0);
  vector<double> sigma(510,0);
  vector<double> sigma_err(510,0);
  vector<double> offset(510,0);  
  vector<double> offset_err(510,0);
  vector<double> binCenters(510,0);
  vector<double> zeros(510,0);

    

  // testing reading in and outputting to csv files.
  const string filepath = "/w/halla-scshelf2102/sbs/msatnik/GRINCH_macros/fitting/csv/"; // Jack change this to your path 
  
  const string inputCSVFile = filepath + inputstring +".csv";
  const string outputCSVFile = filepath + outputstring +".csv";

  

  vector<double> input_offsets;
  ReadCSVFile(inputCSVFile, input_offsets);

 
  cout<<"vector size: "<<input_offsets.size()<<endl;
  if( input_offsets.size() != 510)
    {
      cout << "input vector is not size 510. Check input csv file"<<endl;
    }

  // // print the read in values 
  // for (int i = 0; i< input_offsets.size(); i++){
  //   std::cout<<std::fixed<<std::setprecision(2) << input_offsets[i]<<", ";
  // }
  // cout<<endl;

  //// Check that something was read in. 
  if (!input_offsets.empty()) { 
  
    // Use the function to fit each bin to a gaussian 

    //// displayHist is going to be used for displaying. 
    //// function doubleTH2DRange to double the y range of the histogram.
    TH2D* displayHist = doubleTH2DRange(clonedHist); 
    //// Use the function shiftTH2D to apply the offsets that were read in from the csv file to the displayHist
    shiftTH2D(displayHist, input_offsets); 

    ///// Use the function shfitTH2D to apply the offsets that were read in from the csv file to the original histogram.
    shiftTH2D(clonedHist,input_offsets);


    // Make a histogram of the projection to the Y axis of the shifted histogram. 
    TH1D *h_proj = clonedHist->ProjectionY("h_proj");
    h_proj->SetDirectory(0); //// this prevents the automatic TCanvas that pops up for a projection
    TF1* gaussian_proj =  new TF1("gaussian_proj",gausFit_offset, yMin, yMax,4); //
    gaussian_proj ->SetParameters(1000,0,5,0);
    gaussian_proj ->SetParNames("amp","mean","sigma","offset");
    h_proj->Fit(gaussian_proj, "Q");// "quiet" option
    cout<<"Projecting the histogram to Y-Axis after shifting: "<<endl;
    cout<<"mean = "<<gaussian_proj->GetParameter(1) << " +- "<<gaussian_proj->GetParError(1)<<endl;
    cout<<"sigma = "<<gaussian_proj->GetParameter(2) << " +- "<<gaussian_proj->GetParError(2)<<endl;



    ///// use the function to fit every bin in the TH2D. The function populates the different result vectors. 
    FitBinsGausTH2D_vectors(clonedHist, amp,  amp_err,  mean,  mean_err, sigma,  sigma_err, offset, offset_err, binCenters);
    //// (for some reason, I can't get this function to fit the displayHist that has the larger range. Which I feel would be the more "correct" thing to fit since it's not cutting out any data. )
    
    // convert vectors to arrays suitable for TGraphErrors 
    int nPoints = 510;
    double* xValues =  &binCenters[0];
    double* yValues = &mean[0];
    //double* yErrors = &sigma[0];
    double* yErrors = &zeros[0];
    double* xErrors = &zeros[0];
    
    // make TGraphErrors
    TGraphErrors* graph =  new TGraphErrors(nPoints, xValues, yValues, xErrors, yErrors);
  
    //// draw the histogram before any of the offsets are applied. 
    TCanvas* canvas2 =  new TCanvas("canvas2","2D histo with gaussian fits");
    hist1->Draw("colz");

    TCanvas* canvas3 =  new TCanvas("canvas3","2D histo with gaussian fits");
    h_proj ->Draw();
    //canvas3->SetShowParameters(11111);

    TCanvas* canvas =  new TCanvas("canvas","2D histo with gaussian fits");
    graph ->SetTitle("Vertical Bin Means with Gaussian fit");
    graph ->GetXaxis()->SetTitle("x axis");
    graph ->GetYaxis()->SetTitle("mean");
    graph ->SetMarkerStyle(6);
    displayHist->Draw("colz");
    graph ->Draw("P SAME");

 
    //// loop over the means from the fits to determine the new offset to be read in for the next time. 
    vector<double> result(510,0);
    for (int i = 0; i<510; i++){
      result[i] =input_offsets[i]+mean[i]; // the previous correction plus the new correction 
      //cout<< mean[i]<<", ";
      if(result[i]>= yMax || result[i]<= yMin)
      	{
      	  result[i] = input_offsets[i];
      	}
    }
    cout<<endl;
    

    //// check some results to see what fits might need to be checked manually. 
    /// Jack: looking at the uncertainties as well (or instead) would also probably be good.
    vector<int>pmtcheck;
    int counter = 0;
    cout<<"PMTs that should probaby be checked manually: "<<endl;
    for (int i = 0; i<510; i++){
      if(sigma[i] < 3 || sigma[i] >8 || mean_err[i]>sigma[i] || mean[i]<-10 || mean[i]>10 ||sigma_err[i]>sigma[i])
  	{
  	  cout<<i<<", ";
  	  counter++;
	  pmtcheck.push_back(i);
  	}
      }
    cout<<endl;
    cout<<"Total PMTs to check: "<<counter<<endl;


    ////Write the modified values to a new CSV file
    WriteToCSV(outputCSVFile, result);
    //// Tip: if you overwrite your csv file that's just zeros, replace "result" with "zeros" in this line and run the program again and name the output "zeros". 
    cout << "Modified values written to CSV file successfully." << endl;
  } else {
    cerr << "No data read from the CSV file." << endl;
  }

}// end main



// FUNCTIONCS
Double_t gausFit(Double_t *x, Double_t *par)
{
  Double_t height = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];

  Double_t fit = height * exp(-0.5 * pow((x[0] - mean)/sigma,2));
  return fit;
}
//// end gausFit

Double_t gausFit_offset(Double_t *x, Double_t *par)
{
  Double_t height = par[0];
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t offset = par[3];

  Double_t fit = height * exp(-0.5 * pow((x[0] - mean)/sigma,2)) +offset;
  return fit;
}
//// end gausFit_offset



void FitBinsGausTH2D_vectors(TH2D* orighist,vector<double>& amp, vector<double>& amp_err, vector<double>& mean, vector<double>& mean_err,vector<double>& sigma, vector<double>& sigma_err,vector<double>& offset,  vector<double>& offset_err, vector<double>& binCenters)
//// Fits each vertical bin in a TH2D to a gaussian plus offset. Results returned as vectors. 
{
  // Define parameters
  const int nBinsX = orighist ->GetNbinsX();
  const int nBinsY = orighist ->GetNbinsY();
  const int nBins =  nBinsX;
  const double xmin = orighist ->GetXaxis()->GetXmin();
  const double xmax = orighist ->GetXaxis()->GetXmax();
  const double ymin = orighist ->GetYaxis()->GetXmin();
  const double ymax = orighist ->GetYaxis()->GetXmax();

  // cout<<endl;
  // cout<<"nBins in FitBinsGaus: "<<nBins<<endl;
  // cout<<"ymin "<<ymin<<" ymax: "<<ymax <<endl;
  // cout<<endl;


  //loop through the bins
  for (int i = 0; i<nBins; i++) // 
    {
      // get 1d histo of y projection
      TH1D* binHist =  orighist->ProjectionY("_py",i +1, i+1);
      binHist->SetDirectory(0); //// this prevents the automatic TCanvas that pops up for a projection
      // get bin center
      double binCenter = orighist->GetXaxis()->GetBinCenter(i+1);
      // define fit and gaussian fcn
      TF1* gaussian =  new TF1("gaussian",gausFit_offset,ymin,ymax,4); //
      // give initial guesses and limits on parameters. 
      double mean_prediction = 0;
      gaussian ->SetParameters(1000,0,5,0); 
      gaussian->SetParLimits(1,ymin,ymax); /// limit on the mean
      gaussian ->SetParLimits(2,0.01,10);/// limit on the sigma

      binHist->Fit(gaussian, "Q");// "quiet" option
      // cout<<"bincenter: " <<binCenter<< " mean pred: "<<mean_prediction<<endl;
      //get fit results
      mean[i]= gaussian -> GetParameter(1);
      mean_err[i] =gaussian -> GetParError(1);
      sigma[i] = gaussian ->GetParameter(2);
      sigma_err[i] = gaussian ->GetParError(2);
      amp[i] = gaussian -> GetParameter(0);
      amp_err[i] =  gaussian -> GetParError(0);
      offset[i] = gaussian -> GetParameter(3);
      offset_err[i] = gaussian -> GetParError(3);
      binCenters[i] = binCenter;
     
      //// clean up
      delete gaussian;
    }// end loop over bins

  return;
}
//// end FitBinsGausTH2D_vectors


// Write values from a vector to a CSV file with a new line every 10 values
void WriteToCSV(const string& filename, const vector<double>& data) {
  int counter = 0;
  std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Error: Unable to create file " << filename << endl;
        return;
    }
    //// loop over the vector and write the values to the file 
    for (const auto& value : data) {
      outFile <<std::fixed<<std::setprecision(2)<< value << ",";

      //// comment this if statment out if you don't want the new line every 10 values
      counter++;
      if (counter == 10) 
	{
	  outFile<<endl;
	  counter = 0;
	}
    }

    // Remove the last comma
    outFile.seekp(-1, std::ios_base::end);
    outFile << endl;

    outFile.close();
}// end WriteToCSV


//// Read in values from a csv to a vector
void ReadCSVFile(const string& filename, vector<double>& data) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return;
    }
    //// "tokenizing" the string line-by-line. This will seperate values by the "," token. Should work for any number of lines. 
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string token;
        while (getline(iss, token, ',')) {
            double value;
            istringstream(token) >> value;
            data.push_back(value);
        }
    }
    file.close();
}// end ReadCSVFile


void shiftTH2D(TH2D* hist, vector<double>& shifts)  {
  // Check if the number of shifts matches the number of bins in X axis
  if (shifts.size() != hist->GetNbinsX()) {
    std::cerr << "Number of shifts does not match the number of bins in X axis!" << std::endl;
    return;
  }

  TH2D *tempHist = (TH2D*)(hist->Clone("tempHist"));
  hist->Reset();

  // Loop over each bin in X axis
  for (int i = 1; i <= tempHist->GetNbinsX(); ++i) {
    // Get the content of each bin in the X axis
    for (int j = 1; j <= tempHist->GetNbinsY(); ++j) {
      double content = tempHist->GetBinContent(i, j);
      // Shift the bin index along the y-axis by the corresponding value in the shifts vector
      int new_j = j - shifts[i - 1];
      // If the new index is within the histogram range, set the content to the new index
      if (new_j > 0 && new_j <= tempHist->GetNbinsY()) {
          hist->SetBinContent(i, new_j, content);
      }
    }
  }    
  return;
}
// end shiftTH2D



TH2D* doubleTH2DRange(TH2D* originalHist) {
    // Get original histogram properties
    double xMin = originalHist->GetXaxis()->GetXmin();
    double xMax = originalHist->GetXaxis()->GetXmax();
    int xBins = originalHist->GetXaxis()->GetNbins();
    double yMin = originalHist->GetYaxis()->GetXmin();
    double yMax = originalHist->GetYaxis()->GetXmax();
    int yBins = originalHist->GetYaxis()->GetNbins();

    // Double the range
    double newXMin = 2 * xMin;
    double newXMax = 2 * xMax;
    double newYMin = 2 * yMin;
    double newYMax = 2 * yMax;

    // Create a new histogram with the doubled range
    TH2D* newHist = new TH2D("newHist", originalHist->GetTitle(), xBins, xMin, xMax, 2 * yBins, newYMin, newYMax);

    // Fill the new histogram with content from the original histogram
    for (int i = 1; i <= xBins; ++i) {
        for (int j = 1; j <= yBins; ++j) {
            double content = originalHist->GetBinContent(i, j);
            newHist->Fill(newHist->GetXaxis()->GetBinCenter(i), newHist->GetYaxis()->GetBinCenter(j)+(yMax-yMin)/2, content);
        }
    }

    //cout<<"new yMin: "<<newHist->GetYaxis()->GetXmin()<<endl;
    //cout<<"new yMax: "<<newHist->GetYaxis()->GetXmax()<<endl;

    return newHist;
}/// end doubleTH2DRange
