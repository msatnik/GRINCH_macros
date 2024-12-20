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

#include "/work/halla/sbs/msatnik/GMn/classes/CSVHandler.h"
#include "/work/halla/sbs/msatnik/GMn/classes/Utility.h"
#include "/work/halla/sbs/msatnik/GMn/classes/Utility.cpp"

// global params

// function declarations 
Double_t gausFit(Double_t *x, Double_t *par);
Double_t gausFit_offset(Double_t *x, Double_t *par);
void WriteToCSV(const string& filename, const vector<double>& data);
void ReadCSVFile(const string& filename, vector<double>& data);
void shiftTH2D(TH2D* hist, vector<double>& shifts);
void FitBinsGausTH2D_vectors(TH2D* hist,vector<double>& amp, vector<double>& amp_err, vector<double>& mean, vector<double>& mean_err,vector<double>& sigma, vector<double>& sigma_err,vector<double>& offset,  vector<double>& offset_err, vector<double>& binCenters);
TH2D* doubleTH2DRange(TH2D* originalHist);

// Hard coding the number of entries from the original root file for now
double nEntries =6104665;


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
void GRINCH_LE_align_v2(string inputstring = "zeros", string outputstring = "try0"){ // main
  // cout<<"------------------------------------------------------------------------------------------------- "<<endl;
  // cout<<"Give the name of the csv file you want to read in, followed by the name of the file you want it to write to."<<endl;
  // cout<<"Something like:  .x GRINCH_LE_allign.C(``input'', ``output'')"<<endl;
  // cout<<"You will need to manually set the path the rootfile containing the TH2D you want to fit."<<endl;
  // cout<<"------------------------------------------------------------------------------------------------- "<<endl;

  /// set up some gStyle things to make the plots look nicer
  gStyle->SetNumberContours(255); //smooths out the color gradient
  gStyle->SetOptFit(111111);//// Instead of right-clicking on every stat box and typing a million 1's!!!!
  // gStyle->SetOptStat(0);
  gStyle->SetPalette(kCool);

  // Load rootfiles for histograms
  TFile *f1 = TFile::Open("/lustre24/expphy/volatile/halla/sbs/msatnik/output/grinch_output/grinch_dec_20/grinch_sbs9_dec_20.root");// Load rootfile
  //TFile *f1 = TFile::Open("/lustre24/expphy/volatile/halla/sbs/msatnik/output/grinch_output/sbs8.root"); // Load rootfile
  //TFile *f1 = TFile::Open("../output/sbs9_mar14_size.root"); // Load rootfile
 
 
  //TFile *f1 = TFile::Open("/w/halla-scshelf2102/sbs/cmjackso/GEn/out_files/GRINCH_GEn2_He3.root"); // this need to be manually set to the rootfile that you are reading the TH2D from. 

  // Load Histograms
  TH2D *hist1 = (TH2D*)f1->Get("h_BBgr_hit_time_elemID"); // Jack's parsing script has it named "pmtnum" and Maria's has it "elemID"
  // _bestcluster


  TParameter<int> *counterParam = (TParameter<int>*) f1->Get("normalizer_counter");
  int normalizer_counter = counterParam->GetVal();
  std::cout << "Final counter value: " << normalizer_counter << std::endl;
  

  //// make a clone
  TH2D *clonedHist = (TH2D*)(hist1->Clone("clonedHist"));
  double yMin = clonedHist ->GetYaxis()->GetXmin();
  double yMax = clonedHist ->GetYaxis()->GetXmax();

  TH1D* meanHist = new TH1D("meanHist",";mean (ns)",100,-10,10);
  TH1D* sigmaHist = new TH1D("sigmaHist",";sigma (ns)",100,0,10);

  TH1D* meanHist_good_fits = new TH1D("meanHist_good_fits",";mean LE for good fits only (ns)",100,-5,5);
  TH1D* sigmaHist_good_fits = new TH1D("sigmaHist_good_fits",";sigma for good fits only (ns)",100,0,10);

  TH1D* rateHist = new TH1D("rateHist",";background rate (kHz)",50,0,50);

  TH1D* rateHist_offset = new TH1D("rateHist_offset",";background rate (kHz)",50,0,50);

  
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

    

  // reading in and outputting to csv files.
  const string filepath = "/w/halla-scshelf2102/sbs/msatnik/GRINCH_macros/fitting/csv/"; // Jack change this to your path 
  
  const string inputCSVFile = filepath + inputstring +".csv";
  const string outputCSVFile = filepath + outputstring +".csv";

  

  vector<double> input_offsets;
  CSVHandler::ReadFromCSV(inputCSVFile, input_offsets);

  // class to handle various functions 
  Utility utilityHandler;
  

 
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
    double offset_guess = h_proj->GetBinContent(5); // get the value on bin#5 which is the wings
    double amp_guess = offset_guess; // ehhh they're prob the same order of magnitude. 
    gaussian_proj ->SetParameters(amp_guess,0,5,offset_guess);
    gaussian_proj ->SetParNames("amp","mean","sigma","offset");
    h_proj->Fit(gaussian_proj, "Q");// "quiet" option
    cout<<"Projecting the histogram to Y-Axis after shifting: "<<endl;
    cout<<"mean = "<<gaussian_proj->GetParameter(1) << " +/- "<<gaussian_proj->GetParError(1)<<endl;
    cout<<"sigma = "<<gaussian_proj->GetParameter(2) << " +/- "<<gaussian_proj->GetParError(2)<<endl;



    ///// use the function to fit every bin in the TH2D. The function populates the different result vectors. 
    FitBinsGausTH2D_vectors(clonedHist, amp,  amp_err,  mean,  mean_err, sigma,  sigma_err, offset, offset_err, binCenters);
    //// (for some reason, I can't get this function to fit the displayHist that has the larger range. Which I feel would be the more "correct" thing to fit since it's not cutting out any data. )


    // cout<<"From Gaussian + offset fit for PMT 0: "<<endl;
    // cout<< "amp: "<<amp[0]<< ", mean: "<<mean[0] << ", sigma: "<<sigma[0]<<", offset: "<<offset[0]<<endl<<endl;

    //   cout<<"From Gaussian + offset fit for PMT 1: "<<endl;
    // cout<< "amp: "<<amp[1]<< ", mean: "<<mean[1] << ", sigma: "<<sigma[1]<<", offset: "<<offset[1]<<endl<<endl;
    

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

    
    std::vector<double> three_sigma = utilityHandler.MultiplyVectorByScalar(sigma,3); // make a vector of sigma*3
    std::vector<double> five_sigma = utilityHandler.MultiplyVectorByScalar(sigma,5); // make a vector of sigma*5
    
    // convert vectors to arrays suitable for TGraphErrors 
    int nPoints = 510;
    double* xValues =  &binCenters[0];
    double* yValues = &mean[0];
    double* yErrors = &three_sigma[0];
    //double* yErrors = &sigma[0];
    // double* yErrors = &zeros[0];
    double* xErrors = &zeros[0];
    
    // make TGraphErrors
    TGraphErrors* graph =  new TGraphErrors(nPoints, xValues, yValues, xErrors, yErrors);

    // Make an overlay to change the color of the markers that we may want to check
    TPolyMarker* overlay = utilityHandler.CreatePointColorOverlay(graph, pmtcheck, kRed);

    TGraphErrors* sigmaGraph = new TGraphErrors(nPoints, xValues, &sigma[0], &zeros[0], &zeros[0]);

    // Make an overlay to change the color of the markers that we may want to check
    TPolyMarker* sigmaOverlay = utilityHandler.CreatePointColorOverlay(sigmaGraph, pmtcheck, kRed);
  
    //// draw the histogram before any of the offsets are applied. 
    TCanvas* canvas2 =  new TCanvas("canvas2","2D histo with gaussian fits");
    hist1->Draw("colz");

    TCanvas* canvas3 =  new TCanvas("canvas3","2D histo with gaussian fits");
    h_proj ->Draw();
    //canvas3->SetShowParameters(11111);

    TCanvas* canvas =  new TCanvas("canvas","2D histo with gaussian fits");
    graph ->SetTitle("Vertical Bin Means with Gaussian fit");
    graph ->GetXaxis()->SetTitle("PMT");
    graph ->GetYaxis()->SetTitle("mean (ns)");
    graph ->SetMarkerStyle(20);
    graph ->SetMarkerSize(1);
    displayHist->Draw("colz");
    graph ->Draw("P SAME");
    if (overlay){
      overlay->Draw("same");
    }

    // fill the histogram with the mean values 
    utilityHandler.FillHistogramFromVector(mean, meanHist);
    // fill the histogram with the sigma values
    utilityHandler.FillHistogramFromVector(sigma, sigmaHist);

    // Function to remove elements from a vector based on a list of indices
    std::vector<double> mean_good_fits  = utilityHandler.RemoveElementsByIndices(mean, pmtcheck);
    // fill the histogram with the good mean values 
    utilityHandler.FillHistogramFromVector(mean_good_fits, meanHist_good_fits);
     
    std::vector<double> sigma_good_fits  = utilityHandler.RemoveElementsByIndices(sigma, pmtcheck);
    utilityHandler.FillHistogramFromVector(sigma_good_fits, sigmaHist_good_fits);

     
    TCanvas* canvasMean_good_fits = new TCanvas("canvasMean_good_fits","mean of each channel good fit");
    meanHist_good_fits->Draw();

    TCanvas* canvasSigma_good_fits = new TCanvas("canvasSigma_good_fits","sigma of each channel good fit");
    sigmaHist_good_fits->Draw();


    TCanvas* canvasMean = new TCanvas("canvasMean","mean of each channel");
    meanHist->Draw();

    TCanvas* canvasSigma = new TCanvas("canvasSigma","sigma of each channel");
    sigmaHist->Draw();

     
    TCanvas* canvasSigmaGraph = new TCanvas("canvasSigmaGraph","sigma of each channel");


    utilityHandler.customizeGraphMore(sigmaGraph, 8, 4 , 1, 
				      "sigma from gaussian fit (ns)","PMT", "Sigma (ns)",1,1);
    sigmaGraph ->SetTitle("sigma from gaussian fit (ns)");
    sigmaGraph->Draw("A P");
    sigmaOverlay->Draw("same");

     
    // get the backgorund rates from each histogram by looking at the wings
    std::vector<double> rates_not_normalized = utilityHandler.GetRatesFromTH2D(clonedHist,-30,-20);

    // we now need to normalize this by the total number of events that we processed 
    std::vector<double> rates_ns =utilityHandler.DivideVectorByScalar(rates_not_normalized,normalizer_counter);

    std::vector<double> rates_Hz = utilityHandler.MultiplyVectorByScalar(rates_ns, 1E9); // convert from 1/(ns) to 1/s

    std::vector<double> rates_kHz = utilityHandler.DivideVectorByScalar(rates_Hz, 1E3); // convert from Hz to kHz

    std::vector<double> rates_mHz = utilityHandler.DivideVectorByScalar(rates_Hz, 1E6); // convert from Hz to mHz


    // Alternative way to calculate the rates using the offset from the fit
    std::vector<std::pair<double, double>> ratesAndErrors_not_normalized = utilityHandler.GetRatesAndErrorsFromTH2D(clonedHist,offset,offset_err);
      
    std::vector<std::pair<double, double>> ratesAndErrors_ns = utilityHandler.DividePairVectorByScalar(ratesAndErrors_not_normalized, normalizer_counter);

    std::vector<std::pair<double, double>> ratesAndErrors_Hz = utilityHandler.MultiplyPairVectorByScalar(ratesAndErrors_ns, 1E9); // convert from 1/(ns) to 1/s
     
    std::vector<std::pair<double, double>> ratesAndErrors_kHz = utilityHandler.DividePairVectorByScalar(ratesAndErrors_Hz, 1E3); // convert from Hz to kHz

    std::vector<std::pair<double, double>> ratesAndErrors_mHz = utilityHandler.DividePairVectorByScalar(ratesAndErrors_Hz, 1E6); // convert from Hz to mHz

    // Extract rates (first) and errors (second) into separate vectors
    std::vector<double> rates_kHz_offset = utilityHandler.ExtractPairElement(ratesAndErrors_kHz, true);
    std::vector<double> errors_kHz_offset = utilityHandler.ExtractPairElement(ratesAndErrors_kHz, false);

    utilityHandler.FillHistogramFromVector(rates_kHz_offset, rateHist_offset);

    TGraphErrors* ratesGraph_offset = new TGraphErrors(nPoints, &binCenters[0] , &rates_kHz_offset[0], &zeros[0],  &zeros[0]);
      
    // Make an overlay to change the color of the markers that we may want to check
    TPolyMarker* ratesOverlay_offset = utilityHandler.CreatePointColorOverlay(ratesGraph_offset, pmtcheck, kRed);
      

    TCanvas* canvasRates_offset = new TCanvas("canvasRates_offset","background rates offset method");
    rateHist_offset->Draw();

    TCanvas* ratesGraphCanvas_offset= new TCanvas("ratesGraphCanvas_offset","background rates offset method");
     
    utilityHandler.customizeGraphMore(ratesGraph_offset, 8, 4 , 1, 
				      "GRINCH PMT Background Rates (kHz)","PMT", "Background Rate (kHz)",1,1.5);
    ratesGraph_offset->Draw("A P");
    ratesOverlay_offset->Draw("same");
    // overlay_rates ->Draw("same");
      

    // fill the histogram with the rates
    utilityHandler.FillHistogramFromVector(rates_kHz, rateHist);

    TGraphErrors* ratesGraph = new TGraphErrors(nPoints, &binCenters[0] , &rates_kHz[0], &zeros[0], &zeros[0]);

    TPolyMarker* overlay_rates = utilityHandler.CreatePointColorOverlay(ratesGraph, pmtcheck, kRed);

    TCanvas* ratesGraphCanvas= new TCanvas("ratesGraphCanvas","background rates");
     
    utilityHandler.customizeGraphMore(ratesGraph, 8, 4 , 1, 
				      "GRINCH PMT Background Rates (kHz)","PMT", "Background Rate (kHz)",1,1.5);
    ratesGraph->Draw("A P");
    overlay_rates ->Draw("same");

    TCanvas* canvasRates = new TCanvas("canvasRates","background rates");
    rateHist->Draw();



   

 
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
    



   
  

    ////Write the modified values to a new CSV file
    //WriteToCSV(outputCSVFile, result);
    CSVHandler::WriteToCSV(outputCSVFile, result, 10, 2);
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

void FitBinsGausTH2D_vectors(TH2D* orighist, vector<double>& amp, vector<double>& amp_err, vector<double>& mean, vector<double>& mean_err, vector<double>& sigma, vector<double>& sigma_err, vector<double>& offset, vector<double>& offset_err, vector<double>& binCenters) {
  Utility utilityHandler;

  const int nBins = orighist->GetNbinsX();
  const double ymin = orighist->GetYaxis()->GetXmin();
  const double ymax = orighist->GetYaxis()->GetXmax();

  // Resize vectors to hold fit results
  amp.resize(nBins, 0);
  amp_err.resize(nBins, 0);
  mean.resize(nBins, 0);
  mean_err.resize(nBins, 0);
  sigma.resize(nBins, 0);
  sigma_err.resize(nBins, 0);
  offset.resize(nBins, 0);
  offset_err.resize(nBins, 0);
  binCenters.resize(nBins, 0);
    
  cout<<"nBins "<<nBins<<endl;

  for (int i = 0; i < nBins; i++) {
    // Get 1D projection
    std::string binHistName = "_py_" + std::to_string(i);
    TH1D* binHist = orighist->ProjectionY(binHistName.c_str(), i + 1, i + 1);
    binHist->SetDirectory(0);

    // Define fit function
    std::string funcName = "gaussian_bin_" + std::to_string(i);
    TF1* gaussian = new TF1(funcName.c_str(), gausFit_offset, ymin, ymax, 4);
    gaussian->SetParNames("amp (counts)", "mean (ns)", "sigma (ns)", "p0 (counts)");
    gaussian->SetParameters(1000, 0, 5, 0);
    gaussian->SetParLimits(1, ymin, ymax);
    gaussian->SetParLimits(2, 0.01, 10);

    // Fit and get results
    binHist->Fit(gaussian, "Q");
    mean[i] = gaussian->GetParameter(1);
    mean_err[i] = gaussian->GetParError(1);
    sigma[i] = gaussian->GetParameter(2);
    sigma_err[i] = gaussian->GetParError(2);
    amp[i] = gaussian->GetParameter(0);
    amp_err[i] = gaussian->GetParError(0);
    offset[i] = gaussian->GetParameter(3);
    offset_err[i] = gaussian->GetParError(3);
    binCenters[i] = orighist->GetXaxis()->GetBinCenter(i + 1);

    // Draw the first fit for validation
    if (i == 250) {
      TH1D *clonedBinHist = (TH1D*)(binHist->Clone("clonedBinHist"));
	  
      std::cout << "i = " << i << ", amp[i]=" << amp[i] 
		<< ", mean[i]=" << mean[i] 
		<< ", sigma[i]=" << sigma[i] 
		<< ", offset[i]=" << offset[i] << std::endl;

      std::cout << "amp: " << gaussian->GetParameter(0) 
		<< ", mean: " << gaussian->GetParameter(1)
		<< ", sigma: " << gaussian->GetParameter(2) 
		<< ", offset: " << gaussian->GetParameter(3) 
		<< std::endl;
	    
      utilityHandler.DrawGaussianPlusPol0(clonedBinHist,
					  gaussian->GetParameter(0),
					  gaussian->GetParameter(1),
					  gaussian->GetParameter(2),
					  gaussian->GetParameter(3));

	   
    
      // TCanvas* looptestcanvas = new TCanvas("looptestcanvas", "ProjectionY Fit", 800, 600);
      // clonedBinHist->Draw(); // Explicitly draw the first histogram
      // gaussian->Draw("same"); // Overlay the fit on top
   
    }

    // Remove the fit function from the histogram
    //binHist->GetListOfFunctions()->Clear();

    delete gaussian;
  }
}// end FitBinsGausTH2D_vectors


// void FitBinsGausTH2D_vectors(TH2D* orighist,vector<double>& amp, vector<double>& amp_err, vector<double>& mean, vector<double>& mean_err,vector<double>& sigma, vector<double>& sigma_err,vector<double>& offset,  vector<double>& offset_err, vector<double>& binCenters)
// //// Fits each vertical bin in a TH2D to a gaussian plus offset. Results returned as vectors. 
// {

//   Utility utilityHandler;
  
//   // Define parameters
//   const int nBinsX = orighist ->GetNbinsX();
//   const int nBinsY = orighist ->GetNbinsY();
//   const int nBins =  nBinsX;
//   const double xmin = orighist ->GetXaxis()->GetXmin();
//   const double xmax = orighist ->GetXaxis()->GetXmax();
//   const double ymin = orighist ->GetYaxis()->GetXmin();
//   const double ymax = orighist ->GetYaxis()->GetXmax();

//   // cout<<endl;
//   // cout<<"nBins in FitBinsGaus: "<<nBins<<endl;
//   // cout<<"ymin "<<ymin<<" ymax: "<<ymax <<endl;
//   // cout<<endl;

//   amp.resize(nBins, 0);
//   amp_err.resize(nBins, 0);
//   mean.resize(nBins, 0);
//   mean_err.resize(nBins, 0);
//   sigma.resize(nBins, 0);
//   sigma_err.resize(nBins, 0);
//   offset.resize(nBins, 0);
//   offset_err.resize(nBins, 0);
//   binCenters.resize(nBins, 0);



//   //loop through the bins
//   for (int i = 0; i<nBins; i++) // 
//     {
//       // get 1d histo of y projection
//       TH1D* binHist =  orighist->ProjectionY("_py",i +1, i+1);
//       binHist->SetDirectory(0); //// this prevents the automatic TCanvas that pops up for a projection
//       // get bin center
//       double binCenter = orighist->GetXaxis()->GetBinCenter(i+1);
//       // define fit and gaussian fcn
//       TF1* gaussian =  new TF1("gaussian",gausFit_offset,ymin,ymax,4); //
//       // give initial guesses and limits on parameters. 
//       double mean_prediction = 0;
//       gaussian ->SetParameters(1000,0,5,0); 
//       gaussian->SetParLimits(1,ymin,ymax); /// limit on the mean
//       gaussian ->SetParLimits(2,0.01,10);/// limit on the sigma

//       binHist->Fit(gaussian, "Q");// "quiet" option
//       // cout<<"bincenter: " <<binCenter<< " mean pred: "<<mean_prediction<<endl;
//       //get fit results
//       mean[i]= gaussian -> GetParameter(1);
//       mean_err[i] =gaussian -> GetParError(1);
//       sigma[i] = gaussian ->GetParameter(2);
//       sigma_err[i] = gaussian ->GetParError(2);
//       amp[i] = gaussian -> GetParameter(0);
//       amp_err[i] =  gaussian -> GetParError(0);
//       offset[i] = gaussian -> GetParameter(3);
//       offset_err[i] = gaussian -> GetParError(3);
//       binCenters[i] = binCenter;

//       if (i == 0){
// 	std::cout << "amp[0]: " << gaussian->GetParameter(0) 
//               << ", mean[0]: " << gaussian->GetParameter(1)
//               << ", sigma[0]: " << gaussian->GetParameter(2) 
//               << ", offset[0]: " << gaussian->GetParameter(3) 
//               << std::endl;
// 	TH1D* binHist0 = (TH1D*)(binHist->Clone("binHist0"));
// 	utilityHandler.DrawGaussianPlusPol0(binHist0 ,gaussian -> GetParameter(0), gaussian -> GetParameter(1),gaussian ->GetParameter(2),gaussian -> GetParameter(3));
//       }
     
//       //// clean up
//       delete gaussian;
//     }// end loop over bins

//   return;
// }
// //// end FitBinsGausTH2D_vectors


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
  newHist ->GetXaxis()->SetTitle("PMT");
  newHist->GetYaxis()->SetTitle("mean (ns)");

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
