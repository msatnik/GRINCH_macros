#include <TH1D.h>
#include <TAttFill.h>

// Global Params
Double_t yMin = -7;
Double_t yMax= 10;
Double_t xMin = 8 ;
Double_t xMax= 16;

// Function Declarations 
TGraphErrors* FitBinsGausTH2D(TH2D* hist);


void Plot_tw(){// MAIN
  
  TCanvas *c1 = new TCanvas("c1","GRINCH Hit branches",1000,500);
  c1->Divide(2,1);
  c1->SetLogy();

  //for( Int_t i=0; i<nkine; i++ ){

    c1->cd(1);

    gStyle->SetOptStat(0);
    //Get pre-ecal energy dist
    TFile *f1 = TFile::Open("../output/sbs8_sept22.root");
    TH1D *h1= (TH1D*)f1->Get("h_BBgr_hit_time");
    TH1D *h1_norm =  (TH1D*)(h1 ->Clone("h1_norm"));
    h1_norm ->Scale(1./h1_norm->Integral(),"width");
    //h1_norm->GetXaxis()->SetRange(0,20);
    h1_norm->SetLineWidth(3);
    h1_norm ->SetFillColor(30);
    h1_norm ->SetFillStyle(3003);
    h1_norm->SetLineColor(30);
    //h1->GetXaxis()->SetRange(0,20);
    h1 ->SetLineWidth(2);
    h1->SetFillColor(30);
    h1 ->SetFillStyle(3003);
    h1 ->SetLineColor(30);
    
    
  
    TH1D *h2 = (TH1D*)f1->Get("h_BBgr_hit_time_tw_poly");
    TH1D *h2_norm =  (TH1D*)(h2 ->Clone("h2_norm"));
    h2_norm ->Scale(1./h2_norm->Integral(),"width");
    // h2_norm->GetXaxis()->SetRange(0,20);
    h2_norm->SetLineWidth(2);
    h2_norm->SetLineColor(kMagenta - 6);
    h2_norm ->SetFillColor(kMagenta - 6);
    h2_norm ->SetFillStyle(3003);
    //h2 ->GetXaxis()->SetRange(0,20);
    h2 ->SetLineWidth(2);
    h2 ->SetLineColor(kMagenta - 6);
    h2 ->SetFillColor(kMagenta - 6);
    h2  ->SetFillStyle(3003);

    TH1D *h3 = (TH1D*)f1->Get("h_BBgr_hit_time_tw_trad");
    TH1D *h3_norm =  (TH1D*)(h3 ->Clone("h3_norm"));
    h3_norm ->Scale(1./h3_norm->Integral(),"width");
    // h2_norm->GetXaxis()->SetRange(0,20);
    h3_norm->SetLineWidth(2);
    h3_norm->SetLineColor(kBlue - 6);
    h3_norm ->SetFillColor(kBlue - 6);
    h3_norm ->SetFillStyle(3003);
    //h2 ->GetXaxis()->SetRange(0,20);
    h3 ->SetLineWidth(2);
    h3 ->SetLineColor(kBlue - 6);
    h3 ->SetFillColor(kBlue - 6);
    h3  ->SetFillStyle(3003);
    


    h1->Draw("hist");
    h2->Draw("hist same");
    //h3->Draw("hist same");
   

    // //Add a legend
    // auto legend = new TLegend(0.43,0.7,0.89,0.89);
    // legend->SetTextSize(0.035);
    // legend->SetHeader("SH + PS Energy");
    // legend->AddEntry(h1,"SBS8","l");
    // legend->AddEntry(h2,"SBS9","l");
    // // legend->AddEntry(h3_norm,"No Offsets","l");
    // legend->Draw();


    
      TH2D *origHist = (TH2D*)f1->Get("h_BBgr_hit_time_amp_tw_poly");

      TH2D* clonedHist3 = new TH2D("clonedHist3", origHist->GetTitle(),
                                origHist->GetXaxis()->GetNbins(),
                                origHist->GetXaxis()->GetXmin(),
                                origHist->GetXaxis()->GetXmax(),
                                origHist->GetYaxis()->FindBin(yMax) - origHist->GetYaxis()->FindBin(yMin) + 1,
                                yMin, yMax);

 // Loop over the X and Y bins of the original TH2D
    for (int i = 1; i <= origHist->GetXaxis()->GetNbins(); i++) {
        for (int j = origHist->GetYaxis()->FindBin(yMin); j <= origHist->GetYaxis()->FindBin(yMax); j++) {
            double binContent = origHist->GetBinContent(i, j);
            clonedHist3->SetBinContent(i, j - origHist->GetYaxis()->FindBin(yMin) + 1, binContent);
        }
    }


    TCanvas* canvas =  new TCanvas("canvas","2D histo with gaussian fits");
    // make a default TGraphErrors
    TGraphErrors* graph =  new TGraphErrors();
    graph = FitBinsGausTH2D(clonedHist3); // function fits each bin to a gaussian and returns a graph o f the means with the sigma as the error bars
    // customize graph
    graph ->SetTitle("Vertical Bin Means with Gaussian fit");
    graph ->GetXaxis()->SetTitle("x axis");
    graph ->GetYaxis()->SetTitle("mean");
    graph ->SetMarkerStyle(20);
    graph ->SetMarkerSize(1.2);

    //draw graph
    clonedHist3->Draw("colz");
    graph ->Draw("P SAME");
    canvas ->Update();

}// end MAIN


TGraphErrors* FitBinsGausTH2D(TH2D* hist)
//// Fits each vertical bin in a TH2D to a gaussian and returns a TGraphErros with the means and the sigma as the vertical error bars. 
{
  // Fit each vertical bin to a gaussian and plot the means on top of the 2D histo. 
  // Define parameters
  const int nBinsX = hist ->GetNbinsX();
  const int nBinsY = hist ->GetNbinsY();
  const int nBins =  nBinsX;
  const double xmin = hist ->GetXaxis()->GetXmin();
  const double xmax = hist ->GetXaxis()->GetXmax();

  // cout<<endl;
  // cout<<"nBins: "<<nBins<<endl;
  // cout<<endl;

  // arrays to store fit results
  std::vector<double> binCenters(nBins);
  std::vector<double> means(nBins);
  std::vector<double> sigmas(nBins);
  std::vector<double> half_sigmas(nBins);
  std::vector<double> zeros(nBins);

  //loop through the bins
  for (int i = 0; i<nBins; i++) // 
    {
      // get 1d histo of y projection
      TH1D* binHist =  hist->ProjectionY("_py",i +1, i+1);
      // get bin center
      double binCenter = hist->GetXaxis()->GetBinCenter(i+1);
      // define fit and gaussian fcn
      TF1* gaussian =  new TF1("gaussian","gaus",xmin,xmax); // this uses the default "gaus" from root. Could potentially be worth it to make your own gaus function as sometimes the default can be weird.
      binHist->Fit(gaussian, "Q");// "quiet" option
      //get fit results
      double mean = gaussian -> GetParameter(1);
      double sigma = gaussian ->GetParameter(2);
      //store results 
      binCenters[i] = binCenter;
      means[i] = mean;
      sigmas[i] = sigma;
      half_sigmas[i] = sigma*0.5;
      zeros[i] = 0;
      // cout<<"bin center: "<<binCenter<<endl;
      // cout<<"mean: "<<mean<<endl;
      // cout<<"sigma: "<<sigma<<endl;
      // cout<<endl;
      //// clean up
      delete gaussian;
    }// end loop over bins

  // convert vectors to arrays suitable for TGraphErrors 
  int nPoints = nBins;
  double* xValues =  &binCenters[0];
  double* yValues = &means[0];
  double* yErrors = &half_sigmas[0];
  double* xErrors = &zeros[0];
    
  // make TGraphErrors
  TGraphErrors* graph =  new TGraphErrors(nPoints, xValues, yValues, xErrors, yErrors);
    
  return graph;
} // end TGraphErrors* FitBinsGausTH2D(TH2D* hist)
