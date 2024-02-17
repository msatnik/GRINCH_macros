void Plot_W2(){
  
  TCanvas *c1 = new TCanvas("c1","GRINCH",1000,500);
  c1->Divide(2,1);
  c1->SetLogy();


    c1->cd(1);

    gStyle->SetOptStat(0);

    TFile *f1 = TFile::Open("../output/sbs8_nov24_huge.root");

    
    TH1D *h_W2= (TH1D*)f1->Get("h_W2");  
    TH1D *h_W2_elastic =(TH1D*)f1->Get("h_W2_elastic");  

    TH1D *h_W2_gr_anticut= (TH1D*)f1->Get("h_W2_gr_anticut");   
    TH1D *h_W2_ps_anticut = (TH1D*)f1->Get("h_W2_ps_anticut");
    TH1D *h_W2_gr_ps_anticut = (TH1D*)f1->Get("h_W2_gr_ps_anticut");

    TH1D *h_W2_gr_cut= (TH1D*)f1->Get("h_W2_gr_cut");   
    TH1D *h_W2_ps_cut = (TH1D*)f1->Get("h_W2_ps_cut");
    TH1D *h_W2_gr_ps_cut = (TH1D*)f1->Get("h_W2_gr_ps_cut");

    TH1D *h_W2_allclus_cut= (TH1D*)f1->Get("h_W2_allclus_cut"); 
    TH1D *h_W2_allclus_anticut= (TH1D*)f1->Get("h_W2_allclus_anticut");   

    h_W2 ->SetLineWidth(2);
    h_W2 ->SetLineColor(kRed- 6);

    h_W2_elastic ->SetLineWidth(2);
    h_W2_elastic ->SetLineColor(kViolet -6);

    h_W2_gr_anticut->SetLineWidth(2);
    h_W2_gr_anticut ->SetLineColor(30);
 
    h_W2_ps_anticut ->SetLineWidth(2);
    h_W2_ps_anticut ->SetLineColor(kMagenta - 6);
    
    h_W2_gr_ps_anticut ->SetLineWidth(2);
    h_W2_gr_ps_anticut  ->SetLineColor(kBlue - 6);


    h_W2_gr_cut->SetLineWidth(2);
    h_W2_gr_cut ->SetLineColor(30);
 
    h_W2_ps_cut ->SetLineWidth(2);
    h_W2_ps_cut ->SetLineColor(kMagenta - 6);
    
    h_W2_gr_ps_cut ->SetLineWidth(2);
    h_W2_gr_ps_cut  ->SetLineColor(kBlue - 6);

    h_W2_allclus_cut ->SetLineWidth(2);
    h_W2_allclus_cut ->SetLineColor(kYellow -6); 

    h_W2_allclus_anticut ->SetLineWidth(2);
    h_W2_allclus_anticut ->SetLineColor(kYellow -6); 

    // Double_t maxBin_allclus_anticut, maxBinValue_allclus_anticut,  maxBinXValue_allclus_anitcut,maxBin_gr_anticut,maxBinValue_gr_anticut,maxBinXValue_gr_anitcut ;

    //// Get max bin of the allclus anticut
    double maxBin_allclus_anticut = h_W2_allclus_anticut -> GetMaximumBin();
    double maxBinValue_allclus_anticut = h_W2_allclus_anticut -> GetBinContent(maxBin_allclus_anticut);
    double maxBinXValue_allclus_anticut = h_W2_allclus_anticut->GetXaxis()->GetBinCenter(maxBin_allclus_anticut);
    cout<<"maxBin_allclus_anticut: "<<maxBin_allclus_anticut<<endl;
    cout<<"maxBinValue_allclus_anticut: "<<maxBinValue_allclus_anticut<<endl;
    cout<<"maxBinXValue_allclus_anticut: "<<maxBinXValue_allclus_anticut<<endl;

    //// Get max bin of the gr anticut
    double maxBin_gr_anticut = h_W2_gr_anticut -> GetMaximumBin();
    double maxBinValue_gr_anticut = h_W2_gr_anticut -> GetBinContent(maxBin_gr_anticut);
    double maxBinXValue_gr_anticut = h_W2_allclus_anticut->GetXaxis()->GetBinCenter(maxBin_gr_anticut);
    cout<<"maxBin_gr_anticut: "<<maxBin_gr_anticut<<endl;
    cout<<"maxBinValue_gr_anticut: "<<maxBinValue_gr_anticut<<endl;
    cout<<"maxBinXValue_gr_anticut: "<<maxBinXValue_gr_anticut<<endl;

    Double_t scalefactor_gr_allclus = maxBinValue_gr_anticut/maxBinValue_allclus_anticut ;
    //Double_t scalefactor_gr_allclus = 

    TH1D *h_W2_allclus_anticut_scaled =  (TH1D*)(h_W2_allclus_anticut->Clone("h_W2_allclus_anticut_scaled"));
    h_W2_allclus_anticut_scaled->Scale(scalefactor_gr_allclus);
    cout<<"scalefactor for alllclus anticut "<< scalefactor_gr_allclus <<endl;

     //// Get max bin of the gr ps anticut
    double maxBin_gr_ps_anticut = h_W2_gr_ps_anticut -> GetMaximumBin();
    double maxBinValue_gr_ps_anticut = h_W2_gr_ps_anticut -> GetBinContent(maxBin_gr_ps_anticut);
    double maxBinXValue_gr_ps_anticut = h_W2_gr_ps_anticut->GetXaxis()->GetBinCenter(maxBin_gr_ps_anticut);
    cout<<"maxBin_gr_ps_anticut: "<<maxBin_gr_ps_anticut<<endl;
    cout<<"maxBinValue_gr_ps_anticut: "<<maxBinValue_gr_ps_anticut<<endl;
    cout<<"maxBinXValue_gr_ps_anticut: "<<maxBinXValue_gr_ps_anticut<<endl;

    Double_t scalefactor_gr_ps = maxBinValue_gr_anticut/maxBinValue_gr_ps_anticut;
    
    TH1D *h_W2_gr_ps_anticut_scaled =  (TH1D*)(h_W2_gr_ps_anticut->Clone("h_W2_gr_ps_scaled"));
    h_W2_gr_ps_anticut_scaled->Scale(scalefactor_gr_ps);
    cout<<"scalefactor for gr_ps anticut "<< scalefactor_gr_ps<<endl;
    

    THStack *stack = new THStack("stack", "");
    stack->Add(h_W2_gr_anticut);
    stack->Add(h_W2_ps_anticut);
    stack->Add(h_W2_gr_ps_anticut);
    stack->Add(h_W2_gr_ps_anticut_scaled);
    //stack->Add(h_W2_allclus_anticut);
    //stack->Add(h_W2_allclus_anticut_scaled);

    c1->cd(1);

    stack->Draw("nostack hist");
    

    //Add a legend
    auto legend = new TLegend();
    legend->SetTextSize(0.035);
    legend->SetHeader("W2 with sideband cuts");
    legend->AddEntry(h_W2_gr_anticut,"W2 w/ cut gr_ToT_sum < 20","l");
    legend->AddEntry(h_W2_ps_anticut,"W2 w/ cut ps.e < 0.2","l");
    legend->AddEntry(h_W2_gr_ps_anticut,"W2 w/ both cuts","l");
    legend->AddEntry(h_W2_gr_ps_anticut_scaled,"W2 w/ both cuts scaled","l");
    //legend->AddEntry(h_W2_allclus_anticut,"W2 allclus branch cut","l");
    // legend ->AddEntry(h_W2_allclus_anticut_scaled,"W2 allclus scaled","l");
    legend->Draw();


    THStack *stack2 = new THStack("stack2", "");
    stack2->Add(h_W2_gr_cut);
    stack2->Add(h_W2_ps_cut);
    stack2->Add(h_W2_gr_ps_cut);
    stack2->Add(h_W2);
    stack2->Add(h_W2_elastic);

    c1->cd(2);

    stack2->Draw("nostack");

    auto legend2 = new TLegend();
    legend2->SetTextSize(0.035);
    legend2->SetHeader("W2 with cuts");
    legend2->AddEntry(h_W2,"W2 no cuts","l");
    legend2->AddEntry(h_W2_gr_cut,"W2 w/ cut gr_ToT_sum > 20","l");
    legend2->AddEntry(h_W2_ps_cut,"W2 w/ cut ps.e > 0.2","l");
    legend2->AddEntry(h_W2_gr_ps_cut,"W2 w/ both cuts","l");
    legend2->AddEntry(h_W2_allclus_cut,"W2 w/ allclus branch cut","l");
    legend2->AddEntry(h_W2_elastic, "W2 with elastic cut","l");
    legend2->Draw();

    c1->SetLogy();
    c1->Update();
    
    // TH1D *h4 =  (TH1D*)(h1 ->Clone("h4"));
    // h4 ->Add(h3_scaled,-1);
    // h4 ->SetLineWidth(3);
    // h4 ->SetLineColor(30);
    // h4 ->Draw("hist");

   
    //  auto legend2 = new TLegend(0.43,0.7,0.89,0.89);
    //  legend2->SetTextSize(0.035);
    //  legend2->SetHeader("subtracted");
    //  legend2->AddEntry(h4,"gr cut - background","l");
    //  // legend2->Draw();


     // TH1D *h5 = (TH1D*)f1->Get("h_W2");
     //  h5 ->SetLineWidth(2);
     //  h5 ->SetLineColor(kBlue - 6);
     //  TH1D *h3_scaled_big =  (TH1D*)(h3 ->Clone("h3_scaled_big"));
     //  h3_scaled_big -> Scale(14);
     //  h3_scaled_big ->SetLineWidth(2);
     //  h3_scaled_big ->SetLineColor(kViolet - 6);

     //  c1 ->cd(3);
     //  h5->Draw("hist");
     //  h3_scaled_big ->Draw("hist same");

     //  auto legend3 = new TLegend(0.43,0.7,0.89,0.89);
     // legend3->SetTextSize(0.035);
     // legend3->SetHeader("W2");
     // legend3->AddEntry(h5,"W2","l");
     //  legend3->AddEntry(h3_scaled_big,"scaled background","l");
     // legend3->Draw();


     // c1 ->cd(4);      
     // TH1D *h6 =  (TH1D*)(h5 ->Clone("h6"));
     // h6 ->Add(h3_scaled_big,-1);
     // h6 ->SetLineWidth(3);
     // h6 ->SetLineColor(30);
     // h6 ->Draw("hist");
     
     //  auto legend4= new TLegend(0.43,0.7,0.89,0.89);
     // legend4->SetTextSize(0.035);
     // legend4->SetHeader("subtracted");
     // legend4->AddEntry(h4,"W2 - background","l");
     // // legend4->Draw();


    
    //}

  // c1->SaveAs("/u/home/sseeds/Plots/ecalE_allkine.png");

}
