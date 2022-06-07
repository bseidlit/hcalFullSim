#include "/sphenix/u/bseidlitz/plotstyle/AtlasStyle.C"
#include "/sphenix/u/bseidlitz/plotstyle/AtlasUtils.C"

void PlotResults(){

  SetAtlasStyle();

  string figPath = "figures/";

  TFile* fin = new TFile("rootFiles/trigAna2.root");

  TH2F* h2_twr_template = new TH2F("h2_twr_template","",24,-1.1,1.1,64,0,2*TMath::Pi());

  TH1F* h_truth_mu_E_surf = (TH1F*) fin->Get("h_truth_mu_E_surf");

  float N_mu = h_truth_mu_E_surf->Integral(h_truth_mu_E_surf->FindBin(1.0),h_truth_mu_E_surf->GetNbinsX()+1);
  float N_O_m2st = N_mu/4.*1.5/TMath::Pi()/0.66;
  float time = N_O_m2st/70.;
  cout << "time: " << time << "[s]" << endl;

  // divide by area and bin width
  h_truth_mu_E_surf->Scale(1./4,"width");
  // divid st and account for wide acceptance 
  h_truth_mu_E_surf->Scale(1.5/TMath::Pi()/0.66);
  //  scale by time
  h_truth_mu_E_surf->Scale(1./time);

  TH1F* h_perp_twr   = (TH1F*) fin->Get("h_perp_twr");
  TH1F* h_active_twr = (TH1F*) fin->Get("h_active_twr");

  float total_perp_events   = h_perp_twr->Integral(2,20);
  float total_active_events = h_active_twr->Integral(2,20);
 
  TCanvas* c11 = new TCanvas("c11","c11",600,600);
  h_truth_mu_E_surf->Draw("ex0");
  h_truth_mu_E_surf->SetXTitle("E_#mu [GeV]");
  h_truth_mu_E_surf->SetYTitle("dN/(dE dst dm^{2}) [1/(m^{2} st GeV s)]");

  gPad->SetLogy();
  gPad->SetLogx();

  c11->SaveAs(Form("%s/muon_E.pdf",figPath.c_str()));
  delete c11;

  TH2F* h_fit_counts_perp = (TH2F*) fin->Get("h_fit_counts_perp");

  float perp_countsOtime = h_fit_counts_perp->Integral()/time;

  h_fit_counts_perp->Scale(1./time);

  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  h_fit_counts_perp->Draw("COLZ");
  h_fit_counts_perp->SetXTitle("#eta^{twr}");
  h_fit_counts_perp->SetYTitle("#phi^{twr}");
  h_fit_counts_perp->SetZTitle("counts/s");
  h_fit_counts_perp->GetZaxis()->SetRangeUser(0,3.3);
  gPad->SetRightMargin(0.2);
 
  myText(0.22,0.9,1,Form("events/sec %0.1f",total_perp_events/time));
 
  c2->SaveAs(Form("%s/hcalOut_eta_phi_perpCut.pdf",figPath.c_str()));
  delete c2;


  const int nt_phi = 64;
  const int nt_eta = 24;

  TH2F* h_fit_counts_para = new TH2F("h_fit_counts_perp","",24,-1.1,1.1,64,0,2*TMath::Pi());
  TH2F* h2_counts_E = new TH2F("h2_counts_E","",24,-1.1,1.1,64,0,2*TMath::Pi());
  

   TH1F* tower_eta_phi_E[nt_eta][nt_phi];
   TH1F* tower_eta_phi_E_filt[nt_eta][nt_phi];
   TH1F* tower_eta_phi_E_perp[nt_eta][nt_phi];
   TH1F* tower_eta_phi_E_para[nt_eta][nt_phi];
   for(int i = 0; i<nt_eta;i++){
     for(int j = 0; j<nt_phi; j++){
        tower_eta_phi_E[i][j]      = (TH1F*) fin->Get(Form("tower_eta_phi_E_%d_%d",i,j));
        tower_eta_phi_E_filt[i][j] = (TH1F*) fin->Get(Form("tower_eta_phi_E_filt_%d_%d",i,j));
        tower_eta_phi_E_para[i][j] = (TH1F*) fin->Get(Form("tower_eta_phi_E_para_%d_%d",i,j));
        tower_eta_phi_E_perp[i][j] = (TH1F*) fin->Get(Form("tower_eta_phi_E_perp_%d_%d",i,j));
       h_fit_counts_para->SetBinContent(i+1,j+1,tower_eta_phi_E_para[i][j]->GetEntries());
       h2_counts_E->SetBinContent(i+1,j+1,tower_eta_phi_E[i][j]->GetEntries());
     }
   }


  ////////////////////
  // Fit MIP peak
   TF1* f_MIP = new TF1("f_MIP","[3] *pow( (x-[1])/[2] , (([0]-[1])/[2]) ) *TMath::Exp(-(x-[1])/[2])/ (TMath::Gamma((([0]-[1])/[2])+1)*[2])",20,200);
   //TF1* f_MIP = new TF1("f_MIP","TMath::GammaDist(x,[1],[2],[3])+[0]",6,200);
   //TF1* f_MIP = new TF1("f_MIP","gaus",6,200);
   double par[4] = {40,24,10,500};
   f_MIP->SetParLimits(0,28,55);
   f_MIP->SetParLimits(1,9,30);
   f_MIP->SetParLimits(2,5,50);
   f_MIP->SetParameters(par);
  f_MIP->SetLineColor(kRed); 

   tower_eta_phi_E_perp[11][16]->Fit(f_MIP,"","r",25,85);


  TCanvas* c6 = new TCanvas("c6","c6",600,600);
  tower_eta_phi_E_perp[11][16]->Draw();
  //f_MIP->Draw("same");

   c6->SaveAs("figures/sample_perp_dist.pdf");
   delete c6;

   TH1F* h_precision = new TH1F("h_precision","",100,0,0.1);
   TH2F* h2_perp_peak = (TH2F*) h2_twr_template->Clone("h2_perp_peak");

   float best_precision = 1000;
   float worst_precision = 0;
   for(int i = 0; i<nt_eta;i++){
     for(int j = 0; j<nt_phi; j++){
        int res = tower_eta_phi_E_perp[i][j]->Fit("f_MIP","","r",26,80);
        float precision = f_MIP->GetParError(0)/f_MIP->GetParameter(0);
        if (precision < best_precision) best_precision = precision;
        if (precision > worst_precision) worst_precision = precision;
        h_precision->Fill(precision);
        if (res != 0) std::cout << std::endl << "precision " << precision << "   " << tower_eta_phi_E_perp[i][j]->GetEntries();
        h2_perp_peak->SetBinContent(i+1,j+1,f_MIP->GetParameter(0));
     }
  }

  TGraph *gr_best_precision  = new TGraph();
  TGraph *gr_worst_precision = new TGraph();
  gr_best_precision ->SetName("gr_best_precision"); 
  gr_worst_precision->SetName("gr_worst_precision"); 

  gr_best_precision  ->SetPoint(0,time,best_precision);
  gr_worst_precision ->SetPoint(0,time,worst_precision);

  TF1 *f_oneOn = new TF1("f_oneOn","[0]/TMath::Sqrt(x)",0.1,time*2);
  TF1 *f_oneOn_2 = new TF1("f_oneOn_2","[0]/TMath::Sqrt(x)",0.1,time*2);
  gr_best_precision->Fit(f_oneOn);
  gr_worst_precision->Fit(f_oneOn_2);

  TCanvas* c8 = new TCanvas("c8","c8",600,600);
  TH1F* h_frame = new TH1F("h_frame","",1,0,time*2);
  h_frame->GetYaxis()->SetRangeUser(0.001,1);
  h_frame->GetXaxis()->SetRangeUser(0,time*2);
  h_frame->GetXaxis()->SetTitle("live time [s]");
  h_frame->GetYaxis()->SetTitle("precision (peak error/value)");
  h_frame->Draw("AXIS");

  gr_best_precision->Draw("P");
  gr_worst_precision->Draw("P");
  gr_worst_precision->SetMarkerStyle(23);
  gr_worst_precision->SetMarkerColor(kRed);
  f_oneOn->Draw("same");
  f_oneOn_2->Draw("same");
  f_oneOn_2->SetLineColor(kRed);

  myMarkerLineText(0.5,0.85,1.5,kRed,23,kRed,1,"worst twr peak precision",0.04,true);
  myMarkerLineText(0.5,0.80,1.5,kBlack,20,kBlack,1,"best twr peak precision",0.04,true);
  gPad->SetLogy();
 
  c8->SaveAs(Form("%s/error_time_perp.pdf",figPath.c_str()));
  delete c8;
  


  TCanvas* c10 = new TCanvas("c10","c10",600,600);
  h2_perp_peak->Draw("COLZ");
  h2_perp_peak->SetXTitle("#eta^{twr}");
  h2_perp_peak->SetYTitle("#phi^{twr}");
  h2_perp_peak->SetZTitle("peak [ADC]");
  h2_perp_peak->GetZaxis()->SetRangeUser(0,50);
  gPad->SetRightMargin(0.2);
  
  c10->SaveAs(Form("%s/perp_peak.pdf",figPath.c_str()));
  delete c10;


  //////////////////////
  //  raw counts
  h_fit_counts_para->Scale(1./time);

  TCanvas* c3 = new TCanvas("c3","c3",600,600);
  h_fit_counts_para->Draw("COLZ");
  h_fit_counts_para->SetXTitle("#eta^{twr}");
  h_fit_counts_para->SetYTitle("#phi^{twr}");
  h_fit_counts_para->SetZTitle("counts/s");
  gPad->SetRightMargin(0.2);

  c3->SaveAs(Form("%s/counts_para.pdf",figPath.c_str()));
  delete c3;



  float noCut_total = h2_counts_E->Integral()/time;
  h2_counts_E->Scale(1./time);

  TCanvas* c4 = new TCanvas("c4","c4",600,600);
  h2_counts_E->Draw("COLZ");
  h2_counts_E->SetXTitle("#eta^{twr}");
  h2_counts_E->SetYTitle("#phi^{twr}");
  h2_counts_E->SetZTitle("counts/s");
  h2_counts_E->GetZaxis()->SetRangeUser(0,35);
  gPad->SetRightMargin(0.2);
  
  myText(0.22,0.9,1,Form("events/sec %.1e",total_active_events/time));
  
  c4->SaveAs(Form("%s/hcalOut_eta_phi_noCut.pdf",figPath.c_str()));
  delete c4;


  TCanvas* c9 = new TCanvas("c9","c9",600,600);
  h_precision->Draw();
  h_precision->GetXaxis()->SetNdivisions(505);
  h_precision->GetXaxis()->SetTitle("peak error/vale");
  h_precision->GetYaxis()->SetTitle("Outer HCal towers");


   /////////////////////////////////
   // Inner HCal analisys
   //////////////////////////////////
   TH1F* h_active_hcalIn = (TH1F*) fin->Get("h_active_hcalIn");
   TH1F* h_perp_hcalIn = (TH1F*) fin->Get("h_perp_hcalIn");
   TH2F* h2_hcalIn_active = (TH2F*) h2_twr_template->Clone("h2_hcalIn_active");
   TH2F* h2_hcalIn_perp = (TH2F*) h2_twr_template->Clone("h2_hcalIn_perp");

   TH1F* hcalIn_eta_phi_E_perp[nt_eta][nt_phi];
   TH1F* hcalIn_eta_phi_E[nt_eta][nt_phi];

   for(int i = 0; i<nt_eta;i++){
     for(int j = 0; j<nt_phi; j++){
       hcalIn_eta_phi_E_perp[i][j] = (TH1F*) fin->Get(Form("hcalIn_eta_phi_E_perp_%d_%d",i,j));
       hcalIn_eta_phi_E     [i][j] = (TH1F*) fin->Get(Form("hcalIn_eta_phi_E_%d_%d",i,j));
       h2_hcalIn_active->SetBinContent(i+1,j+1,hcalIn_eta_phi_E     [i][j]->Integral());
       h2_hcalIn_perp  ->SetBinContent(i+1,j+1,hcalIn_eta_phi_E_perp[i][j]->Integral());
     }
   }

   h2_hcalIn_active->Scale(1.0/time);
   h2_hcalIn_perp  ->Scale(1.0/time);
   float active_hcalIn_event_rate = h_active_hcalIn->Integral(2,30)/time;
   float perp_hcalIn_event_rate   = h_perp_hcalIn->Integral(2,30)/time;

   TCanvas* c12 = new TCanvas("c12","c12",600,600);
   h2_hcalIn_active->Draw("COLZ");
   h2_hcalIn_active->SetXTitle("#eta^{twr}");
   h2_hcalIn_active->SetYTitle("#phi^{twr}");
   h2_hcalIn_active->SetZTitle("counts/s");
   h2_hcalIn_active->GetZaxis()->SetRangeUser(0,5);
   gPad->SetRightMargin(0.2);
   
   myText(0.22,0.9,1,Form("events/sec %0.1f",active_hcalIn_event_rate));

   c12->SaveAs("figures/hcalIn_eta_phi_noCuts.pdf");
   

   TCanvas* c13 = new TCanvas("c13","c13",600,600);
   h2_hcalIn_perp->Draw("COLZ");
   h2_hcalIn_perp->SetXTitle("#eta^{twr}");
   h2_hcalIn_perp->SetYTitle("#phi^{twr}");
   h2_hcalIn_perp->SetZTitle("counts/s");
   h2_hcalIn_perp->GetZaxis()->SetRangeUser(0,0.3);
   gPad->SetRightMargin(0.2);
   
  myText(0.22,0.9,1,Form("events/sec %0.1f",perp_hcalIn_event_rate));

  c13->SaveAs("figures/hcalIn_eta_phi_perpCuts.pdf");

  
}
