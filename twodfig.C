void twodfig() {

	TFile* fIn = new TFile ("input/fout_corr_pp13TeV.root", "read");

	TH2D* hHM;
	TH2D* hLM;

	hHM = (TH2D*) fIn->Get("hC_6_0_0_2_11"); 

	hLM = (TH2D*) fIn->Get("hC_0_0_0_4_11"); 
	
	int etalHM = hHM->GetXaxis()->FindBin(1.6);
	int etahHM = hHM->GetXaxis()->FindBin(1.8);
	int etalHMn = hHM->GetXaxis()->FindBin(-1.8);
	int etahHMn = hHM->GetXaxis()->FindBin(-1.6);


	TH1D *hDphiHM;
	TH1D *hDphiLM;
	TH1D *hDphiHMn;
	TH1D *hDphiLMn;
	TH1D *hDphiHMp;
	TH1D *hDphiLMp;

	hDphiHMn = (TH1D*) hHM->ProjectionY("hDphiHM", etalHMn, etahHMn);

	hDphiLMn = (TH1D*) hLM->ProjectionY("hDphiLM", etalHMn, etahHMn);

	hDphiHM = (TH1D*) hHM->ProjectionY("hDphiHM", etalHM, etahHM);

	hDphiLM = (TH1D*) hLM->ProjectionY("hDphiLM", etalHM, etahHM);

	hDphiHM->Add(hDphiHMn, 1.);
	hDphiLM->Add(hDphiLMn, 1.);

	TCanvas *c1 = new TCanvas ("hHM", "hHM", 1);
	c1->Divide(2, 1);

    c1->cd(1);
    hHM->Draw("colz");
    
    c1->cd(2);
    hLM->Draw("colz");
  	
  	TCanvas *c2 = new TCanvas ("hHMP", "hHMP", 1);
	c2->Divide(2, 1);  
    
    c2->cd(1);hDphiHM->Draw();
    c2->cd(2);hDphiLM->Draw();

    TFile* fOut = new TFile ("input/fout_long_range_correlation.root", "recreate");

    hDphiLM->Write("hDphiLM_1");
    hDphiHM->Write("hDphiHM_1");
    double YM_min = hDphiLM->GetBinContent(hDphiLM->GetBinCenter(0.));

}