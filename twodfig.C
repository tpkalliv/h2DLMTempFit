void twodfig() {

	TFile* fIn = new TFile ("input/fout_corr_pp13TeV.root", "read");

	TH2D* hHM;
	TH2D* hLM;

	hHM = (TH2D*) fIn->Get("hC_6_0_0_2_11"); 

	hLM = (TH2D*) fIn->Get("hC_0_0_0_4_11"); 

	double etalow=1.6;
	double etahigh=1.8;
	
	TH1D *hDphiHM; // default negative side only
	TH1D *hDphiHMp;
	TH1D *hDphiLM;
	TH1D *hDphiLMp;

	hDphiHM = (TH1D*) hHM->ProjectionY("hDphiHM", hHM->GetXaxis()->FindBin(-etahigh), hHM->GetXaxis()->FindBin(-etalow));
	hDphiHMp  = (TH1D*) hHM->ProjectionY("hDphiHMp",  hHM->GetXaxis()->FindBin(etalow), hHM->GetXaxis()->FindBin(etahigh));
	hDphiHM->Add(hDphiHMp, 1.);

	hDphiLM = (TH1D*) hLM->ProjectionY("hDphiLM", hLM->GetXaxis()->FindBin(-etahigh), hLM->GetXaxis()->FindBin(-etalow));
	hDphiLMp  = (TH1D*) hLM->ProjectionY("hDphiLMp",  hLM->GetXaxis()->FindBin(etalow), hLM->GetXaxis()->FindBin(etahigh));

	hDphiLM->Add(hDphiLMp, 1.);

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

	TCanvas *c3 = new TCanvas ("hHMPp", "hHMPp", 1);
	c3->Divide(2, 1);  

	c3->cd(1);hDphiHMp->Draw();
	c3->cd(2);hDphiLMp->Draw();

	c1->SaveAs("figs/QA_2D.png");
	c2->SaveAs("figs/QA_alllongrange.png");
	c1->SaveAs("figs/QA_plusetarange.png");

	TFile* fOut = new TFile ("input/fout_long_range_correlation.root", "recreate");

	hDphiLM->Write("hDphiLM_1");
	hDphiHM->Write("hDphiHM_1");
	double YM_min = hDphiLM->GetBinContent(hDphiLM->GetBinCenter(0.));

}
