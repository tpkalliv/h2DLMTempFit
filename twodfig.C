#include "TFile.h"

void twodfig() {

	TFile* fIn = new TFile ("input/fout_corr_pp13TeV.root", "read");

	TH2D* hHM;
	TH2D* hLM;

	hHM = (TH2D*) fIn->Get("hC_6_0_0_2_11"); 

	hLM = (TH2D*) fIn->Get("hC_0_0_0_4_11"); 
	
	int etalHM = hHM->GetYaxis()->FindBin(1.6);
	int etahHM = hHM->GetYaxis()->FindBin(1.8);

	TH1D *hDphiHM;
	TH1D *hDphiLM;

	hDphiHM = (TH1D*) hHM->ProjectionY("hDphiHM", etalHM, etahHM);

	hDphiLM = (TH1D*) hLM->ProjectionY("hDphiLM", etalHM, etahHM);

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

   

}