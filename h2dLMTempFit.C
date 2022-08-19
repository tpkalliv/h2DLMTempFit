#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"


/*

	How this program works:

	Program loads two data sets and finds the best fit for them using Chi2.
	Outputs Chi2 statistical value and also parameters from the best fit.
	
	Background:

	G(fourier) is actually HM flow. Normalization constant is there to make it the actual HM flow.
	So that when we add Y_LM to it, with the assumption that it has only non-flow, we can add to the amount of yield in HM
	when we modify Y_LM non-flow yield with F factor to take into account the relative difference in the two events for non-flow.
	This way it is possible to extract only the HM flow part by subtracting F*Y_LM (HM non-flow) from the HM yield. 
*/


Double_t Chi2(TH1D *hY_a, TF1 *fFit);

// Initializations and constants

const int numbOfFVar = 100; // Number of F values
Double_t factorF[numbOfFVar];
double F_min = 0.9;
double F_max = 2;
TString errNames[] = {"fit_G_err","fit_V2_err ","fit_V3_err "};
TString paramNames[] = {"G", "v22", "v33"};
Int_t NH = 2; // 2-3

void h2dLMTempFit() {

	// F factor values
	Double_t stepsize = (F_max-F_min)/(double) 100;
	for (int i = 0; i <= numbOfFVar; i++) factorF[i] = F_min + (i*stepsize);
	cout << Form("%.2f<F<%.2f, NF= %d, step=%.2f\n",F_min,F_max,numbOfFVar,stepsize) << endl;

 	TH1D* hY_a;
 	Double_t chi2all[numbOfFVar];
	TF1 *fitvn_s[NH];
	TH1D* hY_a_G;
	Double_t vn[NH];
	Double_t vnError[NH];
	Double_t params[sizeof(paramNames)];
	TH1D* hFitTotal;
 	// Loading data
	TFile *fIn = new TFile ("input/fout_long_range_correlation.root", "read");
	
	TH1D* hY = (TH1D*) fIn->Get("hDphiHM_1");
	TH1D* hY_MB = (TH1D*) fIn->Get("hDphiLM_1");
	TH1D *hchiq2 = new TH1D("hchiq2","chi2 of fits",100,0,5);

	//	FIT FUNCTION FOR PARAMETERS (From v11 to v55)
 	string cosine = "[0]*(1";
	for (int i = 0; i < NH; i++) {
		ostringstream app;
		app << "+2*[" << i+1 << "]*TMath::Cos(" << i+2 << "*x)"; // Coefficients are Vn(delta)phi = Vn^2
		string append = app.str();
		cosine = cosine + append;
	}
	cosine = cosine + ")";
	cout << cosine << endl;
	const char* cos = cosine.c_str();
  

	TF1* fFit[sizeof(paramNames)]; 
	for (int j = 0; j < numbOfFVar; j++){
		fFit[j]= new TF1(Form("fFit%d",j), cos, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
		fFit[j]->SetParameter(0, 1);
		for (int i = 0; i < NH; i++) fFit[j]->SetParName(i+1, paramNames[i]); 
		for (int i = 0; i < NH; i++) fFit[j]->SetParameter(i+1, TMath::Power(1.0 - (i*0.06),2)); // Initial Vn values are Vn(delta)phi = Vn^2
	}
	
	// Fit with a given F
	Double_t min_val=-999;
	// PARAMETER EXTRACTION
 	for (int j = 0; j < numbOfFVar; j++) 
 	{
 		hY_a = (TH1D*) hY->Clone(); 
 		hY_a->Add(hY_MB, -factorF[j]);
 		hY_a->Fit(fFit[j], "", "", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
    	chi2all[j] = fFit[j]->GetChisquare()/fFit[j]->GetNDF();//Chi2(hY_a, fFit); // CHI-SQUARE TEST
 		hchiq2->Fill(chi2all[j]);		
 		cout << Form("F[%d] = %.2f, chi2=%0.2f",j,factorF[j],chi2all[j]) << endl;
 	}
 	Double_t chiq_min = chi2all[0];
 	int index = -1;
    // search num in inputArray from index 0 to elementCount-1 
    for (int j = 1; j < numbOfFVar; j++) {
	    if(chi2all[j] < chiq_min){
	        chiq_min = chi2all[j];
            index = j;
        }
    }
    cout << Form("chi2_min=%0.4f, index=%d/%d\n",chi2all[index], index,numbOfFVar);
  
	params[0] = fFit[index]->GetParameter(0);
	for (int l = 1; l < NH; l++) params[l] = fFit[index]->GetParameter(l); // Saving Vn^2 values

    cout << "Saving results in to the ROOT file" << endl;
 	hFitTotal = (TH1D*) hY_MB->Clone();
 	for (int k = 1; k <= hY_MB->GetNbinsX(); k++) 
 	{
 		Double_t ylm = hY_MB->GetBinContent(k); // Taking k'th bin value
 		Double_t x = hY_MB->GetXaxis()->GetBinCenter(k);
 		Double_t tot =  fFit[index]->Eval(x) + (factorF[index]*ylm); // Adding all up
 		hFitTotal->SetBinContent(k, tot);
 	}

	// F*Y_LM + G DISTRIBUTION
    hY_a_G = (TH1D*) hY_MB->Clone(); 
	for (int k = 1; k <= hY_a_G->GetNbinsX(); k++) {
			double value = hY_a_G->GetBinContent(k);
			value = value*factorF[index] + params[0];
			hY_a_G->SetBinContent(k, value);
	}
	hY_a_G->SetMarkerStyle(24);


	// SAVINGS (Signal, Fit, F*Y_LM+G)
	TFile *fOut = new TFile ("output/h2dCorrFit.root", "recreate");
	hY->Write("hDphiHM"); // SIGNAL
	fFit[index]->Write("fFit_best"); 
	hY_a_G->Write("hY_a_G"); // F*Y_LM+G
	hFitTotal->Write("hFitTotal");
	hchiq2->Write();

	// PRODUCING V2 AND V3 HARMONICS AND SAVING 
	Double_t Y_LM_min = hY_MB->GetMinimum(0);
	Double_t ScaleFYmin = factorF[index]*Y_LM_min;
	for (Int_t n=0; n<NH; n++)
	{
		TString formula = Form("[0]*(1 + 2*[1]*TMath::Cos(%d*x)) + [2]",n+2);					
		fitvn_s[n]= new TF1(Form("fit_s_v%d", n+2),formula, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
		vn[n] = fFit[index]->GetParameter(1);																
		fitvn_s[n]->SetParameter(1, vn[n]);
		fitvn_s[n]->SetParameter(0, params[0]);
		fitvn_s[n]->SetParameter(2, ScaleFYmin);
		fitvn_s[n]->Write();
	}

	TCanvas *c2 = new TCanvas ("hFit", "hHM", 1);
	hY->SetMinimum(hY->GetMinimum()*0.98);
	hY->Draw();
	hY_a_G->SetLineColor(2);
	hY_a_G->Draw("same");
	fitvn_s[0]->Draw("same");
	fitvn_s[1]->Draw("same");
	hFitTotal->SetLineColor(4);
	hFitTotal->Draw("same");
	
} // PROGRAM ENDS HERE

/* 

STEPS IN THE ALGORITHM:

1. 	Loads two input histos
2. 	Creates G(fourier) fit function
3. 	Gives fit some initial values
4. 	Multiplies hY_MB histo with F_i value
5. 	Substracts hY_MB from hY histo to create hY' histo
6. 	Fits G(fourier) to hY'
7. 	Calculates Chi2 value
8. 	Compares Chi2 value to the previous best Chi2 value
9. 	Saves histos and fits to .root file
10.	Outputs best chi2 value and associated parameters on screen

*/