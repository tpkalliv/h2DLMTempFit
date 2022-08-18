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


Double_t Chi2(TH1D *hY_a, TF1 *fFit, Double_t *err);

// Initializations and constants
TF1* fFit_best;
const int numbOfFVar = 100; // Number of F values
Double_t factorF[numbOfFVar];
double F_min = 0;
double F_max = 3;
TString errNames[] = {"fit_G_err","fit_V1_err","fit_V2_err ","fit_V3_err ","fit_V4_err","fit_V5_err", "hist_err"};
Double_t err[sizeof(errNames)];
TString paramNames[] = {"G const", "v11", "v22", "v33", "v44", "v55", "F"};
Int_t NH = 5;

void h2dLMTempFit() {

	// F factor values
	Double_t stepsize = (F_max-F_min)/(double) 100;
	for (int i = 0; i <= numbOfFVar; i++) factorF[i] += (i*stepsize);


	Double_t chi2_best;
 	Double_t factorF_best;
 	Int_t indexVal;
 	TH1D* hY_a[numbOfFVar];
	TF1 *fitvn_s[NH];
	TH1D* hY_a_G;
	Double_t vn[NH];
	Double_t vnError[NH];
	Double_t params[sizeof(paramNames)];
	TH1D* Y_periph;

 	// Loading data
	TFile *fIn = new TFile ("input/fout_long_range_correlation.root", "read");
	
	TH1D* hY; 
	hY = (TH1D*) fIn->Get("hDphiHM_1");

	TH1D* hY_MB;
	hY_MB = (TH1D*) fIn->Get("hDphiLM_1");



	//	FIT FUNCTION FOR PARAMETERS (From v11 to v55)
 	string cosine = "[0]*(1";
	for (int i = 1; i <= NH; i++) {
		ostringstream app;
		app << "+2*[" << i << "]*TMath::Cos(" << i << "*(x-[" << i + NH << "]))"; 
		string append = app.str();
		cosine = cosine + append;
	}
	cosine = cosine + ")";
	cout << cosine << endl;
	const char* cos = cosine.c_str();

	TF1* fFit = new TF1("fFit", cos, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);

	fFit->SetParameter(0, 1);

	for (int i = 1; i <= NH; i++) 
	{
		fFit->SetParName(i, paramNames[i]); 
		fFit->SetParameter(i, 1.0 - (i*0.04));
	}


	//	FIT FUNCTION FOR Y($Delta\\varphi$) FIT (Only v22 and v33)
 	cosine = "[0]*(1";
	for (int i=1; i < 3; i++) 
	{
		ostringstream app;
		app << "+2*[" << i << "]*TMath::Cos(" << i+1 << "*x)"; 
		string append = app.str();
		cosine = cosine + append;
	}
	cosine = cosine + ")";
	cout << "Fit for Y($Delta\\varphi$):\n" << cosine << endl;
	const char* fcos = cosine.c_str();

	TF1* fFity = new TF1("fFity", fcos, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);


	// PARAMETER EXTRACTION
 	for (int j = 0; j < numbOfFVar; j++) 
 	{

 		hY_a[j] = (TH1D*) hY->Clone(); 
 	
 		hY_a[j]->Add(hY_MB, -factorF[j]);
 	
 		hY_a[j]->Fit("fFit", "", "", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);

 		Double_t min_val = Chi2(hY_a[j], fFit, err); // CHI-SQUARE TEST

 		if (j == 0) chi2_best = min_val;
 		if (min_val < chi2_best) 
 		{
 			chi2_best = min_val;
 			indexVal = j;

 			params[0] = fFit->GetParameter(0);

			for (int l = 1; l <= NH; l++) params[l] = TMath::Power(fFit->GetParameter(l), 2);
			params[6] = factorF[j]; // Saving F value 
			fFit_best = (TF1*)fFit->Clone();
			hY_a[j]->Write();
 		}	
 	}

 	// Y($Delta\\varphi$) HIST FOR FITTING
 	Y_periph = (TH1D*) hY_MB->Clone();
 	fFity->SetParameter(0, params[0]);
 	for (int i = 1; i <= 2; i++) 
 	{
		fFity->SetParameter(i, fFit->GetParameter(i+1)); // Feeding fFity with initial values for Eval()
	}
 	for (int k = 1; k <= hY_MB->GetNbinsX(); k++) 
 	{
 		Double_t val = hY_MB->GetBinContent(k); // Taking k'th bin value
 		Double_t x = hY_MB->GetXaxis()->GetBinCenter(k);
 		Double_t paramVal = fFity->Eval(x); // Taking fit value
 		Double_t tot = (params[6]*val) + paramVal; // Adding all up
 		Y_periph->SetBinContent(k, tot);
 	}
 	hY->Fit("fFit");


 	// F*Y_LM + G DISTRIBUTION
    hY_a_G = (TH1D*) hY_MB->Clone(); 
	for (int k = 1; k <= hY_a_G->GetNbinsX(); k++) {
			double value = hY_a_G->GetBinContent(k);
			value = value*params[6];
			hY_a_G->SetBinContent(k, value);
	}
	for (int z = 1; z <= hY_a_G->GetNbinsX(); z++) {

		Double_t binV = hY_a_G->GetBinContent(z);
		hY_a_G->SetBinContent(z, binV + params[0]);
	}
	hY_a_G->SetMarkerStyle(24);


	// SAVINGS (Signal, Fit, F*Y_LM+G)
	TFile *fOut = new TFile ("output/h2dCorrFit.root", "recreate");
	hY->Write("hDphiHM"); // SIGNAL
	fFit_best->Write("fFit_best"); 
	hY_a_G->Write("hY_a_G"); // F*Y_LM+G
	fFit->Write("FIT"); // FIT


	// PRODUCING V2 AND V3 HARMONICS AND SAVING 
	Double_t Y_LM_min = hY_MB->GetMinimum(0);
	Double_t ScaleFYmin = params[6]*Y_LM_min;
	for (Int_t n=0; n<NH; n++)
	{
		TString formula = Form("[0]*(1 + 2*[1]*TMath::Cos(%d*x)) + [2]",n+1);							
		fitvn_s[n]= new TF1(Form("fit_s_v%d", n+1),formula, -TMath::Pi()/2.0, 3.0/2.0*TMath::Pi());
		vn[n] = fFit->GetParameter(1);																
		fitvn_s[n]->SetParameter(1, vn[n]);
		fitvn_s[n]->SetParameter(0, params[0]);
		fitvn_s[n]->SetParameter(2, ScaleFYmin);
		fitvn_s[n]->Write();
	}
	

 	// OUTPUTS
 	cout << "\n\n" << "Lowest Chi2: " << chi2_best << "\n" << endl;
 	cout << "PARAMETERS \n" << endl; 
 	for (int j = 0; j < 7; j++) cout << paramNames[j] << ": " << params[j] << "\n" << endl;
	cout << "Index: " << indexVal << "\n\n" << endl;
	cout << "fFity function is " << fcos << endl;
	cout << "fFity v2 " << fFity->GetParameter(1) << endl;
	cout << "fFity v3 " << fFity->GetParameter(2) << endl;
} // PROGRAM ENDS HERE










/*	
	CHI2 TESTING

	Parameters: hY' , fFit -> "Yield and fit function" 
	Returns: Double_t -> "Chi2 statistic value"  

*/
 Double_t Chi2(TH1D *hY_a, TF1 *fFit, Double_t *err)
{
	Double_t chi2 = 0.0;

	for (int i = 1; i <= hY_a->GetNbinsX(); i++) 
	{
		Double_t bincent = hY_a->GetXaxis()->GetBinCenter(i); // x-value for bin center
		Double_t obs = hY_a->GetBinContent(i); // bin value
		Double_t exp = fFit->Eval(bincent); // fit value for each bin center

		
		Double_t total_err;

		// Calculating errors 
		for (int j = 0; j < sizeof(err); j++) 
		{
			if (j != sizeof(err)-1) // Last err, which is histo
			{
				err[j] = TMath::Power(hY_a->GetBinError(i), 2);
				total_err = total_err + err[j];
			} else 
			{
				err[j] = TMath::Power(fFit->GetParError(j), 2);
				total_err = total_err + err[j];
			} 
				
		}

		Double_t val = obs - exp;
		Double_t chi2_temp = 0.0;

		if (total_err != 0) // Exclude zero denominator
		{
			chi2_temp = TMath::Power(val, 2) / total_err;
			chi2 = chi2 + chi2_temp;
		}
	}

	return chi2;
}


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