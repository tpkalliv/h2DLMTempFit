#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"


 // Initializing Chi2 function
 Double_t Chi2(TH1D *hY_a, TF1 *fFit, Double_t *err);


/*
	Program loads two data sets and finds the best fit for them using Chi2.
	Outputs Chi2 statistical value and also parameters from the best fit.
*/
void h2dLMTempFit() {

	// Initializings 
	Int_t numbOfFVar = 100; // Number of F values
	Double_t factorF[numbOfFVar];
	Double_t stepsize = (3-1)/(double) 100;
 	Double_t chi2_best;
 	Double_t factorF_best;
 	Int_t indexVal;
	Int_t NH = 5;
 	TH1D* hY_a[numbOfFVar];
	TF1 *fitvn[NH];
	TH1D* hY_a_G;
	Double_t vn[NH];
	Double_t vnError[NH];
	TString errNames[] = {"fit_G_err","fit_V1_err","fit_V2_err ","fit_V3_err ","fit_V4_err","fit_V5_err", "hist_err"};
	Double_t err[sizeof(errNames)];
	TString paramNames[] = {"G const", "v11", "v22", "v33", "v44", "v55", "F"};
	Double_t params[sizeof(paramNames)];

 	// Opens data 
	TFile *fIn = new TFile ("input/fout_long_range_correlation.root", "read");

	fIn->Print();

	//TH1D* chi2_hist = new TH1D ("chi2", "chi2", 1001, 0, 200);

	cout << "1" << endl;

	
	TH1D* hY; 
	hY = (TH1D*) fIn->Get("hDphiHM_1");

	cout << "2" << endl;

	TH1D* hY_MB;
	hY_MB = (TH1D*) fIn->Get("hDphiLM_1");
	
	cout << "3" << endl;

	//	Fit function
 	string cosine = "[0]*(1";
	for (int i=1; i<=NH; i++) {
		ostringstream app;
		app << "+2*[" << i << "]*TMath::Cos(" << i << "*(x-[" << i + NH << "]))"; 
		string append = app.str();
		cosine = cosine + append;
	}
	cosine = cosine + ")";
	cout << cosine << endl;
	const char* cos = cosine.c_str();

	cout << "4" << endl;

	TF1* fFit = new TF1("fFit", cos, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);

	//fFit->SetParLimits(0, -10, 10);
	fFit->SetParameter(0, 1);

	for (int i = 1; i <= NH; i++) 
	{
		fFit->SetParName(i, paramNames[i]); // Param 1: V22 and Param 2: V33 ...
		//fFit->SetParLimits(i, -1, 1);
		fFit->SetParameter(i, 1.0 - (i*0.04));
	}
	
	cout << "5" << endl;

	// F factors
 	for (int i = 0; i <= numbOfFVar; i++) 
 	{
 		factorF[i] = 1 + (i*stepsize);
 	}	

	cout << "6" << endl;

	TFile *fOut = new TFile ("output/h2dCorrFit.root", "recreate");

	// 	Multiplying, subtracting, fitting and Chi2 testing
 	for (int j = 0; j < numbOfFVar; j++) 
 	{
 		cout << "7" << endl;

 		hY_a[j] = (TH1D*) hY->Clone(); 
 		cout << "8" << endl;
 		hY_a[j]->Add(hY_MB, -factorF[j]);
 		cout << "9" << endl;
 		hY_a_G = (TH1D*) hY_a[j]->Clone(); 
 		cout << "10" << endl;
 		hY_a[j]->Fit("fFit", "", "", -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0);
	
 		Double_t min_val = Chi2(hY_a[j], fFit, err); // Chi-Square testing

 		if (j == 0) chi2_best = min_val;

 		// Saving histos and fits
 		if (min_val < chi2_best) 
 		{
 			chi2_best = min_val;
 			indexVal = j;



 			params[0] = fFit->GetParameter(0);

			for (int l = 1; l <= NH; l++) params[l] = TMath::Power(fFit->GetParameter(l), 2);
			params[7] = factorF[j]; // Saving F value 

			hY_a[j]->Write();
			fFit->Write();

			
			for (int z = 0; z < hY_a_G->GetNbinsX(); z++) {

				Double_t binV = hY_a_G->GetBinContent(z);
				hY_a_G->SetBinContent(z, binV + params[0]);
			}
			hY_a_G->Write("F*Y+G");
			
 		}	
 	}


 	// Saving harmonics
	for (Int_t n=0; n<NH; n++)
	{
		TString formula = Form("[0]*(1 + 2*[1]*TMath::Cos(%d*x))",n+1);									
		fitvn[n]= new TF1(Form("fit_v%d", n+1),formula, -TMath::Pi()/2.0, 3.0/2.0*TMath::Pi());			
		vn[n] = fFit->GetParameter(n+1);																	
		vnError[n] = fFit->GetParError(n+1);																
		fitvn[n]->SetParameter(1,vn[n]);
		fitvn[n]->SetParameter(0, params[0]);
		fitvn[n]->Write();
	}

	Double_t Y_LM_min = hY_MB->GetBinContent(hY_MB->GetBinCenter(0.));
	
	Double_t ScaleFYmin = params[7]*Y_LM_min;

	TF1* fitvn_s[NH];

	// Saving harmonics
	for (Int_t n=0; n<NH; n++)
	{
		TString formula = Form("[0]*(1 + 2*[1]*TMath::Cos(%d*x)) + [2]",n+1);									
		fitvn_s[n]= new TF1(Form("fit_s_v%d", n+1),formula, -TMath::Pi()/2.0, 3.0/2.0*TMath::Pi());																		
		fitvn_s[n]->SetParameter(1,vn[n]);
		fitvn_s[n]->SetParameter(0, params[0]);
		fitvn_s[n]->SetParameter(2,ScaleFYmin);
		fitvn_s[n]->Write();
	}
	
 	// Outputs
 	cout << "\n\n" << "Lowest Chi2: " << chi2_best << "\n" << endl;
 	cout << "PARAMETERS \n" << endl; 
 	for (int j = 0; j < 7; j++) cout << paramNames[j] << ": " << params[j] << "\n" << endl;
	cout << "Index: " << indexVal << "\n\n" << endl;

	// 	fIn->Close();
 	//fOut->Close();
} 

/*	
	Chi2 Test

	Parameters: hY' , fFit -> "Yield and fit function" 
	Returns: Double_t -> "Chi2 statistic value"  

*/
 Double_t Chi2(TH1D *hY_a, TF1 *fFit, Double_t *err)
{
	Double_t chi2 = 0.0;

	for (int i = 1; i < hY_a->GetNbinsX(); i++) 
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