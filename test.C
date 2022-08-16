




void test() {

	TFile *fIn = new TFile ("input/fout_long_range_correlation.root", "read");

	TH1D* hist;

	hist = (TH1D*) fIn->Get("hDphiLM_1");


	if (hist==NULL) cout << "Histogram is NULL" << endl;

	hist->Print();

	fIn->Close();

}