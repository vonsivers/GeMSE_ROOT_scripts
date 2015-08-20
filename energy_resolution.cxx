#include "macros.h"

TString fFileName_spectrum;
std::vector<double> ffitrange_low;
std::vector<double> ffitrange_high;
std::vector<double> famp_st;
std::vector<double> fmean_st;
std::vector<double> fsigma_st;
std::vector<double> fconst_st;
std::vector<double> fslope_st;


// ----------------------------------------------------
// read parameters from user specified file
// ----------------------------------------------------

int read_parameters(TString FileName) {
    
    double fitrange_low, fitrange_high, amp_st, mean_st, sigma_st, const_st, slope_st;
    
    ifstream File;
    File.open(FileName);
    
    if (!File.is_open()) {
        std::cout << "##### ERROR: could not open " << FileName << std::endl;
        return 1;
    }
    
    std::string headerline;
    
    getline(File, headerline);
    File >> fFileName_spectrum;
    getline(File, headerline);
    getline(File, headerline);
    while (true)
    {
        File >> fitrange_low >> fitrange_high >> amp_st >> mean_st >> sigma_st >> const_st >> slope_st;
        ffitrange_low.push_back(fitrange_low);
        ffitrange_high.push_back(fitrange_high);
        famp_st.push_back(amp_st);
        fmean_st.push_back(mean_st);
        fsigma_st.push_back(sigma_st);
        fconst_st.push_back(const_st);
        fslope_st.push_back(slope_st);
        if( File.eof() ) break;
        
    }
    
    File.close();
    
    // print information
    std::cout << "######################################" << std::endl;
    std::cout << "#### Reading "<< FileName << " ..." << std::endl;
    std::cout << "spectrum file name: " << fFileName_spectrum << std::endl;
    std::cout << "fitrange low \t fitrange high \t counts \t mean \t sigma \t const. \t slope" << std::endl;
    for (int i=0; i<ffitrange_low.size(); ++i)
    {
        std::cout << ffitrange_low[i] << "\t" << ffitrange_high[i] << "\t" << famp_st[i] << "\t" << fmean_st[i] << "\t" << fsigma_st[i]<< "\t"  << fconst_st[i]<< "\t"  << fslope_st[i]<< "\t"  << std::endl;
    }
    std::cout << "######################################" << std::endl;
    
    return 0;
}

// ----------------------------------------------------
// determine resolution
// ----------------------------------------------------

int main(int argc, char *argv[]) {
    
    // argc should be 2 for correct execution
    if ( argc != 2 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <parameters_file.txt>" << std::endl;
        return 1;
    }
    
    TString FileName_parameters = argv[1];
    
    // read parameters from file
    if (read_parameters(FileName_parameters)) {
        return 1;
    }
    
    // get spectrum
    TH1D* hist = getspectrum(fFileName_spectrum);
    
    if (hist==0) {
        return 1;
    }
    
    // draw spectrum
    TCanvas* c1 = new TCanvas("c1");
    hist->Draw();
    
    TF1* fit = new TF1();
    std::vector<double> peakpos, peakpos_err, sigma, sigma_err;
    
    int npeaks = ffitrange_low.size();

    // loop over all peaks
    for (int i=0; i<npeaks; ++i) {
        
        // fit peak with Gauss+Pol1
        std::cout << "###################################################" << std::endl;
        std::cout << "Fitting line " << i+1 << " in range " << ffitrange_low[i] << " - " << ffitrange_high[i] << std::endl;
        fit = FitGaussPol1(hist,famp_st[i],fmean_st[i], fsigma_st[i], fconst_st[i], fslope_st[i], ffitrange_low[i], ffitrange_high[i]);
        if (fit==0) {
            return 1;
        }
        peakpos.push_back(fit->GetParameter(1));
        peakpos_err.push_back(fit->GetParError(1));
        sigma.push_back(fit->GetParameter(2));
        sigma_err.push_back(fit->GetParError(2));

    }
    
    // save calibration fits
    c1->SaveAs(fFileName_spectrum+"_resolution_fits.root");
    c1->SaveAs(fFileName_spectrum+"_resolution_fits.pdf");
    
    // plot peak position vs. energy
    TGraphErrors* graph = new TGraphErrors(npeaks, &peakpos.at(0), &sigma.at(0), &peakpos_err.at(0), &sigma_err.at(0) );
    
    // fit with sqrt
    std::cout << "###################################################" << std::endl;
    std::cout << "Fitting resolution curve ..." << std::endl;
    TF1* fit_res = FitSqrt(graph,0.1,0.,0.);
    if (fit_res==0) {
        return 1;
    }
    
    TCanvas* c2 = new TCanvas("c2");
    graph->SetTitle("Energy Resolution");
    graph->SetMarkerStyle(2);
    graph->GetXaxis()->SetTitle("Energy (keV)");
    graph->GetYaxis()->SetTitle("Sigma (keV)");
    graph->Draw("ap");
    
    c2->SaveAs(fFileName_spectrum+"_resolution_function.root");
    c2->SaveAs(fFileName_spectrum+"_resolution_function.pdf");
    
    return 0;
}


