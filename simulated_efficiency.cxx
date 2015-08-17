#include <TH1D.h>
#include <TGraph.h>
#include <TF1.h>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TVectorD.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>

TString fInputFolder;
std::vector<TString> fIsotopeName;
std::vector<double> fEnergy, fWidth_sig, fWidth_bck, fBR;
std::vector<int> fFlag;


// ----------------------------------------------------
// read parameters from user specified file
// ----------------------------------------------------

int read_parameters(TString FileName) {
    
    TString IsotopeName;
    double energy, BR, width_sig, width_bck;
    int flag;
    
    ifstream File;
    File.open(FileName);
    
    if (!File.is_open()) {
        std::cout << "##### ERROR: could not open " << FileName << std::endl;
        return 1;
    }
    
    std::string headerline;
    
    getline(File, headerline);
    File >> fInputFolder;
    getline(File, headerline);
    getline(File, headerline);
    while (true)
    {
        File >> IsotopeName >> energy >> BR >> width_sig >> width_bck >> flag;
        fIsotopeName.push_back(IsotopeName);
        fEnergy.push_back(energy);
        fBR.push_back(BR);
        fWidth_sig.push_back(width_sig);
        fWidth_bck.push_back(width_bck);
        fFlag.push_back(flag);

        if( File.eof() ) break;
        
    }
    
    File.close();
    
    // print information
    std::cout << "######################################" << std::endl;
    std::cout << "#### Reading "<< FileName << " ..." << std::endl;
    std::cout << "folder of simulation files: " << fInputFolder << std::endl;
    std::cout << "Isotope \t Energy \t Width signal region \t Width bkg region" << std::endl;
    for (int i=0; i<fEnergy.size(); ++i)
    {
        std::cout << fIsotopeName[i] << "\t" << fEnergy[i] << "\t" << fWidth_sig[i] << "\t" << fWidth_bck[i] << std::endl;
    }
    std::cout << "######################################" << std::endl;
    
    return 0;
}



// --------------------------------
// get efficiency for single line
// --------------------------------

double GetEfficiency(TString IsotopeName, double energy, double BR, double width_sig, double width_bck, int flag) {
    
    TString FileName;
    
    TString peakname;
    peakname = TString::Format("%1.1fkeV",energy);
    
    if (flag) {
        FileName = fInputFolder +"/simulation_"+peakname+".root";
    }
    else {
        FileName = fInputFolder +"/simulation_"+IsotopeName+".root";
        BR=1;
    }
    
    
    
    // open root file
    TFile* dataFile = TFile::Open(FileName);
    
    if (!dataFile) {
        std::cout << "##### ERROR: could not open " << FileName << std::endl;
        return 0;
    }

    
    // read GeHit tree
    TTree* GeHitTree = (TTree*) dataFile->Get("GeHits");
    
    // read Run tree
    TTree* RunTree = (TTree*) dataFile->Get("RunInfo");
    
    // get number of decays
    int NEvents;
    RunTree->SetBranchAddress("NEvents", &NEvents);
    RunTree->GetEntry(0);
    
    // make histogram
    TH1D* hist_TotEdep = new TH1D("hist_TotEdep","Deposited Energy per Event; Energy (keV); Counts (a.u.)", 30000, 0., 3000.);
    GeHitTree->Draw("TotEdep>>hist_TotEdep","","goff");
    
    
    // counts in peak
    double bin_min_peak = hist_TotEdep->FindBin(energy-width_sig/2.);
    double bin_max_peak = hist_TotEdep->FindBin(energy+width_sig/2.);
    double sum_peak=hist_TotEdep->Integral(bin_min_peak,bin_max_peak);
    double sum_peak_err=sqrt(sum_peak);
    int nBins_peak=bin_max_peak-bin_min_peak+1;
    
    
    // calculate counts in left bck region
    double bin_min_bckleft = hist_TotEdep->FindBin(energy-width_bck/2.-width_sig/2.);
    double bin_max_bckleft = hist_TotEdep->FindBin(energy-width_bck/2.)-1;
    double sum_bckleft=hist_TotEdep->Integral(bin_min_bckleft,bin_max_bckleft);
    double sum_bckleft_err=sqrt(sum_bckleft);
    int nBins_bckleft=bin_max_bckleft-bin_min_bckleft+1;
    
    
    // calculate counts in right bck region
    double bin_min_bckright = hist_TotEdep->FindBin(energy+width_sig/2.)+1;
    double bin_max_bckright = hist_TotEdep->FindBin(energy+width_sig/2.+width_bck/2.);
    double sum_bckright=hist_TotEdep->Integral(bin_min_bckright,bin_max_bckright);
    double sum_bckright_err=sqrt(sum_bckright);
    int nBins_bckright=bin_max_bckright-bin_min_bckright+1;
    
    // substract bck counts from peak
    double f_bckleft = nBins_peak/2. / nBins_bckleft;
    double f_bckright = nBins_peak/2. / nBins_bckright;
    
    double counts = sum_peak - f_bckleft*sum_bckleft - f_bckright*sum_bckright;
    double counts_err = sqrt(sum_peak_err*sum_peak_err + f_bckleft*sum_bckleft_err*sum_bckleft_err + f_bckright*sum_bckright_err*sum_bckright_err);
    
    double efficiency = counts/NEvents;
    
    // draw histogram
    TCanvas* c1 = new TCanvas("c1");
    gStyle->SetOptStat(0);
    c1->SetLogy();
    
    
    TH1D *hpeak = (TH1D*)hist_TotEdep->Clone();
    TH1D *hbckleft = (TH1D*)hist_TotEdep->Clone();
    TH1D *hbckright = (TH1D*)hist_TotEdep->Clone();
    
    hist_TotEdep->GetXaxis()->SetRangeUser(energy-width_sig/2.-width_bck/2.,energy+width_sig/2.+width_bck/2.);
    hist_TotEdep->Draw("hist");
    
    
    hpeak->GetXaxis()->SetRangeUser(energy-width_sig/2.,energy+width_sig/2.);
    hpeak->SetFillColor(2);
    hpeak->Draw("histsame");
    
    
    hbckleft->GetXaxis()->SetRangeUser(energy-width_sig/2.-width_bck/2.,energy-width_sig/2.);
    hbckleft->SetFillColor(3);
    hbckleft->Draw("histsame");
    
    hbckright->GetXaxis()->SetRangeUser(energy+width_sig/2.,energy+width_sig/2.+width_bck/2.);
    hbckright->SetFillColor(3);
    hbckright->Draw("histsame");
    
    char text_peak[80];
    char text_bckleft[80];
    char text_bckright[80];
    char text_amp[80];
    sprintf(text_peak,"peak region: %1.4e#pm%1.4e",sum_peak, sum_peak_err);
    sprintf(text_bckleft,"left bck region: %1.4e#pm%1.4e",sum_bckleft, sum_bckleft_err);
    sprintf(text_bckright,"right bck region: %1.4e#pm%1.4e",sum_bckright, sum_bckright_err);
    sprintf(text_amp,"peak - bck: %1.4e#pm%1.4e",counts, counts_err);
    
    TLegend* leg1 = new TLegend(0.6,0.7,0.9,0.9);
    leg1->AddEntry(hpeak,text_peak,"f");
    leg1->AddEntry(hbckleft,text_bckleft,"f");
    leg1->AddEntry(hbckright,text_bckright,"f");
    leg1->AddEntry(hist_TotEdep,text_amp,"l");
    
    leg1->Draw();
    
    TString Outputfolder = fInputFolder + "/spectra";
    
    if (!gSystem->OpenDirectory(Outputfolder)) {
        if (gSystem->MakeDirectory(Outputfolder)) {
            std::cout << "##### ERROR: could not create directory " << Outputfolder << std::endl;
            return 0;
        }
    }
    
    
    
    c1->SaveAs(Outputfolder+"/spectrum_"+IsotopeName+"_"+peakname+".pdf");
    c1->SaveAs(Outputfolder+"/spectrum_"+IsotopeName+"_"+peakname+".root");
    
    delete c1;
    
    dataFile->Close();
    
    return efficiency*BR;
    
}

// --------------------------------
// main program
// --------------------------------

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
    
    const int npeaks = fIsotopeName.size();

    double efficiency[npeaks];
    
    double Efficiency_times_BR[npeaks];
    
    for (int i=0; i<npeaks; ++i) {
        
        Efficiency_times_BR[i] = GetEfficiency(fIsotopeName[i], fEnergy[i], fBR[i], fWidth_sig[i], fWidth_bck[i], fFlag[i]);
        
        if (Efficiency_times_BR[i]==0) {
            return 1;
        }

    }

    
    TFile* file = TFile::Open(fInputFolder+"/simulated_efficiencies.root","recreate");
    
    TTree* tree = new TTree("tree","tree");
    double tenergy, teff_BR;
    
    tree->Branch("energy",&tenergy);
    tree->Branch("eff_BR",&teff_BR);
    
    for (int i = 0; i<npeaks; ++i) {
        
        efficiency[i]=Efficiency_times_BR[i]/fBR[i];
        
        tenergy=fEnergy[i];
        teff_BR=Efficiency_times_BR[i];
                
        tree->Fill();
    }

    file->cd();
    tree->Write();
    file->Close();
    
    // draw graph
    TGraph* graph = new TGraph(npeaks, &fEnergy.at(0), efficiency);
    graph->SetName("graph_efficiency");
    
    TCanvas* c2 = new TCanvas("c2");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(1);
    graph->SetLineColor(1);
    graph->SetTitle("Simulated Efficiency");
    graph->GetXaxis()->SetTitle("Energy (keV)");
    graph->GetYaxis()->SetTitle("Detection Efficiency");
    graph->Draw("ap");
    
    c2->SaveAs(fInputFolder+"/energy_vs_efficiency.pdf");
    c2->SaveAs(fInputFolder+"/energy_vs_efficiency.root");
    
    return 0;
    
}





