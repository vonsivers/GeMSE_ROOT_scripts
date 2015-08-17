#include <TH1D.h>
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

#include <iostream>
#include <fstream>
#include <stdlib.h>     /* atof */

// erase comma from number
double getNumber(std::string a) {
    
    for(int i=0; i<=a.size(); i++) {  //Loop through the entire string
        if(a[i]==',') { //Check for occurence of comma
            a.erase(i,1); //If comma , remove that single character
        }
    }
    
    return atof(a.c_str()); //Will return double

}

// make root file from ascii spectrum
int make_spectrum(TString FileName, TF1* calibration=0) {
    
    // open ascii file
    ifstream File;
    
    File.open(FileName);
    
    if (!File.is_open()) {
        std::cout << "##### ERROR: could not open " << FileName << std::endl;
        return 1;
    }
    
    double t_live, t_real;
    std::string t_live_string, t_real_string;
    
    const int Nchannels=16382;
    double channel[Nchannels+1];
    int counts[Nchannels];
    
    // read real and live time
    File >> t_real_string >> t_live_string;
    
    // remove commas and convert to double
    t_real = getNumber(t_real_string);
    t_live = getNumber(t_live_string);
    
    
    // read counts for every channel
    for (int i=0; i<Nchannels; i++) {
        File >> counts[i];
        
        // apply energy calibration
        if (calibration!=0) {
            channel[i]=calibration->Eval(i);
        }
        else {
            channel[i]=i;
        }
    }
    if (calibration!=0) {
        channel[Nchannels]=calibration->Eval(Nchannels);
    }
    else {
        channel[Nchannels]=Nchannels;
    }

    
    File.close();
    
    // create root file
    if (calibration!=0) {
        FileName += "_calibrated";
    }
    
    TFile* histFile = TFile::Open(FileName+".root","RECREATE");
    
    // make histogram
    TH1D* hist = new TH1D("hist",";ADC Channel;Counts",Nchannels,channel);
    
    if (calibration!=0) {
        hist->GetXaxis()->SetTitle("Energy (keV)");
    }
    
    for (int i=0; i<Nchannels; i++) {
        
        hist->SetBinContent(i+1,counts[i]);
        
    }
    
    // draw histogram
    TCanvas* c1 = new TCanvas("c1");
    gStyle->SetOptStat(0);
    c1->SetLogy();
    hist->Draw();
    c1->SaveAs(FileName+".pdf");
    
    
    // save histogram to root file
    hist->Write();
    TVectorD v_live(1);
    TVectorD v_real(1);
    v_live[0] = t_live;
    v_real[0] = t_real;
    v_live.Write("t_live");
    v_real.Write("t_real");
    
    histFile->Close();
    
    return 0;
    
    
}

// get spectrum from root file
TH1D* getspectrum(TString FileName) {
    
    // check for root file
    TFile* File = TFile::Open(FileName);
    
    if (!File) {
        
        std::cout << "###### ERROR: could not open " << FileName << std::endl;
        return 0;
    }
    
    TH1D* hist = (TH1D*) File->Get("hist");
    hist->SetDirectory(0);
    File->Close();
    
    return hist;
    
}

// get spectrum from root file
TF1* getcalibration(TString FileName) {
    
    // check for root file
    TFile* File = TFile::Open(FileName);
    
    if (!File) {
        
        std::cout << "###### ERROR: could not open " << FileName << std::endl;
        return 0;
    }
    TCanvas* c = (TCanvas*) File->Get("c2");
    TGraphErrors* graph = (TGraphErrors*) c->GetPrimitive("Graph");
    TF1* fit = (TF1*) graph->GetFunction("fitFunction");
    
    File->Close();
    
    return fit;
    
}

// fit peak only
TF1* FitGauss(TH1D* hist,double amp_st, double mean_st, double sigma_st, double range_min, double range_max) {
    
    TObject* func = gROOT->FindObject("fitFunction");
    
    if (func) {
        std::cout << "#### 'fitFunction' already defined, deleting ..." << std::endl;
        delete func;
    }
    TF1* fitFunction = new TF1("fitFunction", "[0]/([2]*sqrt(2*3.14159265))*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))");
    
    fitFunction->SetParName(0,"amplitude");
    fitFunction->SetParName(1,"mean");
    fitFunction->SetParName(2,"sigma");
    
    fitFunction->SetRange(range_min,range_max);
    
    fitFunction->SetParameter(0,amp_st);
    fitFunction->SetParameter(1,mean_st);
    fitFunction->SetParameter(2,sigma_st);
    
    //fitFunction->SetParLimits(0,0.,1000.*amp_st);
    //fitFunction->SetParLimits(1,range_min,range_max);
    //fitFunction->SetParLimits(2,0.,1000.*sigma_st);
    
    // fit signal
    hist->Fit(fitFunction,"0RMEL");
    
    fitFunction->DrawCopy("LSAME");
    
    
    return fitFunction;
    
}

// fit peak + background
TF1* FitPeak(TH1D* hist,double amp_st, double mean_st, double sigma_st, double const_st, double slope_st, double range_min, double range_max) {
    
    TObject* func = gROOT->FindObject("fitFunction");
    
    if (func) {
        std::cout << "#### 'fitFunction' already defined, deleting ..." << std::endl;
        delete func;
    }
    TF1* fitFunction = new TF1("fitFunction", "[0]/([2]*sqrt(2*3.14159265))*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))+pol1(3)");

    
    fitFunction->SetParName(0,"amplitude");
    fitFunction->SetParName(1,"mean");
    fitFunction->SetParName(2,"sigma");
    fitFunction->SetParName(3,"constant");
    fitFunction->SetParName(4,"slope");
    
    fitFunction->SetRange(range_min,range_max);
    
    fitFunction->SetParameter(0,amp_st);
    fitFunction->SetParameter(1,mean_st);
    fitFunction->SetParameter(2,sigma_st);
    fitFunction->SetParameter(3,const_st);
    
    //fitFunction->SetParameter(4,slope_st);
    
    //fitFunction->SetParLimits(0,0.,1000.*amp_st);
    //fitFunction->SetParLimits(1,range_min,range_max);
    //fitFunction->SetParLimits(2,0.,1000.*sigma_st);
    
    // fit signal
    hist->Fit(fitFunction,"0RMEL");
    
    fitFunction->DrawCopy("LSAME");
    
    
    return fitFunction;

}

// fit 2 peaks + background
TF1* Fit2Peak(TH1D* hist,double amp1_st, double mean1_st, double sigma1_st, double amp2_st, double mean2_st, double sigma2_st, double const_st, double slope_st, double range_min, double range_max) {
    
    TObject* func = gROOT->FindObject("fitFunction");
    
    if (func) {
        std::cout << "#### 'fitFunction' already defined, deleting ..." << std::endl;
        delete func;
    }

    
    TF1* fitFunction = new TF1("fitFunction", "[0]/([2]*sqrt(2*3.14159265))*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))+[3]/([5]*sqrt(2*3.14159265))*exp(-(x-[4])*(x-[4])/(2*[5]*[5]))+pol1(6)");
    
    fitFunction->SetParName(0,"amplitude1");
    fitFunction->SetParName(1,"mean1");
    fitFunction->SetParName(2,"sigma1");
    fitFunction->SetParName(3,"amplitude2");
    fitFunction->SetParName(4,"mean2");
    fitFunction->SetParName(5,"sigma2");
    fitFunction->SetParName(6,"constant");
    fitFunction->SetParName(7,"slope");
    
    fitFunction->SetRange(range_min,range_max);
    
    fitFunction->SetParameter(0,amp1_st);
    fitFunction->SetParameter(1,mean1_st);
    fitFunction->SetParameter(2,sigma1_st);
    fitFunction->SetParameter(3,amp2_st);
    fitFunction->SetParameter(4,mean2_st);
    fitFunction->SetParameter(5,sigma2_st);
    fitFunction->SetParameter(6,const_st);
    fitFunction->SetParameter(7,slope_st);
    
    // fitFunction->SetParLimits(0,0.,1000.);
    
    // fit signal
    hist->Fit(fitFunction,"0RMWL");
    
    fitFunction->DrawCopy("LSAME");
    
    return fitFunction;
    
}

// fit background only
TF1* FitPol1(TH1D* hist,double a_st, double b_st, double range_min, double range_max) {
    
    TF1* fitFunction = new TF1("fitFunction", "pol1");
    
    fitFunction->SetParameter(0,a_st);
    fitFunction->SetParameter(1,b_st);
    
    fitFunction->SetRange(range_min,range_max);
    
    // fit graph
    hist->Fit(fitFunction,"0REMF");
    
    gStyle->SetOptFit(1);
    
    fitFunction->DrawCopy("LSAME");
    
    return fitFunction;
}

// fit energy calibration
TF1* FitPol2(TGraphErrors* graph,double a_st, double b_st, double c_st) {

    TF1* fitFunction = new TF1("fitFunction", "pol2");
    
    fitFunction->SetParameter(0,a_st);
    fitFunction->SetParameter(1,b_st);
    fitFunction->SetParameter(2,c_st);

    // fit graph
    graph->Fit(fitFunction,"EMF");
    
    gStyle->SetOptFit(1);
    
    return fitFunction;
}


// fit resolution
TF1* FitSqrt(TGraphErrors* graph,double a_st, double b_st, double c_st) {
    
    TF1* fitFunction = new TF1("fitFunction", "sqrt([0]+[1]*x+[2]*x*x)");
    
    fitFunction->SetParameter(0,a_st);
    fitFunction->SetParameter(1,b_st);
    fitFunction->SetParameter(2,c_st);
    
    // fit graph
    graph->Fit(fitFunction,"EMF");
    
    gStyle->SetOptFit(1);
    
    return fitFunction;
}


// make root file from list file
void read_listfile(TString FileName) {
    
    
    // check for root file
    TFile* rootFile = new TFile(FileName+".root","recreate");
    
    TTree* dataTree = new TTree("dataTree","dataTree");
    TTree* headerTree =new TTree("headerTree","headerTree");
    
    long long int time, time_raw;
    int pulseheight;
    long long int extras;
    std::string headerline;

    dataTree->Branch("time",&time);
    dataTree->Branch("pulseheight",&pulseheight);
    dataTree->Branch("extras",&extras);
    
    headerTree->Branch("headerline",&headerline);
    
    ifstream File;
    
    //ifstream myFile ("myfile.dat", ios::in | ios::binary);

    
    File.open(FileName);
    
    std::cout << "Reading list file ..." << std::endl;
    
    // start reading from file
    for (int i=0; i<5; ++i) {
        File >> headerline ;
        //std::cout << headerline << std::endl;
        headerTree->Fill();
    }

    // read histogram
    while (true) {
        File >> time_raw >> pulseheight >> extras;
        //std::cout << time_raw << "\t" << pulseheight << "\t" << extras << std::endl;
        time = time_raw*10; // convert time to ns
        dataTree->Fill();
        if( File.eof() ) break;
    }
    
    File.close();
    
    headerTree->Write();
    dataTree->Write();
    
    rootFile->Close();
    
    
}

// make spectrum from list file
TH1D* spectrum_from_list(TString FileName, TString option="") {
    
    
    // check for root file
    TFile* File = new TFile(FileName+".root");
    
    // create root file from ASCII file if necessary
    if (File->IsZombie()) {
        
        File->Close();
        
        std::cout << "creating new root file ..." << std::endl;
        
        read_listfile(FileName);
        
        TFile* File = TFile::Open(FileName+".root");

    }
    
    TTree* dataTree = (TTree*) File->Get("dataTree");
    
    // get real time
    long long int time;
    
    dataTree->SetBranchAddress("time",&time);
    int nEntries = dataTree->GetEntries();
    
    dataTree->GetEntry(0);
    long long int t_start = time;

    dataTree->GetEntry(nEntries-2);
    long long int t_end = time;

    double t_real = (double)(t_end-t_start)/(1.e9);
    
    // get spectrum
    TH1D* hist = new TH1D("hist",";ADC Channel;Counts",32768,0,32768);
    
    dataTree->Draw("pulseheight>>hist","extras==0","goff");
    hist->SetDirectory(0);
    File->Close();
    
    // draw spectrum
    TCanvas* c1 = new TCanvas("c1");
    gStyle->SetOptStat(0);
    c1->SetLogy();
    
    TH1D *hdraw = (TH1D*)hist->Clone("hdraw");
    if (option=="scale") {
        hdraw->Sumw2();
        hdraw->Scale(1./t_real*86400,"width");
        hdraw->GetYaxis()->SetTitle("Counts/bin/day");
    }
    hdraw->Draw();
    
    c1->SaveAs(FileName+"_spectrum.pdf");
    
   
    
    
    // write to file
    TFile* histFile = new TFile(FileName+"_spectrum.root","recreate");
    hist->Write();
    TVectorD v_real(1);
    v_real[0] = t_real;
    v_real.Write("t_real");
    histFile->Close();
    
    return hist;
}

// make time spectrum from list file
TH1D* time_from_list(TString FileName, double binwidth) {
    
    
    // check for root file
    TFile* File = new TFile(FileName+".root");
    
    // create root file from ASCII file if necessary
    if (File->IsZombie()) {
        
        File->Close();
        
        std::cout << "creating new root file ..." << std::endl;
        
        read_listfile(FileName);
        
        TFile* File = TFile::Open(FileName+".root");
        
    }
    
    TTree* dataTree = (TTree*) File->Get("dataTree");
    
    double t_max = (double) dataTree->GetMaximum("time")/1.e9;
    
    int nbins = t_max/binwidth;
    TH1D* hist = new TH1D("hist","",nbins,0,t_max);
    
    dataTree->Draw("time/1.e9>>hist","","goff");
    hist->SetDirectory(0);
    File->Close();
    
    // scale rate to Hz
    hist->Scale(1./binwidth);
    
    // draw spectrum
    TCanvas* c1 = new TCanvas("c1");
    gStyle->SetOptStat(0);
    
    hist->GetXaxis()->SetTitle("Time (s)");
    hist->GetYaxis()->SetTitle("Rate [Hz]");
    
    hist->Sumw2();
    
    hist->Draw();
    
    c1->SaveAs(FileName+"_rate_vs_time.pdf");
    
    
    // write to file
    TFile* histFile = new TFile(FileName+"_rate_vs_time.root","recreate");
    hist->Write();
    histFile->Close();
    
    return hist;
}



