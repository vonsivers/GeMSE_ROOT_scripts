#include <TChain.h>
#include <TString.h>
#include <TCanvas.h>

#include <iostream>

// --------------------------------
// main program
// --------------------------------

int main(int argc, char *argv[]) {
    
    // argc should be > 3 for correct execution
    if ( argc < 3 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <simulated_efficiencies1.root> <simulated_efficiencies2.root> ... <simulated_efficiencies_merged>" << std::endl;
        return 1;
    }
    
    // name of output file
    TString FileName_out = argv[argc-1];
    
    // create TChain for all trees
    TChain* chain = new TChain("tree");
    
    // loop over all files
    for (int i=0; i<argc-2; ++i) {
        
        TString FileName = argv[i+1];
        chain->Add(FileName);
    
    }
    
    // merge chain into single file
    chain->Merge(FileName_out+".root");

    // plot efficiencies
    TCanvas* c = new TCanvas("c");
    chain->Draw("efficiency:energy","","*");
    c->SaveAs(FileName_out+".pdf");
    
    
}





