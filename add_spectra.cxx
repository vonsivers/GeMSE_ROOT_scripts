#include "macros.h"


// ----------------------------------------------------
// add up spectra
// ----------------------------------------------------

int main(int argc, char *argv[]) {
    
    // argc should be >=4 for correct execution
    if ( argc < 4 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <spectrum1.root> <spectrum2.root> ... <results_filename>" << std::endl;
        
        return 1;
    }
    
    std::vector<TString> Files;
    
    for (int i=1; i<argc-1; ++i) {
        Files.push_back(argv[i]);
    }
    
    TString results_filename = argv[argc-1];
    
    if (addspectra(Files,results_filename)) {
        return 1;
    }
    
    return 0;
}


