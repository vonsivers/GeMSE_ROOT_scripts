#include "macros.h"

// ----------------------------------------------------
// make root file from ascii spectrum
// ----------------------------------------------------

int main(int argc, char *argv[]) {
    
    // argc should be 2 for correct execution
    if ( argc != 2 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <spectrum_file.txt>"<< std::endl;
        return 1;
    }
    
    TString FileName = argv[1];
    
    // read parameters from file
    if (make_spectrum(FileName)) {
        return 1;
    }
    
    return 0;
}


