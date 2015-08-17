#include "macros.h"

// ----------------------------------------------------
// calibrate spectrum
// ----------------------------------------------------

int main(int argc, char *argv[]) {
    
    // argc should be 2 for correct execution
    if ( argc != 3 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <spectrum_file.txt>" <<" <calibration_function.root>" << std::endl;
        return 1;
    }
    
    TString FileName_spectrum = argv[1];
    TString FileName_calibration = argv[2];
    
    // get calibration function from root file
    TF1* calibration = getcalibration(FileName_calibration);
    
    if (calibration==0) {
        return 1;
    }
 
    // calibrate spectrum
    if(make_spectrum(FileName_spectrum,calibration)) {
        return 1;
    }
    
    return 0;
}


