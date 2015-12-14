#include "macros.h"

// ----------------------------------------------------
// make root file from list file
// ----------------------------------------------------

int main(int argc, char *argv[]) {
    
    // argc should be 2 for correct execution
    if ( argc < 2 || argc > 3 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <list_file.dat>" <<" <calibration_function.root>" << std::endl;
        return 1;
    }
    
    TString FileName_list = argv[1];
    
    if (argc==3) {
        TString FileName_calibration = argv[2];
        
        // get calibration function from root file
        TF1* calibration = getcalibration(FileName_calibration);
        
        if (calibration==0) {
            std::cout << "###### ERROR: no calibration function found in " << FileName_calibration << std::endl;
            return 1;
        }
        
        // calibrate spectrum
        if(read_listfile(FileName_list,calibration)) {
            return 1;
        }

    }
    
    else {
        if(read_listfile(FileName_list)) {
            return 1;
        }
    }
    
    return 0;
}


