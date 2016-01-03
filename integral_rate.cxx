#include "macros.h"



int main(int argc, char *argv[]) {
    
    // argc should be >=4 for correct execution
    if ( argc != 4 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <spectrum.root> <E_min> <E_max>" << std::endl;
        
        return 1;
    }
    
    TString filename = argv[1];
    
    double E_min = atof(argv[2]);
    double E_max = atof(argv[3]);

    
    integralrate(filename,E_min,E_max);
    
    return 0;
    
    
}



