#include "macros.h"
#include <getopt.h>


// ----------------------------------------------------
// plot rate vs time
// ----------------------------------------------------

int main(int argc, char *argv[]) {
    
    // argc should be >=2 for correct execution
    if ( argc < 2 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <arguments> " << std::endl;
        std::cout<<"arguments:" << std::endl;
        std::cout<<"--file <list_file.root> \t list file" << std::endl;
        std::cout<<"--binwidth <binwidth> \t bin width (s)" << std::endl;
        std::cout<<"--energy \t use energy calibration" << std::endl;
        std::cout<<"--range_min <r0> --range_max <r1> \t select pulseheight/energy range [r0,r1]" << std::endl;

        return 1;
    }
    
    int c = 0;
    TString option="";
    double binwidth=2000.;
    double range_min=100.;
    double range_max=30000.;
    TString FileName;
    struct option longopts[] = {
        { "energy",         no_argument,     0,     'e' },
        { "file",         required_argument, 0,     'f' },
        { "range_min",    required_argument, 0,     'a' },
        { "range_max",    required_argument, 0,     'b' },
        { "binwidth",     required_argument, 0,     'c' },
        { 0, 0, 0, 0 }
    };
    while((c = getopt_long(argc,argv,"ea:b:c:",longopts, NULL)) != -1)
    {
        switch(c)
        {
            case 'e':
                option="energy";
                break;
                
            case 'f':
                FileName = optarg;
                break;
                
            case 'a':
                range_min = atof(optarg);
                break;
                
            case 'b':
                range_max = atof(optarg);
                break;
            
            case 'c':
                binwidth = atof(optarg);
                break;
                
            case ':':   /* missing option argument */
                fprintf(stderr, "%s: option `-%c' requires an argument\n",
                        argv[0], optopt);
                break;

        }
    }

    if (FileName=="") {
        std::cout << "##### ERROR: no file name specified" << std::endl;
        std::cout << "use option --file <list_file.root>" << std::endl;
        return 1;
    }
    
    // read parameters from file
    if (rate_from_list(FileName,binwidth,range_min,range_max,option)) {
        return 1;
    }
    
    return 0;
}


