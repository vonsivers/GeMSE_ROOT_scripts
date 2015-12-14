#include "macros.h"
#include <getopt.h>


// ----------------------------------------------------
// make root file from list file
// ----------------------------------------------------

int main(int argc, char *argv[]) {
    
    // argc should be >=2 for correct execution
    if ( argc < 2 ) {
        
        // We print argv[0] assuming it is the program name
        std::cout<<"usage: "<< argv[0] <<" <arguments>" << std::endl;
        std::cout<<"arguments:" << std::endl;
        std::cout<<"--file <list_file.root> \t list file" << std::endl;
        std::cout<<"--energy \t use energy calibration" << std::endl;
        std::cout<<"--t_min <t0> --t_max <t1> \t select time range [t0,t1]" << std::endl;

        return 1;
    }
    
    int c = 0;
    TString option="";
    TString FileName;
    double t_min=0.;
    double t_max=0.;
    struct option longopts[] = {
        { "energy",     no_argument,     0,     'e' },
        { "file",     required_argument, 0,     'f' },
        { "t_min",    required_argument, 0,     'a' },
        { "t_max",    required_argument, 0,     'b' },
        { 0, 0, 0, 0 }
    };
    while((c = getopt_long(argc,argv,"ea:b:f:",longopts, NULL)) != -1)
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
                t_min = atof(optarg);
                break;
                
            case 'b':
                t_max = atof(optarg);
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

    
    if (spectrum_from_list(FileName,option,t_min,t_max)) {
        return 1;
    }
    
    return 0;
}


