###################################################################
make_rootfile
###################################################################
converts an ASCI spectrum of CAEN DT5781A to a ROOT file
first line in ASCI file must be real time, second line live time
usage
	./make_rootfile <spectrum_file.txt>
output
	<spectrum_file.txt.root>
	contains the spectrum as a TH1D and the live time and real time both as TVectorD

###################################################################
energy_calibration
###################################################################
fits the peaks in a spectrum with a Gaussian
the peak position is then fitted with a second order polynomial to get the energy calibration
requires a parameters txt file (see „example_parameters_calibration.txt“) with the following content:
	-name of spectrum root file
	also for every peak
	-fit range, fit start parameters (amplitude, mean, sigma), literature value of peak energy
usage
	./energy_calibration <parameters_calibration.txt>
output	
	<spectrum_file.txt.root_calibration_fits.root>
	contains the original spectrum with all Gaussian fits to the peaks
	
	<spectrum_file.txt.root_calibration_function.root>
	contains a TGraphErrors with the peak position vs energy and the polynomial fit function to the data


###################################################################
calibrate_spectrum
###################################################################
calibrates a spectrum txt file with a calibration function
usage
	./calibrate_spectrum <spectrum_file.txt> <calibration_function.root>
output
	<spectrum_file.txt_calibrated.root>
	contains the calibrated spectrum as a TH1D and the live time and real time both as TVectorD

###################################################################
simulated_efficiency
###################################################################
analyzes the simulated efficiencies
requires a paramters txt file (see "example_parameters_efficiency_simulation.txt")
usage 
    ./simulated_efficiency <parameters_simulated_efficiency.txt>
output
    simulated_efficiencies.root
    contains a TTree with branches "energy" containing the energy in keV and "eff_BR" containing the product of detection efficiency and branching ratio
    in the folder "/spectra" the spectra of all simulated gamma lines are written as ROOT and pdf files


###################################################################
macros.h
###################################################################
contains all macros which are used by the other scripts