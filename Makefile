###################################################################
# This Makefile was created using the bat-project script.
# bat-project is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://mpp.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CXXFLAGS and LIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# compiler and flags
CXX          = g++
CXXFLAGS     = -g -O2 -Wall -fPIC -Wno-deprecated 
LD           = /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ld
LDFLAGS      = -g -O2 

# ----------------------------------------------------------------------
# The following definitions rely on the script bat-config being
# available in $PATH. If BAT is not installed in the standard system
# directories, update $PATH accordingly.

CXXFLAGS += `bat-config --cflags`
LIBS := `bat-config --libs`

# List of all classes (models) used in the program
# Add classes to the end. Backslash indicates continuation
# on the next line
CXXSRCS_all      = \
	calibrate_spectrum.cxx energy_calibration.cxx energy_resolution.cxx make_rootfile_list.cxx make_rootfile_spectrum.cxx make_spectrum_list.cxx plot_rate.cxx simulated_efficiency.cxx
CXXSRCS1      = \
        calibrate_spectrum.cxx 
CXXSRCS2      = \
	energy_calibration.cxx
CXXSRCS3      = \
	energy_resolution.cxx
CXXSRCS4      = \
	make_rootfile_list.cxx
CXXSRCS5      = \
	make_rootfile_spectrum.cxx
CXXSRCS6      = \
	make_spectrum_list.cxx
CXXSRCS7      = \
	plot_rate.cxx
CXXSRCS8      = \
	simulated_efficiency.cxx
# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#
CXXOBJS_all      = $(patsubst %.cxx,%.o,$(CXXSRCS_all))

CXXOBJS1      = $(patsubst %.cxx,%.o,$(CXXSRCS1))
CXXOBJS2      = $(patsubst %.cxx,%.o,$(CXXSRCS2))
CXXOBJS3      = $(patsubst %.cxx,%.o,$(CXXSRCS3))
CXXOBJS4      = $(patsubst %.cxx,%.o,$(CXXSRCS4))
CXXOBJS5      = $(patsubst %.cxx,%.o,$(CXXSRCS5))
CXXOBJS6      = $(patsubst %.cxx,%.o,$(CXXSRCS6))
CXXOBJS7      = $(patsubst %.cxx,%.o,$(CXXSRCS7))
CXXOBJS8      = $(patsubst %.cxx,%.o,$(CXXSRCS8))

MYPROGS_all     = \
	calibrate_spectrum energy_calibration energy_resolution make_rootfile_list make_rootfile_spectrum make_spectrum_list plot_rate simulated_efficiency

MYPROGS1     = \
        calibrate_spectrum
MYPROGS2     = \
	energy_calibration
MYPROGS3     = \
	energy_resolution
MYPROGS4     = \
	make_rootfile_list
MYPROGS5     = \
	make_rootfile_spectrum
MYPROGS6     = \
	make_spectrum_list
MYPROGS7     = \
	plot_rate
MYPROGS8     = \
	simulated_efficiency

GARBAGE = $(CXXOBJS_all) *~ link.d $(MYPROGS_all)

# targets
all : calibrate_spectrum energy_calibration energy_resolution make_rootfile_list make_rootfile_spectrum make_spectrum_list plot_rate simulated_efficiency

link.d : $(patsubst %.cxx,%.h,$(CXXSRCS_all))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS_all) > link.d;

-include link.d

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	rm -f $(GARBAGE)

calibrate_spectrum : $(CXXOBJS1)
	$(CXX) $(LDFLAGS) $(CXXOBJS1) $(LIBS) -o $@

energy_calibration : $(CXXOBJS2)
	$(CXX) $(LDFLAGS) $(CXXOBJS2) $(LIBS) -o $@

energy_resolution : $(CXXOBJS3)
	$(CXX) $(LDFLAGS) $(CXXOBJS3) $(LIBS) -o $@

make_rootfile_list : $(CXXOBJS4)
	$(CXX) $(LDFLAGS) $(CXXOBJS4) $(LIBS) -o $@

make_rootfile_spectrum : $(CXXOBJS5)
	$(CXX) $(LDFLAGS) $(CXXOBJS5) $(LIBS) -o $@

make_spectrum_list : $(CXXOBJS6)
	$(CXX) $(LDFLAGS) $(CXXOBJS6) $(LIBS) -o $@

plot_rate : $(CXXOBJS7)
	$(CXX) $(LDFLAGS) $(CXXOBJS7) $(LIBS) -o $@

simulated_efficiency : $(CXXOBJS8)
	$(CXX) $(LDFLAGS) $(CXXOBJS8) $(LIBS) -o $@


print :
	@echo compiler  : $(CXX)
	@echo c++ srcs  : $(CXXSRCS_all)
	@echo c++ objs  : $(CXXOBJS_all)
	@echo c++ flags : $(CXXFLAGS)
	@echo libs      : $(LIBS)



