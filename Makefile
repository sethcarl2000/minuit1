# This Makefile builds two stand alone programs that have access to 
# the ROOT libraries.  To use this with your own code just substitute
# the name of your program below.

# here we access the root configuration, include files, and libraries
ROOTCFLAGS=$(shell root-config --cflags)
ROOTINC=$(shell root-config --incdir)
ROOTLIBDIR=$(shell root-config --libdir)
ROOTLIBS=$(shell root-config --libs) -lMinuit
ROOTLDFLAGS=$(shell root-config --ldflags)

ROOTC=$(ROOTCFLAGS) 
#-I$(ROOTINC)
ROOTLINK=-L$(ROOTLIBDIR) $(ROOTLIBS) $(ROOTLDFLAGS)

CPP=g++

default: expFit rootExample fit_distros fit_experiments fit_2d

expFit: expFit.cpp
	$(CPP) -O -Wall $(ROOTC) -o expFit expFit.cpp $(ROOTLINK) 

fit_distros: fit_distros.cpp
	$(CPP) -O -Wall $(ROOTC) -o fit_distros fit_distros.cpp $(ROOTLINK) 

fit_experiments: fit_experiments.cpp
	$(CPP) -O -Wall $(ROOTC) -o fit_experiments fit_experiments.cpp $(ROOTLINK) 

fit_experiments: fit_2d.cpp
	$(CPP) -O -Wall $(ROOTC) -o fit_2d fit_2d.cpp $(ROOTLINK) 

rootExample: rootExample.cpp
	$(CPP) -O -Wall $(ROOTC) -o rootExample rootExample.cpp $(ROOTLINK) 
# note: just replace the -O flag with -g to build a debug version

clean: 
	rm -f expFit rootExample *~ *.d *.so
