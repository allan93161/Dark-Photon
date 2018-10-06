CXX  = g++
CXXFLAGS = -Wall -g -std=c++11 -I/usr/local/include -L/usr/local/lib

all: DfrLum.o utility.o Dark_Photon.o
	$(CXX) $(CXXFLAGS) DfrLum.o utility.o Dark_Photon.o -o DfrLum -lgsl -lgslcblas -lm
	LD_LIBRARY_PATH=/usr/local/lib
	export LD_LIBRARY_PATH

utility.o: utility.C utility.h
	$(CXX) $(CXXFLAGS) -c utility.C -lgsl -lgslcblas -lm

Dark_Photon.o: Dark_Photon.C Dark_Photon.h
	$(CXX) $(CXXFLAGS) -c Dark_Photon.C -lgsl -lgslcblas -lm

DfrLum.o: DfrLum.C
	$(CXX) $(CXXFLAGS) -c DfrLum.C -lgsl -lgslcblas -lm

clean: 
	rm *.o DfrLum
