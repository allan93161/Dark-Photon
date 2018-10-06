Dark_Photon.C: Contains the calculations of all three absorption rates as attributes of the class Dark_Photon.

utility.C: Contains the calculations of all the other quantities, such as optical depth, differential power and luminosity, integrated luminosity, etc.

DfrLum.C: Plots differential power and luminosity as functions of #omega (dark photon energy) as Fig. 5, along with the optical depthes as functions of #omega. This program generates two text files, "DfrLum_1MeV.txt" and "DfrLum_1_1MeV.txt," that can be used to draw the same plots as Fig. 5 by other plotting software.

Before compiling, one needs to set up the path for GSL. Open Makefile, go to the 2nd line, and change the include directory -I<include path> and the library -L<lib path> to wherever GSL is located.
To compile, type "make" and copy the two lines to the command line:

	LD_LIBRARY_PATH=/usr/local/lib
	export LD_LIBRARY_PATH

Then, you can go ahead and run ./DfrLum.
