#include <iostream>
#include <fstream>
#include <cmath>
#include "Constants.h"
#include "utility.h"
using namespace std;

int main(){
	double m1 = 1.0, m2 = 1.1, ep1 = 1e-5, ep2 = 1e-9, R1 = 10e3, R2 = 30e3, L_n = 3e52/* erg/s */*624151/* MeV/s */;
	double range = 10e3/1.0;
	int n_points = 1000;
	double interval = pow(range, 1.0/n_points);
	double eng = 1.0;

	//Input density and temperature profile
	double r[1189], dens[1189], temp[1189];
	ifstream input;
	input.open("../Degenerate/Density_Temperature.txt");
	for(int i = 0; i < 1189; i++){
		input >> r[i] >> dens[i] >> temp[i];
	}
	input.close();	

	ofstream out1;
	out1.open("Dfr_Lum_1MeV.txt");
	for(int i = 0; i <= n_points; i++){
		out1 << eng << "\t"
		     << LumInt_general(eng, R1, m1, ep1, R_far, r, dens, temp, 1189)*1000/L_n << "\t"
		     << LumInt_general(eng, R2, m1, ep1, R_far, r, dens, temp, 1189)*1000/L_n << "\t"
		     << LumInt_general(eng, R1, m1, ep2, R_far, r, dens, temp, 1189)*1000/L_n << "\t"
		     << LumInt_general(eng, R2, m1, ep2, R_far, r, dens, temp, 1189)*1000/L_n << "\t"
		     << PowInt_general(eng, R1, m1, ep2, r, dens, temp, 1189)*1000/L_n << "\t"
		     << PowInt_general(eng, R2, m1, ep2, r, dens, temp, 1189)*1000/L_n << "\t"
		     << 	   Tau(eng, R1, m1, ep2, R_far, r, dens, temp, 1189) << "\t"
		     << 	   Tau(eng, R2, m1, ep2, R_far, r, dens, temp, 1189) << "\t"
		     << endl;
		eng *= interval;
	}
	out1.close();

	eng = 1.1;
	ep1 = 1e-6;
	ep2 = 1e-8;

	ofstream out2;
	out2.open("Dfr_Lum_1_1MeV.txt");
	for(int i = 0; i <= n_points; i++){
		out2 << eng << "\t"
		     << LumInt_general(eng, R1, m2, ep1, R_far, r, dens, temp, 1189)*1000/L_n << "\t"
		     << LumInt_general(eng, R2, m2, ep1, R_far, r, dens, temp, 1189)*1000/L_n << "\t"
		     << LumInt_general(eng, R1, m2, ep2, R_far, r, dens, temp, 1189)*1000/L_n << "\t"
		     << LumInt_general(eng, R2, m2, ep2, R_far, r, dens, temp, 1189)*1000/L_n << "\t"
		     << PowInt_general(eng, R1, m2, ep2, r, dens, temp, 1189)*1000/L_n << "\t"
		     << PowInt_general(eng, R2, m2, ep2, r, dens, temp, 1189)*1000/L_n << "\t"
		     << Tau(eng, R1, m2, ep2, R_far, r, dens, temp, 1189) << "\t"
		     <<	Tau(eng, R2, m2, ep2, R_far, r, dens, temp, 1189) << "\t"
		     << endl;
		eng *= interval;
	}
	out2.close();
	return 0; 
}

