//Constants.h
#if !defined(MYLIB_CONSTANTS_H)
#define MYLIB_CONSTANTS_H 1

const double Alpha = 1/137.036;		//Fine structure constant
const double M_n = 939.565;		//Neutron mass [MeV]
const double M_e = 0.511;		//Electron mass[MeV]
const double Hbar = 6.582e-22;		//Reduced Planck constant [MeV-s]
const double C = 2.998e8;		//Speed of light [m/s]
const double PI = 3.14159265;		//Pi
const double R_n = 40e3;		//Radius inside of which the dark photons are produced
					//~radius of neutrinosphere (40km)
const double R_far = 100e3;		//Radius beyond which the dark photons must travel before decay
					//~shock radius(100km)
const double Y = 0.3;			//Proton fraction
const double Rho_c = 1.69e47;		//Core density [MeV/m^3]
const double T_c = 30;			//Core temperature [MeV]
const double K_rho = 0.2;		//Parameter of density profile
const double K_t = -0.5;		//Parameter of temperature profile
const int Neu = 5;			//Power law of density/temperature profile
const double R_c = 10e3;		//Core radius

struct my_params{double m_; double ep; double R; double* r; double* dens; double* temp; int N;};
struct res_params{double w; double m;};
struct decay_params{double eng; double mu_e; double temp; double m;};

#endif
