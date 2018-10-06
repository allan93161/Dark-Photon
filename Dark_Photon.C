//Implementation of class Dark_Photon
#include<iostream>
#include<cmath>
#include <gsl/gsl_integration.h>
#include "Dark_Photon.h"
#include "Constants.h"
using namespace std;

//Bremsstrahlung
double Dark_Photon::Brem_t(double eng, double temp, double dens){
	if(eng <= m) return 0;
	double Eps;
	if(eps == 1) Eps = 1;
	else Eps = Eps_t(eng, temp, dens);
	double n_n = (1-Y)*dens/M_n;
	double n_p = (Y*dens)/M_n;
	double sig_np = sig_np_avg(temp);
	double pre = (32/(3*PI)) * (Alpha*Eps*Eps*n_n*n_p/pow(eng, 3)) * pow(PI*temp/M_n, 1.5);
	return pre*sig_np*1e-31*pow(Hbar*C, 4);
}

double Dark_Photon::Brem_l(double eng, double temp, double dens){
	if(eng <= m) return 0;
	double Eps;
	if(eps == 1) Eps = 1;
	else Eps = Eps_l(eng, temp, dens);
	double n_n = (1-Y)*dens/M_n;
	double n_p = (Y*dens)/M_n;
	double sig_np = sig_np_avg(temp);
	double pre = (32/(3*PI)) * (Alpha*Eps*Eps*n_n*n_p/pow(eng, 3)) * pow(PI*temp/M_n, 1.5);
	return pre*sig_np*1e-31*pow(Hbar*C, 4) * (m*m/(eng*eng));
}

double Dark_Photon::Brem_tot(double eng, double temp, double dens){
	return 2*Brem_t(eng, temp, dens) + Brem_l(eng, temp, dens);
}

double Dark_Photon::Brem_avg(double eng, double temp, double dens){
	return Brem_tot(eng, temp, dens)/3;
}

//Semi-Compton
double Dark_Photon::Comp_t(double eng, double temp, double dens){
	if(eng <= m) return 0;
	double Eps;
	if(eps == 1) Eps = 1;
	else Eps = Eps_t(eng, temp, dens);
	double n_e = (Y*dens/M_n)*pow(Hbar*C, 3);
	double E_F = sqrt(M_e*M_e + pow(3*PI*PI*n_e, 2.0/3.0));
	double Ome_pl = sqrt(4*PI*Alpha*n_e/E_F);
	return ((8*PI*Alpha*Alpha*Eps*Eps*n_e)/(3*E_F*E_F))*sqrt(Ome_pl/eng);
}

double Dark_Photon::Comp_l(double eng, double temp, double dens){
	if(eng <= m) return 0;
	double Eps;
	if(eps == 1) Eps = 1;
	else Eps = Eps_l(eng, temp, dens);
	double n_e = (Y*dens/M_n)*pow(Hbar*C, 3);
	double E_F = sqrt(M_e*M_e + pow(3*PI*PI*n_e, 2.0/3.0));
	double Ome_pl = sqrt(4*PI*Alpha*n_e/E_F);
	return ((8*PI*Alpha*Alpha*Eps*Eps*n_e)/(3*E_F*E_F))*sqrt(Ome_pl/eng)*(m*m)/(eng*eng);
}

double Dark_Photon::Comp_tot(double eng, double temp, double dens){
	return 2*Comp_t(eng, temp, dens) + Comp_l(eng, temp, dens);
}

double Dark_Photon::Comp_avg(double eng, double temp, double dens){
	return Comp_tot(eng, temp, dens)/3;
}

//Decay
double Dark_Photon::Decay_t(double eng, double temp, double dens){
	if(m <= 2*M_e || eng <= m) return 0;
	double Eps;
	if(eps == 1) Eps = 1;
	else Eps = Eps_t(eng, temp, dens);
	double pre = (Alpha*Eps*Eps*m*m)/sqrt(eng*eng-m*m);
	double n_e = (Y*dens/M_n)*pow(Hbar*C, 3);
	double mu_e = sqrt(M_e*M_e + pow(3*PI*PI*n_e, 2.0/3.0));
	double res, err;
	double xu = 0.5*(1+sqrt((1-(4*M_e*M_e)/(m*m))*(1-(m*m)/(eng*eng))));
	double xl = 0.5*(1-sqrt((1-(4*M_e*M_e)/(m*m))*(1-(m*m)/(eng*eng))));
	struct decay_params params = {eng, mu_e, temp, m};	

	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	gsl_function F;
	F.function = &Integrand_t;
	F.params = &params;
	gsl_integration_qag(&F, xl, xu, 0, 1e-2, 1000, GSL_INTEG_GAUSS61, w, &res, &err);
	gsl_integration_workspace_free(w);
	return pre * res;
}

double Dark_Photon::Decay_l(double eng, double temp, double dens){
	if(m <= 2*M_e || eng <= m) return 0;
	double Eps;
	if(eps == 1) Eps = 1;
	else Eps = Eps_l(eng, temp, dens);
	double pre = (Alpha*Eps*Eps*m*m)/sqrt(eng*eng-m*m);
	double n_e = (Y*dens/M_n)*pow(Hbar*C, 3);
	double mu_e = sqrt(M_e*M_e + pow(3*PI*PI*n_e, 2.0/3.0));
	double res, err;
	double xu = 0.5*(1+sqrt((1-(4*M_e*M_e)/(m*m))*(1-(m*m)/(eng*eng))));
	double xl = 0.5*(1-sqrt((1-(4*M_e*M_e)/(m*m))*(1-(m*m)/(eng*eng))));
	struct decay_params params = {eng, mu_e, temp, m};	

	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	gsl_function F;
	F.function = &Integrand_l;
	F.params = &params;
	gsl_integration_qag(&F, xl, xu, 0, 1e-2, 1000, GSL_INTEG_GAUSS61, w, &res, &err);
	gsl_integration_workspace_free(w);
	return pre * res;
}

double Dark_Photon::Decay_tot(double eng, double temp, double dens){
	return 2*Decay_t(eng, temp, dens) + Decay_l(eng, temp, dens);
}

double Dark_Photon::Decay_avg(double eng, double temp, double dens){
	return Decay_tot(eng, temp, dens)/3;
}

//Total rate of all processes
double Dark_Photon::Tot_rate(double eng, double temp, double dens){
	return Brem_tot(eng, temp, dens) + Comp_tot(eng, temp, dens) + Decay_tot(eng, temp, dens);
}

//Average(over polarizations) rate of all processes
double Dark_Photon::Avg_rate(double eng, double temp, double dens){
	return Brem_avg(eng, temp, dens) + Comp_avg(eng, temp, dens) + Decay_avg(eng, temp, dens);
}

//Effective mixing
double Dark_Photon::Eps_t(double eng, double temp, double dens){
	return eps*pow(pow(1-Re_PI_t(eng, temp, dens)/(m*m), 2) + pow(Im_PI_t(eng, temp, dens)/(m*m), 2) , -0.5);
//	return eps;
}
double Dark_Photon::Eps_l(double eng, double temp, double dens){
	return eps*pow(pow(1-Re_PI_l(eng, temp, dens)/(m*m), 2) + pow(Im_PI_l(eng, temp, dens)/(m*m), 2) , -0.5);
//	return eps;
}
double Dark_Photon::Re_PI_l(double eng, double temp, double dens){
	double n_e = (dens * Y / M_n)*pow(Hbar*C, 3);
	double E_F = sqrt(M_e*M_e + pow(3*PI*PI*n_e, 2.0/3.0));
	double omg_p2 = 4*PI*Alpha*n_e/E_F;
	double v = sqrt(1-(m*m)/(eng*eng));
	return (3*omg_p2/(v*v))*(1-v*v)*( ( 1 / (2*v) ) * log( (1+v) / (1-v) ) - 1 );
}
double Dark_Photon::Re_PI_t(double eng, double temp, double dens){
	double n_e = (dens * Y / M_n)*pow(Hbar*C, 3);
	double E_F = sqrt(M_e*M_e + pow(3*PI*PI*n_e, 2.0/3.0));
	double omg_p2 = 4*PI*Alpha*n_e/E_F;
	double v = sqrt(1-(m*m)/(eng*eng));
	return (3*omg_p2/(2*v*v))*( 1 - ( ( (1 - v*v) / (2*v) ) * log( (1+v) / (1-v) ) ) );
}

double Dark_Photon::Im_PI_t(double eng, double temp, double dens){
	Dark_Photon dark(m, 1);
	return -eng*( 1 - exp(-eng/temp) )*(dark.Decay_t(eng, temp, dens)+dark.Comp_t(eng, temp, dens)+dark.Brem_t(eng, temp, dens));
}

double Dark_Photon::Im_PI_l(double eng, double temp, double dens){
	Dark_Photon dark(m, 1);
	return -eng*( 1 - exp(-eng/temp) )*(dark.Decay_l(eng, temp, dens)+dark.Comp_l(eng, temp, dens)+dark.Brem_l(eng, temp, dens));
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////	Helper functions	//////////////////
//////////////////////////////////////////////////

//Integrand of the formula for e+e- decay
double Integrand_t(double x, void *p){
	struct decay_params *fp = (struct decay_params *)p;
	double num = Mat_ele2_t(fp->eng * x, fp->eng, fp->m);
	double den = 1 + exp((fp->mu_e - fp->eng * x)/fp->temp);
	return num/den;
}

double Mat_ele2_t(double x, double eng, double m){
	return M_e*M_e/(m*m) + Z(x, eng, m);
}

double Z(double x, double eng, double m){
	return ((x*(eng-x))/(m*m))-((x*eng)/(m*m)-0.5)*((eng*(eng-x))/(m*m)-0.5)/((eng*eng)/(m*m)-1);
}

double Integrand_l(double x, void *p){
	struct decay_params *fp = (struct decay_params *)p;
	double num = Mat_ele2_l(fp->eng * x, fp->eng, fp->m);
	double den = 1 + exp((fp->mu_e - fp->eng * x)/fp->temp);
	return num/den;
}

double Mat_ele2_l(double x, double eng, double m){
	return 1 - 2*Z(x, eng, m);
}


//Averaged np cross section for Bremsstrahlung
double Dark_Photon::sig_np_avg(double temp){
	double T[] = { 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
	double np[] = {0, 1.52961e-15, 2.84589e-05, 0.0543506, 2.03531, 16.3667, 62.0162, 154.306, 296.882, 483.044, 700.768, 937.018, 1180.25, 1421.41, 1653.93, 1873.41, 2077.36, 2264.48, 2434.39, 2587.42, 2724.33, 2845.85, 2953.45, 3047.94, 3130.49, 3202.26, 3264.15, 3317.37, 3362.53, 3400.74, 3432.54, 3518.77, 3404.93, 3235.27, 3058.25, 2890.1, 2735.57, 2595.26, 2468.16, 2353.02, 2248.35, 2152.85, 2065.46, 1985.13, 1911.07, 1842.46, 1778.86, 1719.6, 1295.69, 1035.01, 862.633, 739.403, 645.589, 572.587, 513.792, 465.344, 424.69, 390.107, 360.261, 334.268, 311.551, 291.429, 273.532, 257.527, 243.13, 230.161, 148.16, 108.805, 86.6119, 72.7898, 63.5985, 57.1985, 52.5823, 49.1493, 46.5151, 44.4248, 42.7016, 41.2248, 39.9101, 38.7005, 37.5587, 36.4593, 35.387, 34.3337, 34.3337, 24.6488, 17.2683, 12.2203, 8.84129, 6.55374, 4.97057, 3.84882, 3.03574};
	
	int i = 46, j = 0, k = 92;
	while(abs(k-j) > 1){
		if(temp < T[i]){
			int l = i;
			i = (i+j)/2;
			k = l;
		}
		else{
			int l = i;
			i = (i+k)/2;
			j = l;
		}
	}
	if((k-j) == 1)
		return np[j]+(np[k]-np[j])*(temp-T[j])/(T[k]-T[j]);
	else{
		cout<<"Exception! k-j = " << k-j << endl;
			return 0;
	}
}
