#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include "Dark_Photon.h"
#include "Constants.h"
#include "utility.h"
using namespace std;

//Calculate ejecta kinetic energy
double KE_ejecta(double m_, double eps, double r[], double dens[], double temp[], int N, double E_max, double GraBind[], double rad[]){
	const int n_points = 50;
	double EDark[n_points], Extra[n_points];
	for(int j = 0; j < n_points; j++){
		EDark[j] = E_deposit(rad[j], m_, eps, r, dens, temp, N) - E_max;
		Extra[j] = EDark[j] - GraBind[j];
	}
	
	double rad_M = 0, Extra_M = 0;
	FindMax(rad, Extra, rad_M, Extra_M, n_points);
	cout << "Ejecta Kinetic Energy: " << Extra_M << "erg" << endl;
	return Extra_M;
}

//Find maximum for KE_ejecta
void FindMax(double x[], double f[], double& x_max, double& f_max, int n){
	int i_max = 0;	
	for(int i = 0; i < n; i++){
		if(f[i] > f[i_max]) i_max = i;
	}
	x_max = x[i_max];
	f_max = f[i_max];
}

//Energy deposited by dark photons outside R
double E_deposit(double R, double m, double eps, double r[], double dens[], double temp[], int N){
	double time = 10;
	double E_dark = time*Lum(m, eps, R, r, dens, temp, N);
	return E_dark;
}

//Critical radius outside of which there's no resonance
double R_cr(double eng_plasma){
	double R[501], Eng[501];
	ifstream input;
	input.open("Plasma_frequency.txt");
	for(int i = 0; i <= 500; i++){
		input >> R[i] >> Eng[i];
	}

	int i = 250, j = 0, k = 500;
	while(abs(k-j) > 1){
		if(eng_plasma > Eng[i]){
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
		return R[j]+(R[k]-R[j])*(eng_plasma-Eng[j])/(Eng[k]-Eng[j]);
	else{
		cout<<"Exception! k-j = " << k-j << endl;
			return 0;
	}
}


//Luminosity of dark photons as a function of m_ and eps
//R_crit: Largest resonant radius possible
//For r < R_crit + 10e3, use delta function approx.
//For r > R_crit + 10e3, integrate as usual.

//		R_crit + 10e3 > R_n	R_crit + 10e3 < R_n	R_crit = 0
//Non-res	0			1			1
//Res		1			1			0
double Lum(double m_, double ep, double R, double r[], double dens[], double temp[], int N){
	struct my_params params = {m_, ep, R, r, dens, temp, N};
	double R_crit;
	if(m_ > 16.7) R_crit = 0;
	else R_crit = R_cr( m_*sqrt(2.0/3.0) );
//	cout << "R_n = " << R_n << "\tR_cr = " << R_crit << "\tm' = " << m_ << endl;

	//Non-resonant regime
	double res = 0, err = 0;
	if(R_crit+10e3 < R_n || R_crit == 0){
		gsl_monte_function F;
		F.f = &LumInt;
		F.dim = 2;
		F.params = &params;		
		
		double xl[2];
		if(R_crit != 0){
			xl[0] = m_;
			xl[1] = R_crit+10e3;
		}else{
			xl[0] = m_;
			xl[1] = R_crit;	
		}
		double xu[2] = {1e4, R_n};
		const gsl_rng_type *T;
		gsl_rng *rr;
		size_t calls = 5000;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		rr = gsl_rng_alloc(T);
		
		gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);
		gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, rr, s, &res, &err);
		gsl_monte_vegas_free(s);
	}
	cout << "err = " << err/624151 << "\tres = " << res/624151 << endl; 

	//Resonant regime
	double res_ = 0, err_ = 0;
	if(R_crit != 0){
		double xu_;
		if(R_crit + 10e3 < R_n)
			xu_ = R_crit + 10e3;
		else
			xu_ = R_n;
		double xl_ = 0;

		gsl_integration_workspace* w = gsl_integration_workspace_alloc(5000);
		gsl_function G;
		G.function = &LumInt_approx;
		G.params = &params;
		gsl_integration_qag(&G, xl_, xu_, 0, 5e-2, 5000, GSL_INTEG_GAUSS61, w, &res_, &err_);
		gsl_integration_workspace_free(w);
		cout << "err_ = " << err_/624151 << "\tres_ = " << res_/624151 << endl;
	}	
	//print result
//	cout << "Luminosity = " << (res + res_)/624151 << " erg/s" << endl;
	return (res + res_)/624151;
}

//Delta-function approximated differential luminosity ( 4pi*r^2*dL/dV )
double LumInt_approx(double x, void* p){
	struct my_params *fp = (struct my_params *)p;
	double eng_p = Omega_plasma(x, fp->r, fp->dens, fp->temp, fp->N);
	if(fp->m_ > sqrt(3.0/2.0)*eng_p) return 0;
	else{
		double eng_res = Omega_res(Omega_plasma(x, fp->r, fp->dens, fp->temp, fp->N), fp->m_);
		double tau = Tau(eng_res, x, fp->m_, fp->ep, fp->R, fp->r, fp->dens, fp->temp, fp->N);
		double a = pow(fp->m_*fp->ep, 2) * pow((eng_res*eng_res)-(fp->m_*fp->m_), 3.0/2.0);
		double b = 2*PI*(exp(eng_res/Temp(x, fp->r, fp->temp, fp->N))-1);
		double c = 2 + (fp->m_*fp->m_ - 3*eng_p*eng_p)/(eng_res*eng_res);
		if(fp->m_ < eng_p) return exp(-tau)*4*PI*x*x*a / (b*c*Hbar*pow(Hbar*C,3));
		else return exp(-tau)*4*PI*x*x*2*a/(b*c*Hbar*pow(Hbar*C,3));
	}
}

//Plasma frequency(#omega_{p}) as a function of radius
double Omega_plasma(double rad, double r[], double dens[], double temp[], int N){
	double n_e = (Y*Dens(rad, r, dens, N)/M_n) * pow(Hbar*C, 3);
	double E_F = sqrt(M_e*M_e + pow(3*PI*PI*n_e, 2.0/3.0) );
	return sqrt(4*PI*Alpha*n_e/E_F);
}

//Resonant dark photon energy(energy at which it hits a resonance)
//as a function of local plasma energy and dark photon mass
double Omega_res(double omega_plasma, double m){
	if(omega_plasma < m*sqrt(2.0/3.0) ) return 0;
	else{	
		double v_init = 0.001;	
		double v = v_init;
		double int_v = 0.1;
		struct res_params p = {omega_plasma, m};
		bool pos_i = (Res(v_init, &p) > 0); 
		bool pos_f = pos_i;
		v += int_v;
		for(int i = 0; i < 6; i++){
			while(pos_i == pos_f && v < 1 && v > 0){
				pos_i = pos_f;	
				pos_f = (Res(v, &p) > 0);
				v += int_v;
			}
			int_v *= -0.1;
			if(v >= 1) v = 1 - abs(int_v);
			if(v <= 0) v = abs(int_v);
			pos_i = (Res(v, &p) > 0); 
			pos_f = pos_i;
		}
		return m/sqrt(1-v*v);
	}
}

//Evaluate Re(PI)-m^2
double Res(double x/*resonant v*/, void *p){
	struct res_params *fp = (struct res_params *)p;
	double answer = 0;
	if(fp->m < fp->w /*&& x > 0 && x < 1*/){//longitudinal
		answer = ( 1 / (2*x) )*log( (1+x) / (1-x) ) - 1;
		answer *= ((3*fp->w*fp->w)/(x*x))*(1-x*x);
		answer = answer - fp->m*fp->m;
		return answer;
	}else if (fp->w < fp->m && fp->w > sqrt(2.0/3.0)*fp->m /*&& x > 0 && x < 1*/){//transverse
		answer = ( (1 - x*x) / (2*x) ) * log( (1+x) / (1-x) );
		answer = ((3*fp->w*fp->w)/(2*x*x))*(1 - answer);
		answer = answer - (fp->m*fp->m);
		return answer;
	}
	else return 0;
}	

//The differential luminosity(4pi*r^2*dL/dVdw) as a function of x[0] = energy, x[1] = radius
//For gsl integration
double LumInt(double *x, size_t dim, void *p){
	struct my_params *fp = (struct my_params *)p;
	return exp(-Tau(x[0],x[1],fp->m_, fp->ep, fp->R, fp->r, fp->dens, fp->temp, fp->N))*DfrPow(x[0],x[1],fp->m_, fp->ep, fp->r, fp->dens, fp->temp, fp->N)*4*PI*x[1]*x[1];
}

//The differential luminosity (4pi*r^2*dL/dVdw)
//For general use
double LumInt_general(double eng, double rad, double m_, double ep, double R /* R_far */, double r[], double dens[], double temp[], double N){
	return exp(-Tau(eng, rad, m_, ep, R, r, dens, temp, N))*DfrPow(eng, rad, m_, ep, r, dens, temp, N)*4*PI*rad*rad;
}

//The differential power (4pi*r^2*dP/dVdw)
//For general use
double PowInt_general(double eng, double rad, double m_, double ep, double r[], double dens[], double temp[], double N){
	return DfrPow(eng, rad, m_, ep, r, dens, temp, N)*4*PI*rad*rad;
}

//Optical depth
double Tau(double eng, double rad, double m_, double ep, double R, double r[], double dens[], double temp[], int N){
	double pre = 1 - rad*(rad-R_c)/(2*R_n*R_n);

//This is for scheme 1
/*	Dark_Photon dark(m_, ep);
	double tau_far = (R-R_n)*(dark.Decay_avg(eng, Temp(R_n, r, temp, N), Dens(R_n, r, dens, N)) + dark.Comp_avg(eng, Temp(R_n, r, temp, N), Dens(R_n, r, dens, N)) + 0.1*dark.Brem_avg(eng, Temp(R_n, r, temp, N), Dens(R_n, r, dens, N)));


	double res_, err_;
	gsl_integration_workspace* w_ = gsl_integration_workspace_alloc(5000);
	gsl_function F_;
	F_.function = &Tau_int;
	struct my_params params_ = {m_, ep, eng, r, dens, temp, N};
	F_.params = &params_;
	gsl_integration_qag(&F_, rad, R_n, 0, 1e-3, 5000, GSL_INTEG_GAUSS15, w_, &res_, &err_);
	gsl_integration_workspace_free(w_);
*/

//This is for scheme 2
	double res, err;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(5000);
	gsl_function F;
	F.function = &Tau_int;
	struct my_params params = {m_, ep, eng, r, dens, temp, N};
	F.params = &params;
	gsl_integration_qag(&F, rad, R_far, 0, 1e-3, 5000, GSL_INTEG_GAUSS15, w, &res, &err);
	gsl_integration_workspace_free(w);
//	cout << "res = " << res/(Hbar*C) << "\terr = " << err/(Hbar*C) << "\tTau = " << (pre*res + tau_far)/(Hbar*C) << endl;
	return pre*res/(Hbar*C);
}

//Integrand of optical depth
double Tau_int(double x, void* p){
	struct my_params *fp = (struct my_params *)p;
	Dark_Photon dark(fp->m_, fp->ep);
	return dark.Avg_rate(fp->R, Temp(x, fp->r, fp->temp, fp->N), Dens(x, fp->r, fp->dens, fp->N));
}

//Differential power (dP/dVdw)
double DfrPow(double eng, double rad, double m_, double ep, double r[], double dens[], double temp[], int N){
	double pre = ( pow( eng , 3) * pow( 1 - m_*m_/(eng*eng) , 1.5) / (2*PI*PI) )*exp(-eng/Temp(rad, r, temp, N));
	Dark_Photon dark(m_, ep);
	return pre * (dark.Brem_tot(eng, Temp(rad, r, temp, N), Dens(rad, r, dens, N)) + dark.Comp_tot(eng, Temp(rad, r, temp, N), Dens(rad, r, dens, N)) ) /(Hbar*pow(Hbar*C, 3));
}

//Temperature profile
double Temp(double rad, double r[], double temp[], int N){
	int i = (N-1)/2, j = 0, k = N-1;
	if(rad*100 > 7e13) return temp[N-1];
	while(abs(k-j) > 1){
		if(rad*100 < r[i]){
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
	if((k-j) == 1){
		return temp[j]+(temp[k]-temp[j])*(100*rad-r[j])/(r[k]-r[j]);	
	}
	else{
		cout<<"Exception! k-j = " << k-j << endl;
		return 0;
	}
}

//Density profile
double Dens(double rad, double r[], double dens[], int N){
	int i = (N-1)/2, j = 0, k = N-1;
	if(rad*100 > 7e13) return dens[N-1]/pow(Hbar*C, 3);
	while(abs(k-j) > 1){
		if(rad*100 < r[i]){
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
	if((k-j) == 1){
		return ( dens[j]+(dens[k]-dens[j])*(100*rad-r[j])/(r[k]-r[j]) )/pow(Hbar*C, 3);
	}
	else{
		cout<<"Exception! k-j = " << k-j << endl;
		return 0;
	}
}


