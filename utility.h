#ifndef UTILITY_H
#define UTILITY_H

double KE_ejecta(double m_, double eps, double r[], double dens[], double temp[], int N, double E_max, double GraBind[], double rad[]);
void FindMax(double x[], double f[], double& x_max, double& f_max, int n);
double E_deposit(double R, double m, double eps, double r[], double dens[], double temp[], int N);
double R_cr(double eng_plasma);
double Lum(double m_, double ep, double R, double r[], double dens[], double temp[], int N);
double LumInt_approx(double x, void* p);
double Omega_plasma(double rad, double r[], double dens[], double temp[], int N);
double Omega_res(double omega_plasma, double m);
double Res(double x/*resonant v*/, void *p);
double LumInt(double *x, size_t dim, void *p);
double LumInt_general(double eng, double rad, double m_, double ep, double R /* R_far */, double r[], double dens[], double temp[], double N);
double PowInt_general(double eng, double rad, double m_, double ep, double r[], double dens[], double temp[], double N);
double Tau(double eng, double rad, double m_, double ep, double R, double r[], double dens[], double temp[], int N);
double Tau_int(double x, void* p);
double DfrPow(double eng, double rad, double m_, double ep, double r[], double dens[], double temp[], int N);
double Temp(double rad, double r[], double temp[], int N);
double Dens(double rad, double r[], double dens[], int N);

#endif
