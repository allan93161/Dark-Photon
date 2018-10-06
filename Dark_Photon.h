//Definition of class Dark_Photon
#ifndef DARK_PHOTON_H
#define DARK_PHOTON_H
class Dark_Photon{
	public:
		Dark_Photon(double m_, double ep){m = m_; eps = ep;};
		double Brem_l(double eng, double temp, double dens);
		double Brem_t(double eng, double temp, double dens);
		double Brem_tot(double eng, double temp, double dens);
		double Brem_avg(double eng, double temp, double dens);
		double Comp_l(double eng, double temp, double dens);
		double Comp_t(double eng, double temp, double dens);
		double Comp_tot(double eng, double temp, double dens);
		double Comp_avg(double eng, double temp, double dens);
		double Decay_l(double eng, double temp, double dens);
		double Decay_t(double eng, double temp, double dens);
		double Decay_tot(double eng, double temp, double dens);
		double Decay_avg(double eng, double temp, double dens);
		double Tot_rate(double eng, double temp, double dens);
		double Avg_rate(double eng, double temp, double dens);
		double Eps_t(double eng, double temp, double dens);
		double Eps_l(double eng, double temp, double dens);
		double Re_PI_t(double eng, double temp, double dens);
		double Re_PI_l(double eng, double temp, double dens);
		double Im_PI_t(double eng, double temp, double dens);
		double Im_PI_l(double eng, double temp, double dens);

	private:
		double m;
		double eps;
		double sig_np_avg(double temp);

};

//Helper functions
//double sig_intg(double x, void* params);
//double sig_np(double x);
double Integrand_t(double x, void *p);
double Integrand_l(double x, void *p);
double Mat_ele2_t(double x, double eng, double m);
double Mat_ele2_l(double x, double eng, double m);
double Z(double x, double eng, double m);


#endif
