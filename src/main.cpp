#include <iostream>
#include "Rates.h"
#include "sample_methods.h"
#include <fstream>

double test(double * x, size_t n_dims, void * params){
	double * p = static_cast<double*>(params);
	double result = 0.0;
	for (size_t i=0; i<n_dims; ++i){
		result += std::pow(x[i]/p[i], 2);
	}
	return std::exp(-0.5*result);
	//	return std::exp(-(100.*std::pow(x2-x1*x1, 2)+std::pow(1.-x1, 2))/20.);
}

// x: k, p4, phi4k, cos4 // both sin4 and sin* > 0
double M2_Qq2Qqg(double * x_, size_t n_dims_, void * params_){
	// unpack parameters
	double * params = static_cast<double*>(params_);
	double s = params[0], sqrts = std::sqrt(params[0]);
	double T2 = params[1]*params[1];
	double M2 = params[2]*params[2];
	// unpack variables
	double k = x_[0], p4 = x_[1], phi4k = x_[2], cos4 = x_[3];
	double cos_star = ((s-M2)-2.*sqrts*(p4+k))/(2.*p4*k) +1.;
	// check integration range	
	if (phi4k <= 0. || phi4k >= 2.*M_PI || cos4 <= -1. || cos4 >= 1.
		|| k <= 0. || p4 < 0. || (p4+k) > sqrts || cos_star <= -1. || cos_star >= 1.)return 0.0;
	// more useful variables
	double t = -0.5*(sqrts - M2/sqrts)*p4*(1.+cos4);
	double M2_Qq2Qq = dX_Qq2Qq_dt(t, params);
	double sin_star = std::sqrt(1. - cos_star*cos_star), sin4 = std::sqrt(1. - cos4*cos4);
	double cos_4k = std::cos(phi4k), sin_4k = std::sin(phi4k);
	double kx = k*(sin_star*cos_4k*cos4 - sin4*cos_star), ky = sin_star*sin_4k,
		   kz = k*(sin_star*cos_4k*sin4 + cos4*cos_star);
	double kt2 = kx*kx + ky*ky;
	double qx = -p4*sin4;
	double x = (k+kz)/sqrts, xbar = (k+std::abs(kz))/sqrts;
	double alpha_rad = alpha_s(kt2);
	
	double iD1 = 1./(kt2 + x*x*M2 + 0.2*alpha_rad *pf_g*T2), iD2 = 1./(kt2 + qx*qx - 2*qx*kx  + x*x*M2 + 0.2*alpha_rad *pf_g*T2);
	double Pg = std::pow(1.-xbar*xbar, 2)*
			( kt2*std::pow(iD1+iD2, 2) + std::pow(qx*iD2, 2) - 2.*kx*qx*(iD1+iD2)*iD2 );
	return 48.*M_PI*alpha_rad*Pg*M2_Qq2Qq;
}

int main(){
	//rejection_1d sampler1d;
	//double * p = new double[2]; p[0] = 0.1; p[1] = 0.1;
	//for(int i=0; i< 1000000; i++) std::cout << sampler1d.sample(test, 0., 1., p) << std::endl;

	//--------test MC sampler------------
	double * p = new double[3]; //s, T, M
	double M = 1.3;
	double sqrts = M*1.2;
	double s = sqrts*sqrts;
	p[0] = s; p[1] = 0.4;  p[2] = M;
	size_t n_dims = 4;
	double * guess = new double[n_dims];
	guess[0] = sqrts/2. - M/2.;
	guess[1] = guess[0];
	guess[2] = M_PI;
	guess[3] = 0.1;
	AiMS sampler;
	sampler.sample(M2_Qq2Qqg, n_dims, p, guess);	

	//Xsection_2to2 XQg2Qg(&dX_Qg2Qg_dt, &approx_XQg2Qg, 1.3);
	//Xsection_2to2 XQq2Qq(&dX_Qq2Qq_dt, &approx_XQq2Qq, 1.3);
	//rates<Xsection_2to2> RQg2Qg(&XQg2Qg, 8*2);
	//rates<Xsection_2to2> RQq2Qq(&XQq2Qq, 3*4);

	return 0;
}
