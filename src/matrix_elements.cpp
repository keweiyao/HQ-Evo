#include <cmath>
#include <iostream>

#include "utility.h"
#include "matrix_elements.h"

#include <boost/math/tools/roots.hpp>



//=============running coupling=================================================
double alpha_s(double Q2){
    if (Q2 < Q2cut_l)
        return alpha0 / std::log( -Q2/Lambda2 );
    else if (Q2 <= Q2cut_h)
        return 1.0;
    else
        return alpha0 * ( .5 - std::atan( std::log(Q2/Lambda2)/M_PI ) / M_PI );
}

//=============Debye mass=================================================
Debye_mass * t_channel_mD2 = NULL;

Debye_mass::Debye_mass(const unsigned int _type, const double _mDTc,
						   const double _mDslope, const double _mDcurv, 
						   const double _Tc)
:	TL(0.1), TH(1.0), NT(100), dT((TH-TL)/(NT-1.)), 
	mDTc(_mDTc), mDslope(_mDslope), mDcurv(_mDcurv), Tc(_Tc),
	type(_type), mD2(new double[NT])
{
	
	if (type==0) {
		std::cout << "self-consistent Debye mass" << std::endl;
		// type==0 use self-consistent Debye mass
		for (size_t i=0; i<NT; i++){
			double T = TL+dT*i;
			size_t maxiter=100;
			boost::math::tools::eps_tolerance<double> tol{
	   		 (std::numeric_limits<double>::digits * 3) / 4};
			try{
				auto result = boost::math::tools::toms748_solve(
		 			[&T](double x) {return pf_g*alpha_s(-x)*T*T - x;},
		  			0.01, 20., tol, maxiter);
		 		mD2[i] = .5*(result.first + result.second);
			}
			catch (const std::domain_error&) {
				throw std::domain_error{
		  		"unable to calculate mD2"};
			}
		}
	}
	if (type==1) {
		std::cout << "parameterized Debye mass" << std::endl;
		double mD = 0.;
		// use parameterized Debye mass
		for (size_t i=0; i<NT; i++){
			double T = TL+dT*i;
			mD = mDTc*(1. + mDslope*(T-Tc)/Tc*std::pow(T/Tc, mDcurv) );
		 	mD2[i] = mD*mD;
		}
	}
}

double Debye_mass::get_mD2(double T){
	if (T<TL) T=TL;
	if (T>=TH-dT) T=TH-dT;
	double x = (T-TL)/dT;
	size_t index = std::floor(x);
	double r = x-index;
	return (1.-r)*mD2[index] + r*mD2[index+1];
}

void initialize_Debye_mass(const unsigned int type, const double mDTc,
						   const double mDslope, const double mDcurv, 
						   const double Tc){
	t_channel_mD2 = new Debye_mass(type, mDTc, mDslope, mDcurv, Tc);
}



//=============Soft regulator==================================================
double f_IR(double x){
	double a = std::sqrt(1. + 4.*x);
	return (1. + 2.*x)/a*std::log((a+1.)/(a-1.)) - 1.;
}

//=============Baisc function for Q+q --> Q+q==================================
double M2_Qq2Qq(double t, void * params){
	// unpacking parameters
	double * p = static_cast<double*>(params);
	double s = p[0], T2 = p[1]*p[1], M2 = p[2]*p[2];
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = t, Q2u = M2 - s - t;
	// define coupling constant for each channel
	double At = alpha_s(Q2t);
	// define Deybe mass for each channel
	double mt2 = 0.2*t_channel_mD2->get_mD2(p[1]);
	double result = c64d9pi2*At*At*(Q2u*Q2u + Q2s*Q2s + 2.*M2*Q2t)/std::pow(Q2t - mt2, 2);
	if (result < 0.) return 0.;
	else return result;
}

double dX_Qq2Qq_dPS(double * PS, size_t n_dims, void * params){
	(void)n_dims;
	// unpacking parameters
	double t = PS[0];
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qq2Qq(t, params)/c16pi/std::pow(s-M2, 2);
}

//=============Baisc function for Q+g --> Q+g==================================
double M2_Qg2Qg(double t, void * params) {
	// unpacking parameters
	double * p = static_cast<double *>(params);
	double s = p[0], T2 = p[1]*p[1], M2 = p[2]*p[2];
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = t, Q2u = M2 - s - t;
	// define coupling constant for each channel
	double At = alpha_s(Q2t), Au = alpha_s(Q2u), As = alpha_s(Q2s);
	// define Deybe mass for each channel
	double mt2 = 0.2*t_channel_mD2->get_mD2(p[1]), 
		   mu2 = Au*pf_q*T2, ms2 = As*pf_q*T2;
	double result = 0.0;
	// t*t
	result += 2.*At*At * Q2s*(-Q2u)/std::pow(Q2t - mt2, 2);
		
	// s*s
	result += c4d9*As*As * ( Q2s*(-Q2u) + 2.*M2*(Q2s + 2.*M2) ) / std::pow(Q2s + ms2, 2);
	// u*u
	result += c4d9*Au*Au * ( Q2s*(-Q2u) + 2.*M2*(Q2u + 2.*M2) ) / std::pow(-Q2u + mu2, 2);
	// s*u
	result += c1d9*As*Au * M2*(4.*M2 - Q2t) / (Q2s + ms2) / (-Q2u + mu2);
	// t*s
	result += At*As * ( Q2s*(-Q2u) + M2*(Q2s - Q2u) ) / (Q2t - mt2) / (Q2s + ms2);
    // t*u
	result += -At*Au * ( Q2s*(-Q2u) - M2*(Q2s - Q2u) ) / (Q2t - mt2) / (-Q2u + mu2);
	if (result < 0.) return 0.;
	return result*c16pi2;
}

double M2_Qg2Qg_only_t(double t, void * params) {
	// unpacking parameters
	double * p = static_cast<double *>(params);
	double s = p[0], T2 = p[1]*p[1], M2 = p[2]*p[2];
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = t, Q2u = M2 - s - t;
	// define coupling constant for each channel
	double At = alpha_s(Q2t);
	// define Deybe mass for each channel
	double mt2 = 0.2*t_channel_mD2->get_mD2(p[1]);
	double result = c16pi2*2.*At*At * Q2s*(-Q2u)/std::pow(Q2t - mt2, 2);
	if (result < 0.) return 0.;
	else return result;
}

double dX_Qg2Qg_dPS(double * PS, size_t n_dims, void * params){
	(void)n_dims;
	double t = PS[0];
	double * p = static_cast<double *>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qg2Qg(t, params)/c16pi/std::pow(s-M2, 2);	
}

//=============Basic for 2->3===========================================
// x: k, p4, phi4k, cos4 // both sin4 and sin* > 0
double M2_Qq2Qqg(double * x_, size_t n_dims_, void * params_){
	(void) n_dims_;
	// unpack parameters
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T2 = params[1]*params[1];
	double M = params[2];
	double M2 = M*M;
	double dt = params[3]; // separation time between this and the last scattering, in CoM frame [GeV-1]
	double pmax =  0.5*sqrts*(1.-M2/s);
	// unpack variables
	double k = 0.5*(x_[0]+x_[1]), p4 = 0.5*(x_[0]-x_[1]), phi4k = x_[2], cos4 = x_[3];
	double cos_star = ((s-M2)-2.*sqrts*(p4+k))/(2.*p4*k) +1.;
	// check integration range	
	if ( phi4k <= 0. || phi4k >= 2.*M_PI || cos4 <= -1. || cos4 >= 1.) return 0.0;
	if ( p4 <= 0. || k <= 0. || p4 >= pmax || k >= pmax ) return 0.0;
	if ( (p4+k) > sqrts || cos_star <= -1. || 1. <= cos_star) return 0.0;
	// more useful variables
	double sin_star = std::sqrt(1. - cos_star*cos_star), sin4 = std::sqrt(1. - cos4*cos4);
	double cos_4k = std::cos(phi4k), sin_4k = std::sin(phi4k);
	// k-vec	
	double kx = k*(sin_star*cos_4k*cos4 + sin4*cos_star), 
		   ky = k*sin_star*sin_4k,
		   kz = k*(-sin_star*cos_4k*sin4 + cos4*cos_star);
	double kt2 = kx*kx + ky*ky;
	
	double x = std::max((k+kz)/sqrts, min_xfrac), xbar = (k+std::abs(kz))/sqrts;
	double tauk = 2.*k/(kt2+x*x*M2);


	// here u is the ratio of the mean-free-path over the formation length
	// mean-free-path \sim mean-free-time*v_HQ, 
	// v_HQ = p/E = (s - M^2)/(s + M^2)
	// formation length = tau_k*v_k = tau_k
	double u = dt*(s-M2)/(s+M2)/tauk;
	double LPM = u*u/(1.+u*u);

	// q-perp-vec
	double qx = -p4*sin4;
	double alpha_rad = alpha_s(kt2);
	double mD2 = alpha_rad *pf_g*T2;

	double x2M2 = x*x*M2;
	double qx2Mm = qx*qx + x2M2 + mD2;
	
	// 2->2
	double t = -(sqrts - M2/sqrts)*p4*(1.+cos4);
	double the_M2_Qq2Qq = M2_Qq2Qq(t, params); 

	//double sudakov = std::exp(alpha_rad/M_PI*f_IR(-M2/t)*std::log(-k*k/t));
	// 1->2
	double iD1 = 1./(kt2 + x2M2), iD2 = 1./(kt2 - 2.*qx*kx  + qx2Mm);
	double Pg = alpha_rad*std::pow(1.-xbar, 2) 
				*LPM
				//*sudakov	
				*( (qx2Mm+x2M2)*iD1*iD2 - x2M2*iD1*iD1 - (x2M2 + mD2)*iD2*iD2 );

	// 2->3 = 2->2 * 1->2
	return c48pi*the_M2_Qq2Qq*Pg;
}



double M2_Qg2Qgg(double * x_, size_t n_dims_, void * params_){
	(void) n_dims_;
	// unpack parameters
	double * params = static_cast<double*>(params_);
	double s = params[0];
	double sqrts = std::sqrt(s);
	double T2 = params[1]*params[1];
	double M = params[2];
	double M2 = M*M;
	double dt = params[3]; // separation time between this and the last scattering, in CoM frame [GeV-1]
	double pmax =  0.5*sqrts*(1.-M2/s);
	// unpack variables
	double k = 0.5*(x_[0]+x_[1]), p4 = 0.5*(x_[0]-x_[1]), phi4k = x_[2], cos4 = x_[3];
	double cos_star = ((s-M2)-2.*sqrts*(p4+k))/(2.*p4*k) +1.;
	// check integration range	
	if ( phi4k <= 0. || phi4k >= 2.*M_PI || cos4 <= -1. || cos4 >= 1.) return 0.0;
	if ( p4 <= 0. || k <= 0. || p4 >= pmax || k >= pmax ) return 0.0;
	if ( (p4+k) > sqrts || cos_star <= -1. || 1. <= cos_star ) return 0.0;
	// more useful variables
	double sin_star = std::sqrt(1. - cos_star*cos_star), sin4 = std::sqrt(1. - cos4*cos4);
	double cos_4k = std::cos(phi4k), sin_4k = std::sin(phi4k);
	// k-vec	
	double kx = k*(sin_star*cos_4k*cos4 + sin4*cos_star), 
		   ky = k*sin_star*sin_4k,
		   kz = k*(-sin_star*cos_4k*sin4 + cos4*cos_star);
	double kt2 = kx*kx + ky*ky;
	double x = std::max((k+kz)/sqrts, min_xfrac), xbar = (k+std::abs(kz))/sqrts;
	double tauk = 2.*k/(kt2+x*x*M2);

	// here u is the ratio of the mean-free-path over the formation length
	// mean-free-path \sim mean-free-time*v_HQ, 
	// v_HQ = p/E = (s - M^2)/(s + M^2)
	// formation length = tau_k*v_k = tau_k
	double u = dt*(s-M2)/(s+M2)/tauk;
	double LPM = u*u/(1.+u*u);

	// q-perp-vec
	double qx = -p4*sin4;
	double alpha_rad = alpha_s(kt2);
	double mD2 = alpha_rad *pf_g*T2;

	double x2M2 = x*x*M2;
	double qx2Mm = qx*qx + x2M2 + mD2;

	// 2->2
	double t = -(sqrts - M2/sqrts)*p4*(1.+cos4);
	double the_M2_Qg2Qg = M2_Qg2Qg_only_t(t, params);

	//double sudakov = std::exp(alpha_rad/M_PI*f_IR(-M2/t)*std::log(-k*k/t));
	// 1->2
	double iD1 = 1./(kt2 + x2M2), iD2 = 1./(kt2 - 2.*qx*kx  + qx2Mm);
	double Pg = alpha_rad*std::pow(1.-xbar, 2)
				*LPM
				//*sudakov
				*( (qx2Mm+x2M2)*iD1*iD2 - x2M2*iD1*iD1 - (x2M2 + mD2)*iD2*iD2 );

	// 2->3 = 2->2 * 1->2
	return c48pi*the_M2_Qg2Qg*Pg;
}



//=============Basic for 3->2===========================================
double Ker_Qqg2Qq(double * x_, size_t n_dims_, void * params_){
	(void) n_dims_;
	// unpack variables costheta42 = x_[0]
	double costheta24 = x_[0], phi24 = x_[1];
	if (costheta24<=-1. || costheta24>=1. || phi24 <=0. || phi24 >=2.*M_PI) return 0.;
	double sintheta24 = std::sqrt(1. - costheta24*costheta24), sinphi24 = std::sin(phi24), cosphi24 = std::cos(phi24);
	// unpack parameters
	double * params = static_cast<double*>(params_); // s12k, T, M, 2*E2*E4
	double E2 = params[3];
	double E4 = params[4];
	double TwoE2E4 = 2.*E2*E4;
	double kt2 = params[5];
	double costheta2 = params[6];
	double sintheta2 = std::sqrt(1. - costheta2*costheta2);
	double x2M2 = params[7];
	double mD2 = params[8];

	// 2->2
	double t = TwoE2E4 * (costheta24 - 1);
	double the_M2 = M2_Qq2Qq(t, params);

	// 1->2
	double qx = E2*sintheta2 - E4*(costheta2*sintheta24*cosphi24 + sintheta2*costheta24),
		   qy = -E4*sintheta24*sinphi24;
	double qxkx = -std::sqrt(kt2)*qx;
	double qt2 = qx*qx + qy*qy;
	double D1 = kt2 + x2M2;
	double D2 = kt2 + qt2 + qxkx*2. + x2M2 + mD2;
	double Pg = kt2/D1/D1 + (kt2 + qt2 + qxkx*2.)/D2/D2 - 2.*(kt2 + qxkx)/D1/D2;

	// 2->3 = 2->2 * 1->2
	return the_M2*Pg/16.;
}

double Ker_Qgg2Qg(double * x_, size_t n_dims_, void * params_){
(void) n_dims_;
	// unpack variables costheta42 = x_[0]
	double costheta24 = x_[0], phi24 = x_[1];
	if (costheta24<=-1. || costheta24>=1. || phi24 <=0. || phi24 >=2.*M_PI) return 0.;
	double sintheta24 = std::sqrt(1. - costheta24*costheta24), sinphi24 = std::sin(phi24), cosphi24 = std::cos(phi24);
	// unpack parameters
	double * params = static_cast<double*>(params_); // s12k, T, M, 2*E2*E4
	double E2 = params[3];
	double E4 = params[4];
	double TwoE2E4 = 2.*E2*E4;
	double kt2 = params[5];
	double costheta2 = params[6];
	double sintheta2 = std::sqrt(1. - costheta2*costheta2);
	double x2M2 = params[7];
	double mD2 = params[8];

	// 2->2
	double t = TwoE2E4 * (costheta24 - 1);
	double the_M2 = M2_Qg2Qg_only_t(t, params);

	// 1->2
	double qx = E2*sintheta2 - E4*(costheta2*sintheta24*cosphi24 + sintheta2*costheta24),
		   qy = -E4*sintheta24*sinphi24;
	double qxkx = -std::sqrt(kt2)*qx;
	double qt2 = qx*qx + qy*qy;
	double D1 = kt2 + x2M2;
	double D2 = kt2 + qt2 + qxkx*2. + x2M2 + mD2;
	double Pg = kt2/D1/D1 + (kt2 + qt2 + qxkx*2.)/D2/D2 - 2.*(kt2 + qxkx)/D1/D2;

	// 2->3 = 2->2 * 1->2
	return the_M2*Pg/16.;
}
