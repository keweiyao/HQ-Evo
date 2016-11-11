#include "matrix_elements.h"
#include "constants.h"
#include <cmath>

//=============running coupling=================================================
double inline alpha_s(double Q2){
    if (Q2 < Q2cut_l)
        return alpha0 / std::log( -Q2/Lambda2 );
    else if (Q2 <= Q2cut_h)
        return 1.0;
    else
        return alpha0 * ( .5 - std::atan( std::log(Q2/Lambda2)/M_PI ) / M_PI );
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
	double mt2 = 0.2*At*pf_g*T2;
	double result = c64d9pi2*At*At*(Q2u*Q2u + Q2s*Q2s + 2*M2*Q2t)/std::pow(Q2t - mt2, 2);
	return result;
}

double dX_Qq2Qq_dPS(double * PS, size_t n_dims, void * params){
	(void)n_dims;
	// unpacking parameters
	double t = PS[0];
	double * p = static_cast<double*>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qq2Qq(t, params)/c16pi/std::pow(s-M2, 2);
}

double approx_XQq2Qq(double s, double Temp, double M){
	(void)s;
	(void)Temp;
	(void)M;
	return 1.0;
}

//=============Baisc function for Q+g --> Q+g==================================
double M2_Qg2Qg(double t, void * params) {
	// unpacking parameters
	double * p = static_cast<double *>(params);
	double s = p[0], T2 = p[1]*p[1], M2 = p[2]*p[2];
	// define energy scales for each channel
	double Q2s = s - M2, Q2t = t, Q2u = M2 - s - t;
	// define coupling constant for each channel
	double As = alpha_s(Q2s), At = alpha_s(Q2t), Au = alpha_s(Q2u);
	// define Deybe mass for each channel
	double mt2 = 0.2*At*pf_g*T2, mu2 = Au*pf_q*T2, ms2 = As*pf_q*T2;
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
	return result*c16pi2;
}
double dX_Qg2Qg_dPS(double * PS, size_t n_dims, void * params){
	(void)n_dims;
	double t = PS[0];
	double * p = static_cast<double *>(params);
	double s = p[0], M2 = p[2]*p[2];
	return M2_Qg2Qg(t, params)/c16pi/std::pow(s-M2, 2);	
}
double approx_XQg2Qg(double s, double Temp, double M){	
	double Q2s = s-M*M, T2 = Temp*Temp, M2 = M*M;
	double abstmax = Q2s*Q2s/s;
	return 1./T2 - 1./(T2+abstmax) + 10.*M2/Q2s/Q2s;
}

//=============Basic for 2->2 + 1->2===========================================
// x: k, p4, phi4k, cos4 // both sin4 and sin* > 0
double M2_Qq2Qqg(double * x_, size_t n_dims_, void * params_){
	(void) n_dims_;
	// unpack parameters
	double * params = static_cast<double*>(params_);
	double s = params[0], sqrts = std::sqrt(params[0]);
	double T2 = params[1]*params[1];
	double M2 = params[2]*params[2];
	// unpack variables
	double k = 0.5*(x_[0]+x_[1]), p4 = 0.5*(x_[0]-x_[1]), phi4k = x_[2], cos4 = x_[3];
	double cos_star = ((s-M2)-2.*sqrts*(p4+k))/(2.*p4*k) +1.;
	// check integration range	
	if ( (p4+k) > sqrts || cos_star <= -1. || 1. <= cos_star ) return 0.0;
	// more useful variables
	double sin_star = std::sqrt(1. - cos_star*cos_star), sin4 = std::sqrt(1. - cos4*cos4);
	double cos_4k = std::cos(phi4k), sin_4k = std::sin(phi4k);
	// k-vec	
	double kx = k*(sin_star*cos_4k*cos4 + sin4*cos_star), 
		   ky = k*sin_star*sin_4k,
		   kz = k*(-sin_star*cos_4k*sin4 + cos4*cos_star);
	double kt2 = kx*kx + ky*ky;
	double x = (k+kz)/sqrts, xbar = (k+std::abs(kz))/sqrts;
	// q-perp-vec
	double qx = -p4*sin4;

	// 2->2
	double t = -(sqrts - M2/sqrts)*p4*(1.+cos4);
	double the_M2_Qq2Qq = M2_Qq2Qq(t, params);
	// 1->2
	double alpha_rad = alpha_s(kt2);
	double iD1 = 1./(kt2 + x*x*M2 + alpha_rad *pf_g*T2), iD2 = 1./(kt2 + qx*qx - 2*qx*kx  + x*x*M2 + alpha_rad *pf_g*T2);
	double Pg = alpha_rad*std::pow(1.-xbar, 2)*
			( kt2*std::pow(iD1-iD2, 2) + std::pow(qx*iD2, 2) - 2.*kx*qx*(iD1-iD2)*iD2 );

	// 2->3 = 2->2 * 1->2
	return 48*M_PI*the_M2_Qq2Qq*Pg;
}

double approx_XQq2Qqg(double s, double Temp, double M){
	(void)s;
	(void)Temp;
	(void)M;
	return 1.0;
}

double M2_Qg2Qgg(double * x_, size_t n_dims_, void * params_){
	(void) n_dims_;
	// unpack parameters
	double * params = static_cast<double*>(params_);
	double s = params[0], sqrts = std::sqrt(params[0]);
	double T2 = params[1]*params[1];
	double M2 = params[2]*params[2];
	// unpack variables
	double k = 0.5*(x_[0]+x_[1]), p4 = 0.5*(x_[0]-x_[1]), phi4k = x_[2], cos4 = x_[3];
	double cos_star = ((s-M2)-2.*sqrts*(p4+k))/(2.*p4*k) +1.;
	// check integration range	
	if ( (p4+k) > sqrts || cos_star <= -1. || 1. <= cos_star ) return 0.0;
	// more useful variables
	double sin_star = std::sqrt(1. - cos_star*cos_star), sin4 = std::sqrt(1. - cos4*cos4);
	double cos_4k = std::cos(phi4k), sin_4k = std::sin(phi4k);
	// k-vec	
	double kx = k*(sin_star*cos_4k*cos4 + sin4*cos_star), 
		   ky = k*sin_star*sin_4k,
		   kz = k*(-sin_star*cos_4k*sin4 + cos4*cos_star);
	double kt2 = kx*kx + ky*ky;
	double x = (k+kz)/sqrts, xbar = (k+std::abs(kz))/sqrts;
	// q-perp-vec
	double qx = -p4*sin4;

	// 2->2
	double t = -(sqrts - M2/sqrts)*p4*(1.+cos4);
	double the_M2_Qg2Qg = M2_Qg2Qg(t, params);
	// 1->2
	double alpha_rad = alpha_s(kt2);
	double iD1 = 1./(kt2 + x*x*M2 + alpha_rad *pf_g*T2), iD2 = 1./(kt2 + qx*qx - 2*qx*kx  + x*x*M2 + alpha_rad *pf_g*T2);
	double Pg = alpha_rad*std::pow(1.-xbar, 2)*
			( kt2*std::pow(iD1-iD2, 2) + std::pow(qx*iD2, 2) - 2.*kx*qx*(iD1-iD2)*iD2 );

	// 2->3 = 2->2 * 1->2
	return 48*M_PI*the_M2_Qg2Qg*Pg;
}

double approx_XQg2Qgg(double s, double Temp, double M){
	(void)s;
	(void)Temp;
	(void)M;
	return 1.0;
}

