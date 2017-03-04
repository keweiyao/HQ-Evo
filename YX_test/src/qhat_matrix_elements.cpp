#include <cmath>
#include <iostream>

#include "utility.h"
#include "matrix_elements.h"


//=============Baisc function for Q+q --> Q+q==================================

double dqhat_Qq2Qq_dPS(double *PS, size_t n_dims, void *params)
{
        (void) n_dims;
        double t = PS[0];
        double* p = static_cast<double*>(params);
        // for 2->2 process p = {s, temp, M, index};
        double s = p[0], M2=p[2]*p[2];
        int index = int(p[3]+0.5); // floor double p[3] into integer

        double p10_CoM = (s + M2)/(2*std::sqrt(s));
        double p20_CoM = (s - M2)/(2*std::sqrt(s));
        double cos_theta13 = 1. + t/(2*p20_CoM*p20_CoM);
        double sin_theta13 = std::sqrt(1-cos_theta13*cos_theta13);

        double q1_CoM = -p20_CoM * sin_theta13;
        double q3_CoM = p20_CoM * (1 - cos_theta13);
        double q11_CoM = q1_CoM*q1_CoM;
        double q33_CoM = q3_CoM*q3_CoM;
        double q13_CoM = q1_CoM*q3_CoM;

        double vecq[] = {q1_CoM, q3_CoM, q11_CoM, q33_CoM, q13_CoM};

        return M2_Qq2Qq(t, params)/c16pi/std::pow(s-M2,2) * vecq[index] ;
}




//============Basic function for Q+g->Q+g
double dqhat_Qg2Qg_dPS(double * PS, size_t n_dims, void * params)
{
	(void)n_dims;
	double t = PS[0];
	double * p = static_cast<double *>(params);
	double s = p[0], M2 = p[2]*p[2];
        int index = int(p[3]+0.5); // floor double p[3] into integer

        double p10_CoM = (s + M2)/(2*std::sqrt(s));
        double p20_CoM = (s - M2)/(2*std::sqrt(s));
        double cos_theta13 = 1. + t/(2*p20_CoM*p20_CoM);
        double sin_theta13 = std::sqrt(1-cos_theta13*cos_theta13);

        double q1_CoM = -p20_CoM * sin_theta13;
        double q3_CoM = p20_CoM * (1 - cos_theta13);
        double q11_CoM = q1_CoM*q1_CoM;
        double q33_CoM = q3_CoM*q3_CoM;
        double q13_CoM = q1_CoM*q3_CoM;

        double vecq[] = {q1_CoM, q3_CoM, q11_CoM, q33_CoM, q13_CoM};
	return M2_Qg2Qg(t, params)/c16pi/std::pow(s-M2, 2) * vecq[index];	
}




/*

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
	
	double x = (k+kz)/sqrts, xbar = (k+std::abs(kz))/sqrts;
	double tauk = k/(kt2+x*x*M2);

	double u = dt/tauk;
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
	double iD1 = 1./(kt2 + x2M2 + Lambda2), iD2 = 1./(kt2 - 2.*qx*kx  + qx2Mm);
	double Pg = alpha_rad*std::pow(1.-xbar, 2) 
				*LPM
				//*sudakov	
				*( (qx2Mm+x2M2)*iD1*iD2 - x2M2*iD1*iD1 - (x2M2 + mD2)*iD2*iD2 );

	// 2->3 = 2->2 * 1->2
	return c48pi*the_M2_Qq2Qq*Pg;
}

double approx_XQq2Qqg(double * arg, double M){
	double Temp = arg[1], dt = arg[2];
	(void)M;	
	return std::pow(dt/Temp, 2);
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
	double x = (k+kz)/sqrts, xbar = (k+std::abs(kz))/sqrts;
	double tauk = k/(kt2+x*x*M2);

	double u = dt/tauk;
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
	double iD1 = 1./(kt2 + x2M2 + Lambda2), iD2 = 1./(kt2 - 2.*qx*kx  + qx2Mm);
	double Pg = alpha_rad*std::pow(1.-xbar, 2)
				*LPM
				//*sudakov
				*( (qx2Mm+x2M2)*iD1*iD2 - x2M2*iD1*iD1 - (x2M2 + mD2)*iD2*iD2 );

	// 2->3 = 2->2 * 1->2
	return c48pi*the_M2_Qg2Qg*Pg;
}

double approx_XQg2Qgg(double * arg, double M){
	double Temp = arg[1], dt = arg[2];
	(void)M;	
	return std::pow(dt/Temp, 2);
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
	double D1 = kt2 + x2M2 + Lambda2;
	double D2 = kt2 + qt2 + qxkx*2. + x2M2 + mD2;
	double Pg = kt2/D1/D1 + (kt2 + qt2 + qxkx*2.)/D2/D2 - 2.*(kt2 + qxkx)/D1/D2;

	// 2->3 = 2->2 * 1->2
	return the_M2*Pg/16.;
}

double approx_XQqg2Qq(double * arg, double M){
	double s = arg[0], Temp = arg[1], a1 = arg[2], a2 = arg[3];
	double sqrts = std::sqrt(s);
	double xk = 0.5*(a1*a2 + a1 - a2);
	double x2 = 0.5*(-a1*a2 + a1 + a2);
	double M2 = M*M;
	double A = (2.*xk/x2 - 1.), B = -2.*sqrts*(1. + xk/x2), C = s - M2;
	double E2 = (-B - std::sqrt(B*B-4.*A*C))/2./A;
	double k = xk/x2*E2;
	double p1 = (1. - x2 - xk)/x2*E2;
	double E1 = std::sqrt(p1*p1 + M2);
	double cosk = (E2*E2-k*k-p1*p1)/2./p1/k;
	double kz = k*cosk;
	double kt2 = k*k - kz*kz;
	double frac = (k + kz)/(E1 + p1);
	double x2M2 = frac*frac*M2;
	double mD2t = alpha_s(kt2)*pf_g*Temp*Temp;
	double D1 = kt2 + x2M2 + Lambda2;
	double D2 = kt2 + x2M2 + mD2t;
	double prop2 = kt2/D1/D1 + mD2t/D2/D2;
	return (s - M*M)/Temp/Temp/x2*prop2;
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
	double D1 = kt2 + x2M2 + Lambda2;
	double D2 = kt2 + qt2 + qxkx*2. + x2M2 + mD2;
	double Pg = kt2/D1/D1 + (kt2 + qt2 + qxkx*2.)/D2/D2 - 2.*(kt2 + qxkx)/D1/D2;

	// 2->3 = 2->2 * 1->2
	return the_M2*Pg/16.;
}

double approx_XQgg2Qg(double * arg, double M){
	double s = arg[0], Temp = arg[1], a1 = arg[2], a2 = arg[3];
	double sqrts = std::sqrt(s);
	double xk = 0.5*(a1*a2 + a1 - a2);
	double x2 = 0.5*(-a1*a2 + a1 + a2);
	double M2 = M*M;
	double A = (2.*xk/x2 - 1.), B = -2.*sqrts*(1. + xk/x2), C = s - M2;
	double E2 = (-B - std::sqrt(B*B-4.*A*C))/2./A;
	double k = xk/x2*E2;
	double p1 = (1. - x2 - xk)/x2*E2;
	double E1 = std::sqrt(p1*p1 + M2);
	double cosk = (E2*E2-k*k-p1*p1)/2./p1/k;
	double kz = k*cosk;
	double kt2 = k*k - kz*kz;
	double frac = (k + kz)/(E1 + p1);
	double x2M2 = frac*frac*M2;
	double mD2t = alpha_s(kt2)*pf_g*Temp*Temp;
	double D1 = kt2 + x2M2 + Lambda2;
	double D2 = kt2 + x2M2 + mD2t;
	double prop2 = kt2/D1/D1 + mD2t/D2/D2;
	return (s - M*M)/Temp/Temp/x2*prop2;
}

*/
