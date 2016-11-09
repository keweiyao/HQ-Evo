#ifndef RATE_H
#define RATE_H

#include <iostream>
#include <random>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <functional>
#include <vector>

//=============constants=======================================================
const double c4d9 = 4./9.;
const double c1d9 = 1./9.;
const double c16pi = 16.*M_PI;
const double c16pi2 = 16.*M_PI*M_PI;
const double c64d9pi2 = 64./9.*M_PI*M_PI;
const int Nc = 3, nf = 3;
const double pf_g = 4.*M_PI/3.*(Nc + nf/2.); // prefractor for gluon self energy^2 
const double pf_q = M_PI/2.*(Nc*Nc - 1)/2./Nc; // prefractor for quark self energy^2 
const double alpha0 = 4.*M_PI/(11. - 2./3.*nf); // alpha_s(Q2 = e*Lambda2)
const double Lambda2 = 0.2*0.2; // [GeV^2] Lambda QCD squared
const double Q2cut_l = -Lambda2*exp(alpha0), Q2cut_h = 0.; // [GeV^2] ranges within which alphas > 1 and will be cut

//======running coupling=======================================================
double alpha_s(double Q2);

//=============Baisc function for Q+q --> Q+q==================================
double dX_Qq2Qq_dt(double t, void * params);
double approx_XQq2Qq(double s, double Temp, double M);
//=============Baisc function for Q+g --> Q+g==================================
double dX_Qg2Qg_dt(double t, void * params);
double approx_XQg2Qg(double s, double Temp, double M);

//=============================================================================
class Xsection_2to2{
private:
	void tabulate_st(void);
	double (*dXdt)(double t, void * params);
	double (*approx_X)(double s, double Temp, double M);
	double M1;
	// For interpolation
	// Cross-section changes abruptly near the threshhold s>M2,
	// Use different grid size for M2 < s < 5*M2 amd 5*M2 < s < 400*M2
	std::vector< std::vector<double> > Xtab;
	const size_t Ns, NT;
	const double sL, sM, sH, ds1, ds2;
	const double TL, TH, dT; 
public:
    Xsection_2to2(double (*dXdt_)(double t, void * params), double (*approx_X_)(double s, double Temp, double M), double M1);
    double calculate(double s, double Temp);
	double interpX(double s, double Temp);
	double sample_dXdt(double s, double Temp);
};

struct Vegas_params{
	std::function<double(double, double)> interp;   // store a call to member funtion (interp) and object
	double * params;  // other double type parameters
};

double Vegas_func_wrapper(double * var, long unsigned int n_dims, void *params);

template <class T>
class rates{
private:
	Vegas_params * params;
	T * Xprocess;
	void tabulate_E1_T(size_t T_start, size_t dnT);
	const int degeneracy;
	const size_t NE1, NT;
	std::vector< std::vector<double> > Rtab;
public:
	rates(T * Xprocess_, int degeneracy_);
	double calculate(double E1, double Temp);
	void sample_initial(double E1, double Temp, double &E2, double &s);
};



#endif
