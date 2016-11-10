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
