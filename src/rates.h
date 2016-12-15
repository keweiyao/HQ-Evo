#ifndef RATE_H
#define RATE_H

#include <iostream>
#include <random>
#include <functional>
#include <vector>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include "constants.h"
#include "Xsection.h"

double inline f_0(double x, double xi);
double fy_wrapper22(double y, void * params_);
double fx_wrapper22(double x, void * px_);
double fy_wrapper23(double y, void * params_);
double fx_wrapper23(double x, void * px_);

struct integrate1d_params{
	std::function<double(double*)> interp;   // store a call to member funtion (interp) and object
	double * params;  // other double type parameters
};

double Vegas_func_wrapper(double * var, long unsigned int n_dims, void *params);

class rates{
protected:
	std::random_device rd;
    std::mt19937 gen;
    std::gamma_distribution<double> dist_x;
	std::uniform_real_distribution<double> dist_norm_y;
	std::uniform_real_distribution<double> dist_reject;
	virtual void tabulate_E1_T(size_t T_start, size_t dnT) = 0;
public:
	rates(std::string name_);
	virtual double calculate(double * arg) = 0;
	virtual double interpR(double * arg) = 0;
	virtual void sample_initial(double * arg_in, double * arg_out) = 0;
};

class rates_2to2 : public rates{
private:
	Xsection_2to2 * Xprocess;
	const double M;
	const int degeneracy;
	const size_t NE1, NT;
	const double E1L, E1H, TL, TH;
	const double dE1, dT;
	std::vector< std::vector<double> > Rtab;
	void tabulate_E1_T(size_t T_start, size_t dnT);
public:
	rates_2to2(Xsection_2to2 * Xprocess_, int degeneracy_, std::string name_);
	double calculate(double * arg);
	double interpR(double * arg);
	void sample_initial(double * arg_in, double * arg_out);
};

class rates_2to3 : public rates{
private:
	Xsection_2to3 * Xprocess;
	const double M;
	const int degeneracy;
	const size_t NE1, NT, Ndt;
	const double E1L, E1H, TL, TH, dtL, dtH;
	const double dE1, dT, ddt;
	std::vector< std::vector< std::vector<double> > > Rtab;
	void tabulate_E1_T(size_t T_start, size_t dnT);
public:
	rates_2to3(Xsection_2to3 * Xprocess_, int degeneracy_, std::string name_);
	double calculate(double * arg);
	double interpR(double * arg);
	void sample_initial(double * arg_in, double * arg_out);
};



#endif
