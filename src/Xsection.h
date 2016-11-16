#ifndef XSECTION_H
#define XSECTION_H

#include <cstdlib>
#include <vector>
#include <string>
#include <random>
#include "matrix_elements.h"
#include "sample_methods.h"

struct integrate_params{
	double (*dXdPS)(double * PS, size_t n_dims, void * params);
	double * params;
};

double gsl_1dfunc_wrapper(double x, void * params_);

//=============Xsection base class===================================================
// This is the base class for 2->2 and 2->3 cross-sections.
// It takes care of the tabulating details and the tabulating routines, also the interpolation process
// The actually total Xsection calcuate function and final state sample function are virtual functions, because 2->2 and 2->3 uses quite different techniques to do these jobs.
class Xsection{
protected:
	void tabulate_s_Temp(size_t T_start, size_t dnT);
	double (*dXdPS)(double * PS, size_t n_dims, void * params);
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
	Xsection(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double, double, double), double M1_, std::string name_);
	double get_M1(void) {return M1;};
	double interpX(double s, double Temp);
	virtual double calculate(double s, double Temp) = 0;
	virtual void sample_dXdPS(double s, double Temp, std::vector<std::vector<double> > & fs) = 0;
};

//============Derived 2->2 Xsection class============================================
class Xsection_2to2 : public Xsection{
private:
	rejection_1d sampler1d;
	std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist_phi3;
public:
    Xsection_2to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double, double, double), double M1_, std::string name_);
    double calculate(double s, double Temp);
	void sample_dXdPS(double s, double Temp, std::vector< std::vector<double> > & final_states);
};

//============Derived 2->3 Xsection class============================================
class Xsection_2to3 : public Xsection{
private:
	AiMS sampler;
	std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist_phi4;
	
public:
    Xsection_2to3(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double, double, double), double M1_, std::string name_);
    double calculate(double s, double Temp);
	void sample_dXdPS(double s, double Temp, std::vector< std::vector<double> > & final_states);
};

#endif
