#ifndef XSECTION_H
#define XSECTION_H

#include <cstdlib>
#include <vector>
#include <string>
#include <random>
  
#include "sample_methods.h"

/* all the differential Xsection function are declared by type "double f(double * arg, size_t n_dims, void * params)"
	This type is natural for gsl vegas integration, but does not match the type for gsl 1d integration which requires double f(double x, void * params)
	Therefore, for 1d integration, use the struct below. It stores both the function call and the numerical parameters:
	
	double f_1d(double x, void * p_) {
		Mygsl_integration_params * p = static_cast<Mygsl_integration_params *> p;
		return p->f(x, 1, p->params)
	}
*/
struct Mygsl_integration_params{
	double (*f)(double * arg, size_t n_dims, void * params);
	double * params;
};

//=============Xsection base class===================================================
// This is the base class for 2->2 and 2->3 cross-sections.
// It takes care of the tabulating details and the tabulating routines, also the interpolation process
// The actually total Xsection calcuate function and final state sample function are virtual functions, because 2->2 and 2->3 uses quite different techniques to do these jobs.
class Xsection{
protected:
	virtual void tabulate(size_t T_start, size_t dnT) = 0;
	double (*dXdPS)(double * PS, size_t n_dims, void * params);
	double (*approx_X)(double * arg, double M);
	double M1;
public:
	Xsection(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double *, double), double M1_, std::string name_);
	double get_M1(void) {return M1;};
	// arg = [s, T] fot X22, arg = [s, T, dt] for X23, arg = [s, T, s1k, s2k] for f32
	virtual double interpX(double * arg) = 0; 
	virtual double calculate(double * arg) = 0;
	virtual void sample_dXdPS(double * arg, std::vector<std::vector<double> > & fs) = 0;
};

//============Derived 2->2 Xsection class============================================
class Xsection_2to2 : public Xsection{
private:
	rejection_1d sampler1d;
	std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist_phi3;
	std::vector< std::vector<double> > Xtab;
	void tabulate(size_t T_start, size_t dnT);
	const size_t Nsqrts, NT;
	const double sqrtsL, sqrtsM, sqrtsH, dsqrts1, dsqrts2;
	const double TL, TH, dT;
public:
    Xsection_2to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double *, double), double M1_, std::string name_);
	double interpX(double * arg);
    double calculate(double * arg);
	void sample_dXdPS(double * arg, std::vector< std::vector<double> > & final_states);
};

//============Derived 2->3 Xsection class============================================
class Xsection_2to3 : public Xsection{
private:
	AiMS sampler;
	std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist_phi4;
	std::vector< std::vector< std::vector<double> > > Xtab;
	void tabulate(size_t T_start, size_t dnT);
	const size_t Nsqrts, NT, Ndt;
	const double sqrtsL, sqrtsH, dsqrts,
				 TL, TH, dT,
				 dtL, dtH, ddt;
public:
    Xsection_2to3(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double *, double), double M1_, std::string name_);
	double interpX(double * arg);
    double calculate(double * arg);
	void sample_dXdPS(double * arg, std::vector< std::vector<double> > & final_states);
};

//============Derived 3->2 Xsection class============================================
class f_3to2 : public Xsection{
private:
	std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist_phi4;
	const size_t Nsqrts, NT, Nkx, Nkz;
	const double sqrtsL, sqrtsH, dsqrts,
				 TL, TH, dT,
				 kxL, kxH, dkx,
				 kzL, kzH, dkz;
	std::vector<std::vector< std::vector< std::vector<double> > > > Xtab;
	void tabulate(size_t T_start, size_t dnT);
public:
    f_3to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double *, double), double M1_, std::string name_);
	double interpX(double * arg);
    double calculate(double * arg);
	void sample_dXdPS(double * arg, std::vector< std::vector<double> > & final_states);
};


#endif
