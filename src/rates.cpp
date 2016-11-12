#include <cmath>
#include <vector>
#include <thread> 
#include <fstream>
#include <string>
#include "rates.h"

using std::placeholders::_1;
using std::placeholders::_2;

//=============Thernalized Distribution funtion=================================
// xi = 1: Fermi Dirac; xi = -1 Bose Einsterin; xi = 0, Maxwell-Boltzmann
double inline f_0(double x, double xi){
    if (x<1e-9) x=1e-9;
    return 1./(std::exp(x)+xi);
}

//=============function wrapper for GSL vegas integration======================
double Vegas_func_wrapper(double * var, size_t n_dims, void * params_)
{
	(void) n_dims;
	// unpack variabel x = E2/T, y = (s-M2)/2/E2/E1
	double x = var[0];
	double y = var[1];
	
	// unpack Vegas params
	Vegas_params * params = static_cast<Vegas_params *>(params_);
	double * p = static_cast<double*>(params->params);
	double E1 = p[0];
	double Temp = p[1];
	double M2 = p[2]*p[2];
	double E2 = x*Temp;
	double s = M2 + 2*E1*E2*y;
	
	if (x < 0 || x > 5.) std::cout << "?" << std::endl;
	double Xsection = params->interp(s, Temp);
	
	return x*x*f_0(x, p[3])*y*Xsection;
} 
//=======================Scattering Rate================================

template <class T>
rates<T>::rates(T * Xprocess_, int degeneracy_, std::string name_)
:	Xprocess(Xprocess_),
	degeneracy(degeneracy_),
	NE1(50),
	NT(40)
{
	//Parallel tabulating scattering rate (each core is resonpible for several temperatures)
	// for the first n-1 cores, each takes care of m Temps.
	// the last core could take less jobs
	std::cout << __func__ << " " << name_ << std::endl;
	Rtab.resize(NE1);
	for (auto&& R : Rtab){
		R.resize(NT);
	}
	std::vector<std::thread> threads;
	size_t Ncores = std::thread::hardware_concurrency();
	size_t call_per_core = std::ceil(NT*1./Ncores);
	size_t call_for_last_core = NT - call_per_core*(Ncores-1);
	for (size_t i=0; i< Ncores ; i++)
	{	
		size_t Nstart = i*call_per_core;
		size_t dN = (i==Ncores-1)? call_for_last_core : call_per_core;
		auto code = [this](size_t NTstart_, size_t dNT_) { this->tabulate_E1_T(NTstart_, dNT_); };
		threads.push_back( std::thread(code, Nstart, dN) );
	}
	
	for (std::thread& t : threads)	t.join();
	
	std::ofstream file(name_);
	for (auto roll : Rtab) {
		for (auto item : roll) {
			file << item << " ";
		}
		file << std::endl;
	}
}

template <class T>
void rates<T>::tabulate_E1_T(size_t T_start, size_t dnT){
	for (size_t i=0; i<NE1; i++){
		double E1 = 1.301+1*i;
		for (size_t j=T_start; j<(T_start+dnT); j++){
			double Temp = 0.1+0.02*j;		
			double result = calculate(E1, Temp);
			Rtab[i][j] = result;
		}
	}
}

template <class T>
double rates<T>::calculate(double E1, double Temp)
{
	double result, error;
		
	double M = 1.3; 
	double p1 = std::sqrt(E1*E1-M*M);
	const gsl_rng_type * Tr = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(Tr);
	
	// wrap parameters and interpolating function to Vega_params and Vega_func_wrapper
	Vegas_params * params = new Vegas_params;
	params->interp = std::bind( &Xsection_2to2::interpX, Xprocess, _1, _2);
	double *p = new double[4];
	p[0] = E1; p[1] = Temp; p[2] = M; p[3] = 0.0;
	params->params = p;
	
	gsl_monte_function G;
	G.f = &Vegas_func_wrapper;
	G.dim = 2;
	G.params = params;
	
	// limits of the integration
	// variables:x = E2/T,  y = (s-M2)/2/E2/E1
	double xl[2], xu[2];
	xl[0] = 0.0; xu[0] = 5.;
	xl[1] = 1. - p1/E1; xu[1] = 1 + p1/E1; 
	
	// Actuall integration, require the Xi-square to be close to 1,  (0.5, 1.5) 
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);
	gsl_monte_vegas_integrate(&G, xl, xu, 2, 10000, r, s, &result, & error);
	while(std::abs(gsl_monte_vegas_chisq(s)-1.0)>0.5)
	{
		gsl_monte_vegas_integrate(&G, xl, xu, 2, 10000, r, s, &result, & error);
	}
	gsl_monte_vegas_free(s);
	gsl_rng_free(r);
	delete params;
	delete[] p;
	return result*std::pow(Temp, 3)*4./c16pi2*E1/p1*degeneracy;
}

template <class T>
void rates<T>::sample_initial(double E1, double Temp, double &E2, double &s){
}

template class rates<Xsection_2to2>;
template class rates<Xsection_2to3>;


