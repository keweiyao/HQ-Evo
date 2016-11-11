#include <iostream>
#include <fstream>
#include <cmath>
#include <thread>
#include <vector>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include "Xsection.h"
#include "constants.h"

double gsl_1dfunc_wrapper(double x, void * params_){
	integrate_params * params = static_cast<integrate_params*>(params_);
	return params->dXdPS(&x, 1, params->params);
}

//=============Xsection base class===================================================
// this is the base class for 2->2 and 2->3 cross-sections
Xsection::Xsection(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double, double, double), double M1_, std::string name_)
: dXdPS(dXdPS_), approx_X(approx_X_), M1(M1_), Ns(50), NT(40), 
	sL(M1*M1*1.01), sM(25.*M1*M1), sH(400*M1*M1), 
	ds1((sM-sL)/(Ns-1.)), ds2((sH-sM)/(Ns-1.)),
	TL(0.1), TH(1.0), dT((TH-TL)/(NT-1.))
{
	std::cout << __func__<< " " << name_  << std::endl;
	Xtab.resize(Ns*2);
	for (auto&& ele : Xtab) ele.resize(NT);
}

void Xsection::tabulate_s_Temp(size_t T_start, size_t dnT){
	double s, Temp;
	for (size_t i=0; i<2*Ns; ++i) {
		if (i<Ns) s = sL + i*ds1;
		else s = sM + (i-Ns)*ds2;
		for (size_t j=T_start; j<(T_start+dnT); j++) {
			Temp = TL + j*dT;
			Xtab[i][j] = calculate(s, Temp)/approx_X(s, Temp, M1);
		}
	}
}

double Xsection::interpX(double s, double Temp){
	if (Temp < TL) Temp = TL;
	if (Temp >= TH) Temp = TH-dT;
	if (s < sL) s = sL;
	if (s > sH) s = sH-ds2;
	double xT, rT, xs, rs, ds, smin;
	size_t iT, is, Noffsets;
	xT = (Temp-TL)/dT;	iT = floor(xT); rT = xT - iT;
	if (s < sM) {ds = ds1; smin=sL; Noffsets=0;}
	else {ds = ds2; smin=sM; Noffsets=Ns;}
	xs = (s - smin)/ds; is = floor(xs); rs = xs - is; is += Noffsets;
	
	return approx_X(s, Temp, M1)*(Xtab[is][iT]*(1.-rs)*(1.-rT)
			+Xtab[is+1][iT]*rs*(1.-rT)
			+Xtab[is][iT+1]*(1.-rs)*rT
			+Xtab[is+1][iT+1]*rs*rT);
}

//============Derived 2->2 Xsection class===================================
Xsection_2to2::Xsection_2to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double, double, double), double M1_, std::string name_)
:	Xsection(dXdPS_, approx_X_, M1_, name_)
{
	std::vector<std::thread> threads;
	size_t Ncores = std::thread::hardware_concurrency();
	size_t call_per_core = std::ceil(NT*1./Ncores);
	size_t call_for_last_core = NT - call_per_core*(Ncores-1);
	for (size_t i=0; i< Ncores ; i++)
	{	
		size_t Nstart = i*call_per_core;
		size_t dN = (i==Ncores-1)? call_for_last_core : call_per_core;
		auto code = [this](size_t NTstart_, size_t dNT_) { this->tabulate_s_Temp(NTstart_, dNT_); };
		threads.push_back( std::thread(code, Nstart, dN) );
	}
	for (std::thread& t : threads)	t.join();

	std::ofstream file(name_);
	for (auto roll : Xtab) {
		for (auto item : roll) {
			file << item << " ";
		}
		file << std::endl;
	}
}

double Xsection_2to2::calculate(double s, double Temp){
	double result, error, tmin, tmax;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	integrate_params * params = new integrate_params;
	params->dXdPS = dXdPS;
	double * p = new double[3];
	p[0] = s;
	p[1] = Temp;
	p[2] = M1;
	params->params = p;

    gsl_function F;
	F.function = gsl_1dfunc_wrapper;
	F.params = params;
	tmax = 0.0;
	tmin = -pow(s-M1*M1, 2)/s;
	gsl_integration_qag(&F, tmin, tmax, 0, 1e-4, 1000, 6, w, &result, &error);

	delete[] p;
	delete params;

    return result;
}

double Xsection_2to2::sample_dXdPS(double s, double Temp){
	(void)s;
	(void)Temp;
	return 1.0;
}

//============Derived 2->3 Xsection class===================================
Xsection_2to3::Xsection_2to3(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double, double, double), double M1_, std::string name_)
:	Xsection(dXdPS_, approx_X_, M1_, name_)
{
	std::vector<std::thread> threads;
	size_t Ncores = std::thread::hardware_concurrency();
	size_t call_per_core = std::ceil(NT*1./Ncores);
	size_t call_for_last_core = NT - call_per_core*(Ncores-1);
	for (size_t i=0; i< Ncores ; i++)
	{	
		size_t Nstart = i*call_per_core;
		size_t dN = (i==Ncores-1)? call_for_last_core : call_per_core;
		auto code = [this](size_t NTstart_, size_t dNT_) { this->tabulate_s_Temp(NTstart_, dNT_); };
		threads.push_back( std::thread(code, Nstart, dN) );
	}
	
	for (std::thread& t : threads)	t.join();

	std::ofstream file(name_);
	for (auto roll : Xtab) {
		for (auto item : roll) {
			file << item << " ";
		}
		file << std::endl;
	}
}

double Xsection_2to3::calculate(double s, double Temp){
	double result, error;

	const gsl_rng_type * Tr = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(Tr);
	
	double * params = new double[3];
	params[0] = s; params[1] = Temp; params[2] = M1;
	
	gsl_monte_function G;
	G.f = dXdPS; 
	G.dim = 4; // k, p4, phi4k, cos4
	G.params = params;
	
	// limits of the integration
	double sqrts = std::sqrt(s), M2 = M1*M1;
	double xl[4], xu[4]; // k+p4, k-p4, phi4k, cos4
	xl[0] = 0.5*sqrts*(1.-M2/s); xu[0] = sqrts - M1;
	xl[1] = -0.5*sqrts*(1.-M2/s); xu[1] = 0.5*sqrts*(1.-M2/s);
	xl[2] = 0.0; xu[2] = 2.*M_PI;
	xl[3] = -1.; xu[3] = 1.;
	double Jacobian = 0.5;
	
	// Actuall integration, require the Xi-square to be close to 1,  (0.5, 1.5) 
	gsl_monte_vegas_state * sv = gsl_monte_vegas_alloc(4);
	gsl_monte_vegas_integrate(&G, xl, xu, 4, 10000, r, sv, &result, & error);
	while(std::abs(gsl_monte_vegas_chisq(sv)-1.0)>0.5)
	{
		gsl_monte_vegas_integrate(&G, xl, xu, 4, 10000, r, sv, &result, & error);
	}
	gsl_monte_vegas_free(sv);
	gsl_rng_free(r);
	delete params;
	return result/c256pi4/(s-M2)*Jacobian;
}

double Xsection_2to3::sample_dXdPS(double s, double Temp){
	(void)s;
	(void)Temp;
	return 1.0;
}


