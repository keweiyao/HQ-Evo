#include "Rates.h"
#include <cmath>
#include <vector>
#include <thread> 
#include <fstream>
using std::placeholders::_1;
using std::placeholders::_2;
//=============Thernalized Distribution funtion=================================
// xi = 1: Fermi Dirac; xi = -1 Bose Einsterin; xi = 0, Maxwell-Boltzmann
double inline f_0(double x, double xi){
    if (x<1e-6) x=1e-64;
    return 1./(std::exp(x)+xi);
}

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
double dX_Qq2Qq_dt(double t, void * params){
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
	return result/c16pi/std::pow(s-M2, 2);
}
double approx_XQq2Qq(double s, double t, double M){
	(void)s;
	(void)t;
	(void)M;
	return 1.0;
}

//=============Baisc function for Q+g --> Q+g==================================
double dX_Qg2Qg_dt(double t, void * params){
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
	return result*c16pi2/c16pi/std::pow(s-M2, 2);
}
double approx_XQg2Qg(double s, double Temp, double M){	
	double Q2s = s-M*M, T2 = Temp*Temp, M2 = M*M;
	double abstmax = Q2s*Q2s/s;
	return 1./T2 - 1./(T2+abstmax) + 10.*M2/Q2s/Q2s;
}
//=============Xsection_2to2===================================================
// Integrate 2->2 type cross-section and tabulate at grid points, and provide
// fast interpolator for 2->2 type cross-section

// initialize the class with the specifiv 2->2 process you want
// by the t-differential cross-section, Mass of the massive particle
Xsection_2to2::Xsection_2to2(double (*dXdt_)(double t, void * params), double (*approx_X_)(double s, double Temp, double M), double M1_)
:	dXdt(dXdt_), approx_X(approx_X_), M1(M1_), Ns(100), NT(50), 
	sL(M1*M1*1.001), sM(4.*M1*M1), sH(400*M1*M1), 
	ds1((sM-sL)/(Ns-1.)), ds2((sH-sM)/(Ns-1.)),
	TL(0.1), TH(1.0), dT((TH-TL)/(NT-1.))
{
	// making table
	tabulate_st();
}

void Xsection_2to2::tabulate_st(void){
	std::cout << "Making X-section table" << std::endl;
	double s, Temp;
	Xtab.resize(Ns*2);	
	for (size_t i=0; i<2*Ns; ++i) {
		if (i<Ns) s = sL + i*ds1;
		else s = sM + (i-Ns)*ds2;
		Xtab[i].resize(NT);
		for (size_t j=0; j<NT; ++j) {
			Temp = TL + j*dT;
			Xtab[i][j] = calculate(s, Temp)/approx_X(s, Temp, M1);
		}
	}
	std::cout << "X-section done" << std::endl;
}

double Xsection_2to2::calculate(double s, double Temp){
	double result, error, tmin, tmax;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	double * p = new double[3];
	p[0] = s;
	p[1] = Temp;
	p[2] = M1;
    gsl_function F;
	F.function = dXdt;
	F.params = p;
	tmax = 0.0;
	tmin = -pow(s-M1*M1, 2)/s;
	gsl_integration_qag(&F, tmin, tmax, 0, 1e-4, 1000, 6, w, &result, &error);
    return result;
}

double Xsection_2to2::interpX(double s, double Temp){
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

//=============X-section Class function wrapper for GSL vegas integration======================
double Vegas_func_wrapper(double * var, size_t n_dims, void * params_)
{
	(void) n_dims;
	// unpack variabel s and E2
	double s = var[0];
	double E2 = var[1];
	
	// unpack Vegas params
	Vegas_params * params = static_cast<Vegas_params *>(params_);
	double * p = static_cast<double*>(params->params);
	double kp = p[0];
	double km = p[1];
	double Temp = p[2];
	double M1_2 = p[3];

	double Xsection = params->interp(s, Temp);
	double s_max = M1_2 + kp*E2;
	double s_min = M1_2 + km*E2;
	
	if (s <= s_min || s >= s_max) return 0.0;
	else
	    return f_0(E2/Temp, p[4])*(s-M1_2)*Xsection;
} 
//=======================Scattering Rate================================

template <class T>
rates<T>::rates(T * Xprocess_, int degeneracy_)
:	Xprocess(Xprocess_),
	degeneracy(degeneracy_),
	NE1(50),
	NT(40)
{
	//Parallel tabulating scattering rate (each core is resonpible for several temperatures)
	// for the first n-1 cores, each takes care of m Temps.
	// the last core could take less jobs
	
	Rtab.resize(NE1);
	for (auto&& R : Rtab){
		R.resize(NT);
	}
	std::vector<std::thread> threads;
	size_t Ncores = std::thread::hardware_concurrency();
	std::cout << "Available cores" << Ncores << std::endl;
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
	
	std::ofstream file("data_rate.dat");
	for (auto roll : Rtab) {
		for (auto item : roll) {
			file << item << " ";
		}
		file << std::endl;
	}
}

template <class T>
void rates<T>::tabulate_E1_T(size_t T_start, size_t dnT){
	std:: cout << "ID " << std::this_thread::get_id()  << std::endl;
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
		
	double M = 1.3; double p1 = std::sqrt(E1*E1-M*M);
	const gsl_rng_type * Tr = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(Tr);
	
	// wrap parameters and interpolating function to Vega_params and Vega_func_wrapper
	Vegas_params * params = new Vegas_params;
	params->interp = std::bind( &Xsection_2to2::interpX, Xprocess, _1, _2);
	double *p = new double[5];
	p[0] = 2*(E1+p1); p[1] = 2*(E1-p1); p[2] = Temp; p[3] = M*M; p[4] = 0.;
	params->params = p;
	
	gsl_monte_function G;
	G.f = &Vegas_func_wrapper;
	G.dim = 2;
	G.params = params;
	
	// limits of the integration
	double xl[2], xu[2];
	xl[0] = M*M; xu[0] = M*M + 2.*(E1+p1)*5.*Temp; 
	xl[1] = 0.0; xu[1] = 5.*Temp;
	
	// Actuall integration, require the Xi-square to be close to 1,  (0.5, 1.5) 
	gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);
	gsl_monte_vegas_integrate(&G, xl, xu, 2, 20000, r, s, &result, & error);
	while(std::abs(gsl_monte_vegas_chisq(s)-1.0)>0.5)
	{
		gsl_monte_vegas_integrate(&G, xl, xu, 2, 40000, r, s, &result, & error);
	}
	gsl_monte_vegas_free(s);
	gsl_rng_free(r);
	return result/c16pi2/E1/p1*degeneracy;
}

template class rates<Xsection_2to2>;


