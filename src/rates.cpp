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

//=============function wrapper for GSL integration======================
double fy_wrapper(double y, void * params_){
	// unpack Vegas params
	integrate1d_params * params = static_cast<integrate1d_params *>(params_);
	double coeff = params->params[0];
	double Temp = params->params[1];
	double M2 = params->params[2];
	double s = M2 + coeff*y;
	double Xsection = params->interp(s, Temp);
	//delete[] params;
	return y*Xsection;
} 

double fx_wrapper(double x, void * px_){
	integrate1d_params * px = static_cast<integrate1d_params *>(px_);
	double E1 = px->params[0];
	double v1 = px->params[1];
	double Temp = px->params[2];
	double M2 = px->params[3];
	double zeta = px->params[4];

	double result, error, ymin, ymax;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
	integrate1d_params * py = new integrate1d_params;
	py->interp = px->interp;
	py->params = new double[3];
	py->params[0] = 2.*E1*x*Temp;
	py->params[1] = Temp;
	py->params[2] = M2;

    gsl_function F;
	F.function = fy_wrapper;
	F.params = py;
	ymax = 1.+v1;
	ymin = 1.-v1;
	gsl_integration_qag(&F, ymin, ymax, 0, 1e-4, 2000, 6, w, &result, &error);

	delete py;
	//delete px;
	gsl_integration_workspace_free(w);
	return x*x*f_0(x, zeta)*result;
}

//=======================Scattering Rate================================

template <class T>
rates<T>::rates(T * Xprocess_, int degeneracy_, std::string name_)
:	Xprocess(Xprocess_), M(Xprocess->get_M1()), degeneracy(degeneracy_),
	NE1(50), NT(40), E1L(M*1.01), E1H(M*50), TL(0.1), TH(1.0),
	dE1((E1H-E1L)/(NE1-1.)), dT((TH-TL)/(NT-1.)),
	rd(), gen(rd()),
	dist_x(3.0, 1.0), dist_norm_y(-1.0, 1.0), dist_reject(0.0, 1.0)
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
	file.close();
}

template <class T>
void rates<T>::tabulate_E1_T(size_t T_start, size_t dnT){
	for (size_t i=0; i<NE1; i++){
		double E1 = E1L + i*dE1;
		for (size_t j=T_start; j<(T_start+dnT); j++){
			double Temp = TL + j*dT;		
			double result = calculate(E1, Temp);
			Rtab[i][j] = result;
		}
	}
}

template <class T>
double rates<T>::interpR(double E1, double Temp){
	if (Temp < TL) Temp = TL;
	if (Temp >= TH) Temp = TH-dT;
	if (E1 < E1L) E1 = E1L;
	if (E1 >= E1H) E1 = E1H-dE1;
	double xT, rT, xE1, rE1;
	size_t iT, iE1;
	xT = (Temp-TL)/dT;	iT = floor(xT); rT = xT - iT;
	xE1 = (E1 - E1L)/dE1; iE1 = floor(xE1); rE1 = xE1 - iE1;
	return   Rtab[iE1][iT]*(1.-rE1)*(1.-rT)
			+Rtab[iE1+1][iT]*rE1*(1.-rT)
			+Rtab[iE1][iT+1]*(1.-rE1)*rT
			+Rtab[iE1+1][iT+1]*rE1*rT;
}

template <class T>
double rates<T>::calculate(double E1, double Temp)
{
	double p1 = std::sqrt(E1*E1-M*M);
	double result, error, xmin, xmax;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
	integrate1d_params * px = new integrate1d_params;
	px->interp = std::bind( &Xsection_2to2::interpX, Xprocess, _1, _2);
	px->params = new double[5];
	px->params[0] = E1;
	px->params[1] = p1/E1;
	px->params[2] = Temp;
	px->params[3] = M*M;
	px->params[4] = 0.0;

    gsl_function F;
	F.function = fx_wrapper;
	F.params = px;
	xmax = 6.0;
	xmin = 0.0;
	gsl_integration_qag(&F, xmin, xmax, 0, 1e-4, 2000, 6, w, &result, &error);

	gsl_integration_workspace_free(w);
	delete px;
	return result*std::pow(Temp, 3)*4./c16pi2*E1/p1*degeneracy;
}

template <class T>
void rates<T>::sample_initial(double E1, double Temp, double &E2, double &s){
	// this function samples x = E2/T and y = (s-M^2)/(2*E2*E1) from the distribution:
	// P(x, y) ~ x^2*exp(-x) * y*sigma(M^2 + 2*E1*T*x*y, T)
	// We first generate X from gamma distribution Gamma(x; 3,1) ~ x^3*exp(-x) (cut off x < 20. )
	// and uniform sample y within (1-v1, 1+v1)
	// and finally rejected with P_rej(x,y) = y*sigma(M^2 + 2*E1*T*x*y, T);
	double M2 = M*M, x, y, max, smax, stemp, coeff = 2.*E1*Temp;
	double v1 = std::sqrt(E1*E1 - M2)/E1;
	smax = M2 + coeff*20.*(1.+v1);
	if (smax < 2.*M2) smax = 2.*M2;
	max = (1.+v1)*Xprocess->interpX(smax, Temp);
	double Nd=0.;
	do{
		do{ x = dist_x(gen); }while(x>20.);
		y = 1. + dist_norm_y(gen)*v1;
		stemp = M2 + coeff*x*y;
		Nd += 1.;
	}while( y*Xprocess->interpX(stemp, Temp) <= max*dist_reject(gen) );
	E2 = x*Temp;
	s = M2 + coeff*x*y;
}

template class rates<Xsection_2to2>;
template class rates<Xsection_2to3>;


