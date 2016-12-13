#include <cmath>
#include <vector>
#include <thread> 
#include <fstream>
#include <string>
#include "rates.h"

using std::placeholders::_1;

//=============Thernalized Distribution funtion=================================
// xi = 1: Fermi Dirac; xi = -1 Bose Einsterin; xi = 0, Maxwell-Boltzmann
double inline f_0(double x, double xi){
    if (x<1e-9) x=1e-9;
    return 1./(std::exp(x)+xi);
}

//=============function wrapper for GSL integration R22======================
double fy_wrapper22(double y, void * params_){
	// unpack Vegas params
	integrate1d_params * params = static_cast<integrate1d_params *>(params_);
	double coeff = params->params[0];
	double Temp = params->params[1];
	double M2 = params->params[2];
	double v1 = params->params[3];
	double s = M2 + coeff*(1.-v1*y);
	double * arg = new double[2];
	arg[0] = s; arg[1] = Temp;
	double Xsection = params->interp(arg);
	//delete[] params;
	return (1.-v1*y)*Xsection;
} 

double fx_wrapper22(double x, void * px_){
	integrate1d_params * px = static_cast<integrate1d_params *>(px_);
	double E1 = px->params[0];
	double v1 = px->params[1];
	double Temp = px->params[2];
	double M2 = px->params[3];
	double zeta = px->params[4];

	double result, error, ymin, ymax;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
	integrate1d_params * py = new integrate1d_params;
	py->interp = px->interp;
	py->params = new double[4];
	py->params[0] = 2.*E1*x*Temp;
	py->params[1] = Temp;
	py->params[2] = M2;
	py->params[3] = v1;

    gsl_function F;
	F.function = fy_wrapper22;
	F.params = py;
	ymax = 1.;
	ymin = -1.;
	gsl_integration_qag(&F, ymin, ymax, 0, 1e-3, 10000, 6, w, &result, &error);

	delete py;
	//delete px;
	gsl_integration_workspace_free(w);
	return x*x*f_0(x, zeta)*result;
}

//=============function wrapper for GSL integration R23======================
//=============function wrapper for GSL integration R23======================
double fy_wrapper23(double y, void * params_){
	// unpack Vegas params
	integrate1d_params * params = static_cast<integrate1d_params *>(params_);
	double coeff = params->params[0];
	double Temp = params->params[1];
	double M2 = params->params[2];
	double v1 = params->params[3];
	double dt = params->params[4];// dt in the Cell Frame
	double E2 = params->params[5];// E2 in the Cell Frame
	double E1 = params->params[6];// E1 in the Cell Frame
	double s = M2 + coeff*(1.-v1*y);// s variable
	// transform time separation in to CoM frame
	double costheta2 = (M2 + 2.*E1*E2 - s)/(2.*v1*E1*E2);
	double vz_com = (v1*E1 + E2*costheta2)/(E1+E2);
	double gamma_com = (E1+E2)/std::sqrt(s);
	double dtp = gamma_com*(1. - vz_com*v1)*dt; // dt in the CoM Frame
	double * arg = new double[3];
	arg[0] = s; arg[1] = Temp; arg[2] = dtp; // dt in the CoM Frame
	double Xsection = params->interp(arg);
	//delete[] params;
	return (1.-v1*y)*Xsection;
} 

double fx_wrapper23(double x, void * px_){
	integrate1d_params * px = static_cast<integrate1d_params *>(px_);
	double E1 = px->params[0];
	double v1 = px->params[1];
	double Temp = px->params[2];
	double M2 = px->params[3];
	double zeta = px->params[4];
	double dt = px->params[5]; // dt in the Cell Frame

	double result, error, ymin, ymax;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(10000);
	integrate1d_params * py = new integrate1d_params;
	py->interp = px->interp;
	py->params = new double[7];
	py->params[0] = 2.*E1*x*Temp;
	py->params[1] = Temp;
	py->params[2] = M2;
	py->params[3] = v1;
	py->params[4] = dt; // dt in the Cell Frame
	py->params[5] = x*Temp; // E2 in the Cell Frame
	py->params[6] = E1; // E1 in the Cell Frame

    gsl_function F;
	F.function = fy_wrapper23;
	F.params = py;
	ymax = 1.;
	ymin = -1.;
	gsl_integration_qag(&F, ymin, ymax, 0, 1e-3, 10000, 6, w, &result, &error);

	delete py;
	//delete px;
	gsl_integration_workspace_free(w);
	return x*x*f_0(x, zeta)*result;
}

//=======================Rates abstract class==================================
rates::rates(std::string name_)
:	rd(), gen(rd()),
	dist_x(3.0, 1.0), dist_norm_y(-1.0, 1.0), dist_reject(0.0, 1.0)
{
	std::cout << __func__ << " " << name_ << std::endl;
}

//=======================Derived Scattering Rate class 2 to 2================================
rates_2to2::rates_2to2(Xsection_2to2 * Xprocess_, int degeneracy_, std::string name_)
:	rates(name_), Xprocess(Xprocess_), M(Xprocess->get_M1()), degeneracy(degeneracy_),
	NE1(100), NT(16), E1L(M*1.01), E1H(M*100), TL(0.12), TH(0.8),
	dE1((E1H-E1L)/(NE1-1.)), dT((TH-TL)/(NT-1.))
{
	//Parallel tabulating scattering rate (each core is resonpible for several temperatures)
	// for the first n-1 cores, each takes care of m Temps.
	// the last core could take less jobs
	Rtab.resize(NE1);
	for (auto&& dim1 : Rtab){
		dim1.resize(NT);
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
	for (auto dim1 : Rtab)
		for (auto dim2 : dim1)
			file << dim2 << " ";
	file.close();
}

void rates_2to2::tabulate_E1_T(size_t T_start, size_t dnT){
	double * arg = new double[2];
	for (size_t i=0; i<NE1; i++){
		arg[0] = E1L + i*dE1;
		for (size_t j=T_start; j<(T_start+dnT); j++){
			arg[1] = TL + j*dT;		
			Rtab[i][j] = calculate(arg);
		}
	}
}

double rates_2to2::interpR(double * arg){
	double E1 = arg[0], Temp = arg[1];
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

double rates_2to2::calculate(double * arg)
{
	double E1 = arg[0], Temp = arg[1];
	double p1 = std::sqrt(E1*E1-M*M);
	double result, error, xmin, xmax;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
	integrate1d_params * px = new integrate1d_params;
	px->interp = std::bind( &Xsection_2to2::interpX, Xprocess, _1);
	px->params = new double[5];
	px->params[0] = E1;
	px->params[1] = p1/E1;
	px->params[2] = Temp;
	px->params[3] = M*M;
	px->params[4] = 0.0;

    gsl_function F;
	F.function = fx_wrapper22;
	F.params = px;
	xmax = 5.0;
	xmin = 0.0;
	gsl_integration_qag(&F, xmin, xmax, 0, 1e-2, 2000, 6, w, &result, &error);

	gsl_integration_workspace_free(w);
	delete px;
	return result*std::pow(Temp, 3)*4./c16pi2*degeneracy;
}

void rates_2to2::sample_initial(double * arg_in, double * arg_out){
	// this function samples x = E2/T and y = cos(theta2) from the distribution:
	// P(x, y) ~ x^2*exp(-x) * (1-v1*y) * sigma(M^2 + 2*E1*T*x - 2*p1*T*x*y, T)
	// We first generate X from gamma distribution Gamma(x; 3,1) ~ x^3*exp(-x) (cut off x < 20. )
	// and uniform sample y within (-1., 1.)
	// and finally rejected with P_rej(x,y) = (1-v1*y) * sigma(M^2 + 2*E1*T*x - 2*p1*T*x*y, T);
	double E1 = arg_in[0], Temp = arg_in[1];
	double * Xarg = new double[2]; Xarg[1] = Temp;
	double M2 = M*M, x, y, max, smax, stemp;
	double v1 = std::sqrt(E1*E1 - M2)/E1;
	double intersection = M*M, coeff1 = 2.*E1*Temp, coeff2 = -2.*E1*Temp*v1;
	smax = M2 + coeff1*20. + coeff2*20.;
	if (smax < 2.*M2) smax = 2.*M2;
	Xarg[0] = smax;
	max = (1.+v1)*Xprocess->interpX(Xarg);
	do{
		do{ x = dist_x(gen); }while(x>20.);
		y = dist_norm_y(gen);
		stemp = intersection + coeff1*x + coeff2*x*y;
		Xarg[0] = stemp;
	}while( (1.-v1*y)*Xprocess->interpX(Xarg) <= max*dist_reject(gen) );
	arg_out[0] = x*Temp;
	arg_out[1] = intersection + coeff1*x + coeff2*x*y;
}




//=======================Derived Scattering Rate class 2 to 3================================
rates_2to3::rates_2to3(Xsection_2to3 * Xprocess_, int degeneracy_, std::string name_)
:	rates(name_), Xprocess(Xprocess_), M(Xprocess->get_M1()), degeneracy(degeneracy_),
	NE1(100), NT(16), Ndt(10), E1L(M*1.01), E1H(M*100), TL(0.12), TH(0.8), dtL(0.1), dtH(10.0),
	dE1((E1H-E1L)/(NE1-1.)), dT((TH-TL)/(NT-1.)), ddt((dtH-dtL)/(Ndt-1.))
{
	//Parallel tabulating scattering rate (each core is resonpible for several temperatures)
	// for the first n-1 cores, each takes care of m Temps.
	// the last core could take less jobs
	Rtab.resize(NE1);
	for (auto&& dim1 : Rtab){
		dim1.resize(NT);
		for (auto&& dim2 : dim1){
			dim2.resize(Ndt);
		}
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
	for (auto dim1 : Rtab)
		for (auto dim2 : dim1)
			for (auto dim3 : dim2) 
				file << dim3 << " ";
	file.close();
}

void rates_2to3::tabulate_E1_T(size_t T_start, size_t dnT){
	double * arg = new double[2];
	for (size_t i=0; i<NE1; i++){
		arg[0] = E1L + i*dE1;
		for (size_t j=T_start; j<(T_start+dnT); j++){
			arg[1] = TL + j*dT;
			for (size_t k=0; k<Ndt; k++){
				arg[2] = dtL + k*ddt;
				Rtab[i][j][k] = calculate(arg);
			}
		}
	}
}

double rates_2to3::interpR(double * arg){
	double E1 = arg[0], Temp = arg[1], dt = arg[2];
	if (Temp < TL) Temp = TL;
	if (Temp >= TH) Temp = TH-dT;
	if (dt < dtL) dt = dtL;
	if (dt >= dtH) dt = dtH-ddt;
	if (E1 < E1L) E1 = E1L;
	if (E1 >= E1H) E1 = E1H-dE1;
	double xT, rT, xE1, rE1, xdt, rdt;
	size_t iT, iE1, idt;
	xT = (Temp-TL)/dT;	iT = floor(xT); rT = xT - iT;
	xdt = (dt-dtL)/ddt;	idt = floor(xdt); rdt = xdt - idt;
	xE1 = (E1 - E1L)/dE1; iE1 = floor(xE1); rE1 = xE1 - iE1;
	return   Rtab[iE1][iT][idt]*(1.-rE1)*(1.-rT)*(1.-rdt) + Rtab[iE1][iT][idt+1]*(1.-rE1)*(1.-rT)*rdt
			+Rtab[iE1+1][iT][idt]*rE1*(1.-rT)*(1.-rdt) + Rtab[iE1+1][iT][idt+1]*rE1*(1.-rT)*rdt
			+Rtab[iE1][iT+1][idt]*(1.-rE1)*rT*(1.-rdt) + Rtab[iE1][iT+1][idt+1]*(1.-rE1)*rT*rdt
			+Rtab[iE1+1][iT+1][idt]*rE1*rT*(1.-rdt) + Rtab[iE1+1][iT+1][idt+1]*rE1*rT*rdt;
}

double rates_2to3::calculate(double * arg)
{
	double E1 = arg[0], Temp = arg[1], dt = arg[2]; // dt in the Cell Frame
	double p1 = std::sqrt(E1*E1-M*M);
	double result, error, xmin, xmax;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
	integrate1d_params * px = new integrate1d_params;
	px->interp = std::bind( &Xsection_2to3::interpX, Xprocess, _1);
	px->params = new double[6];
	px->params[0] = E1;
	px->params[1] = p1/E1;
	px->params[2] = Temp;
	px->params[3] = M*M;
	px->params[4] = 0.0;
	px->params[5] = dt; // dt in the Cell Frame

    gsl_function F;
	F.function = fx_wrapper23;
	F.params = px;
	xmax = 5.0;
	xmin = 0.0;
	gsl_integration_qag(&F, xmin, xmax, 0, 1e-2, 2000, 6, w, &result, &error);

	gsl_integration_workspace_free(w);
	delete px;
	return result*std::pow(Temp, 3)*4./c16pi2*degeneracy;
}

void rates_2to3::sample_initial(double * arg_in, double * arg_out){
	// this function samples x = E2/T and y = cos(theta2) from the distribution:
	// P(x, y) ~ x^2*exp(-x) * (1-v1*y) * sigma(M^2 + 2*E1*T*x - 2*p1*T*x*y, T)
	// We first generate X from gamma distribution Gamma(x; 3,1) ~ x^3*exp(-x) (cut off x < 20. )
	// and uniform sample y within (-1., 1.)
	// and finally rejected with P_rej(x,y) = (1-v1*y) * sigma(M^2 + 2*E1*T*x - 2*p1*T*x*y, T);
	double E1 = arg_in[0], Temp = arg_in[1], dt = arg_in[2];
	double * Xarg = new double[3]; Xarg[1] = Temp; Xarg[2] = dt; // dt in Cell Frame
	double M2 = M*M, x, y, max, smax, stemp;
	double v1 = std::sqrt(E1*E1 - M2)/E1;
	double intersection = M*M, coeff1 = 2.*E1*Temp, coeff2 = -2.*E1*Temp*v1;
	smax = M2 + coeff1*20. + coeff2*20.;
	if (smax < 2.*M2) smax = 2.*M2;
	Xarg[0] = smax;
	max = (1.+v1)*Xprocess->interpX(Xarg);
	do{
		do{ x = dist_x(gen); }while(x>20.);
		y = dist_norm_y(gen);
		stemp = intersection + coeff1*x + coeff2*x*y;
		Xarg[0] = stemp; 
	}while( (1.-v1*y)*Xprocess->interpX(Xarg) <= max*dist_reject(gen) );
	arg_out[0] = x*Temp;
	arg_out[1] = intersection + coeff1*x + coeff2*x*y;
}

