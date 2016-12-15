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

double interpolate2d(std::vector<std::vector<double> > & A, int ni, int nj, double ri, double rj){
	double wi[2] = {1.-ri, ri}, wj[2] = {1.-rj, rj};
	double result = 0.;
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++){
			result += A[ni+i][nj+j]*wi[i]*wj[j];
		}
	}
	return result;
}

double interpolate3d(std::vector<std::vector<std::vector<double> > > & A, int ni, int nj, int nk, double ri, double rj, double rk){
	double wi[2] = {1.-ri, ri}, wj[2] = {1.-rj, rj}, wk[2] = {1.-rk, rk};
	double result = 0.;
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++){
			for (int k=0; k<2; k++){
				result += A[ni+i][nj+j][nk+k]*wi[i]*wj[j]*wk[k];
			}
		}
	}
	return result;
}

double interpolate4d(std::vector<std::vector<std::vector<std::vector<double> > > > & A, int ni, int nj, int nk, int nt, double ri, double rj, double rk, double rt){
	double wi[2] = {1.-ri, ri}, wj[2] = {1.-rj, rj}, wk[2] = {1.-rk, rk}, wt[3] = {1.-rt, rt};
	double result = 0.;
	for (int i=0; i<2; i++){
		for (int j=0; j<2; j++){
			for (int k=0; k<2; k++){
				for (int t=0; t<2; t++){
					result += A[ni+i][nj+j][nk+k][nt+t]*wi[i]*wj[j]*wk[k]*wt[t];
				}
			}
		}
	}
	return result;
}



double gsl_1dfunc_wrapper(double x, void * params_){
	integrate_params * params = static_cast<integrate_params*>(params_);
	return params->dXdPS(&x, 1, params->params);
}

//=============Xsection base class===================================================
// this is the base class for 2->2 and 2->3 cross-sections
Xsection::Xsection(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double *, double), double M1_, std::string name_)
: dXdPS(dXdPS_), approx_X(approx_X_), M1(M1_), Ns(50), NT(16), 
	sL(M1*M1*1.01), sM(M1*M1*25.), sH(M1*M1*900.), 
	ds1((sM-sL)/(Ns-1.)), ds2((sH-sM)/(Ns-1.)),
	TL(0.12), TH(0.8), dT((TH-TL)/(NT-1.))
{
	std::cout << __func__<< " " << name_  << std::endl;
}


//============Derived 2->2 Xsection class===================================
Xsection_2to2::Xsection_2to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double *, double), double M1_, std::string name_)
:	Xsection(dXdPS_, approx_X_, M1_, name_), rd(), gen(rd()), dist_phi3(0.0, 2.0*M_PI)
{
	Xtab.resize(Ns*2);
	for (auto&& dim1 : Xtab) dim1.resize(NT);

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
	double * arg = new double[2];
	for (size_t i=0; i<2*Ns; i++) {
		if (i<Ns) arg[0] = sL + i*ds1;
		else arg[0] = sM + (i-Ns)*ds2;
		for (size_t j=0; j<NT; j++) {
			arg[1] = TL + j*dT;
			file << Xtab[i][j] * approx_X(arg, M1) << " ";
		}
	}
}

void Xsection_2to2::tabulate_s_Temp(size_t T_start, size_t dnT){
	double * arg = new double[2];
	for (size_t i=0; i<2*Ns; ++i) {
		if (i<Ns) arg[0] = sL + i*ds1;
		else arg[0] = sM + (i-Ns)*ds2;
		for (size_t j=T_start; j<(T_start+dnT); j++) {
			arg[1] = TL + j*dT;
			Xtab[i][j] = calculate(arg)/approx_X(arg, M1);
		}
	}
}

double Xsection_2to2::interpX(double * arg){
	double s = arg[0], Temp = arg[1];
	if (Temp < TL) Temp = TL;
	if (Temp >= TH) Temp = TH-dT;
	if (s < sL) s = sL;
	if (s >= sH) s = sH-ds2;
	double xT, rT, xs, rs, ds, smin;
	size_t iT, is, Noffsets;
	xT = (Temp-TL)/dT;	iT = floor(xT); rT = xT - iT;
	if (s < sM) {ds = ds1; smin=sL; Noffsets=0;}
	else {ds = ds2; smin=sM; Noffsets=Ns;}
	xs = (s - smin)/ds; is = floor(xs); rs = xs - is; is += Noffsets;
	return approx_X(arg, M1)*(Xtab[is][iT]*(1.-rs)*(1.-rT)
			+Xtab[is+1][iT]*rs*(1.-rT)
			+Xtab[is][iT+1]*(1.-rs)*rT
			+Xtab[is+1][iT+1]*rs*rT);
}


double Xsection_2to2::calculate(double * arg){
	double s = arg[0], Temp = arg[1];
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
	gsl_integration_workspace_free(w);

    return result;
}

void Xsection_2to2::sample_dXdPS(double * arg, std::vector< std::vector<double> > & final_states){
	double s = arg[0], Temp = arg[1];
	double * p = new double[3]; //s, T, M
	p[0] = s; p[1] = Temp;  p[2] = M1;
	double psq = std::pow(s-M1*M1, 2)/4./s;
	double pQ = std::sqrt(psq);
	double t = sampler1d.sample(dXdPS, -std::pow(s-M1*M1, 2)/s, 0.0, p);
	double costheta3 = 1. + t/psq/2.;
	double sintheta3 = std::sqrt(1. - costheta3*costheta3);
	double phi3 = dist_phi3(gen);
	double cosphi3 = std::cos(phi3), sinphi3 = std::sin(phi3);
	final_states.resize(1);
	final_states[0].resize(4);
	final_states[0][0] = std::sqrt(psq + M1*M1);
	final_states[0][1] = pQ*sintheta3*cosphi3;
	final_states[0][2] = pQ*sintheta3*sinphi3;
	final_states[0][3] = pQ*costheta3;
	delete[] p;
}

//============Derived 2->3 Xsection class===================================
Xsection_2to3::Xsection_2to3(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double *, double), double M1_, std::string name_)
:	Xsection(dXdPS_, approx_X_, M1_, name_), rd(), gen(rd()), dist_phi4(0.0, 2.0*M_PI), Ndt(10), dtL(0.1), dtH(5.0), ddt((dtH-dtL)/(Ndt-1.))
{
	Xtab.resize(Ns*2);
	for (auto&& dim1 : Xtab) {
		dim1.resize(NT);
		for (auto && dim2 : dim1) 
			dim2.resize(Ndt);
	}

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
	double * arg = new double[3];
	for (size_t i=0; i<2*Ns; i++) {
		if (i<Ns) arg[0] = sL + i*ds1; else arg[0] = sM + (i-Ns)*ds2;
		for (size_t j=0; j<NT; j++) {
			arg[1] = TL + j*dT;
			for (size_t k=0; k<Ndt; k++) {
				arg[2] = dtL + k*ddt;
				file << Xtab[i][j][k] * approx_X(arg, M1) << " ";
			}
		}
	}
}

void Xsection_2to3::tabulate_s_Temp(size_t T_start, size_t dnT){
	double * arg = new double[3]; // s, T, dt
	for (size_t i=0; i<2*Ns; i++) {
		if (i<Ns) arg[0] = sL + i*ds1; else arg[0] = sM + (i-Ns)*ds2;
		for (size_t j=T_start; j<(T_start+dnT); j++) {
			arg[1] = TL + j*dT;
			for (size_t k=0; k<Ndt; k++) {
				arg[2] = dtL + k*ddt;
				Xtab[i][j][k] = calculate(arg)/approx_X(arg, M1);
			}
		}
	}
}

double Xsection_2to3::interpX(double * arg){
	double s = arg[0], Temp = arg[1], dt = arg[2];
	if (Temp < TL) Temp = TL;
	if (Temp >= TH) Temp = TH-dT;
	if (dt < dtL) dt = TL;
	if (dt >= dtH) dt = dtH-ddt;
	if (s < sL) s = sL;
	if (s >= sH) s = sH-ds2;
	double xT, rT, xs, rs, ds, smin, xdt, rdt;
	size_t iT, is, Noffsets, idt;
	xT = (Temp-TL)/dT;	iT = floor(xT); rT = xT - iT;
	xdt = (dt-dtL)/ddt;	idt = floor(xdt); rdt = xdt - idt;
	if (s < sM) {ds = ds1; smin=sL; Noffsets=0;}
	else {ds = ds2; smin=sM; Noffsets=Ns;}
	xs = (s - smin)/ds; is = floor(xs); rs = xs - is; is += Noffsets;
	return approx_X(arg, M1)*( Xtab[is][iT][idt]*(1.-rs)*(1.-rT)*(1.-rdt)	+ Xtab[is][iT][idt+1]*(1.-rs)*(1.-rT)*rdt
								 + Xtab[is+1][iT][idt]*rs*(1.-rT)*(1.-rdt) 		+ Xtab[is+1][iT][idt+1]*rs*(1.-rT)*rdt
								 + Xtab[is][iT+1][idt]*(1.-rs)*rT*(1.-rdt) 		+ Xtab[is][iT+1][idt+1]*(1.-rs)*rT*rdt
								 + Xtab[is+1][iT+1][idt]*rs*rT*(1.-rdt)			+ Xtab[is+1][iT+1][idt+1]*rs*rT*rdt	
								);
}

double Xsection_2to3::calculate(double * arg){
	double s = arg[0], Temp = arg[1], dt = arg[2];
	double result, error;

	const gsl_rng_type * Tr = gsl_rng_default;
	gsl_rng * r = gsl_rng_alloc(Tr);
	
	double * params = new double[4];
	params[0] = s; params[1] = Temp; params[2] = M1; params[3] = dt;
	
	gsl_monte_function G;
	G.f = dXdPS; 
	G.dim = 4; // k, p4, phi4k, cos4
	G.params = params;
	
	// limits of the integration
	double sqrts = std::sqrt(s), M2 = M1*M1;
	double xl[4], xu[4]; // (k+p4), k-p4, phi4k, cos4
	xl[0] = 0.5*sqrts*(1.-M2/s); xu[0] = sqrts-M1;
	xl[1] = -0.5*sqrts*(1.-M2/s); xu[1] = 0.5*sqrts*(1.-M2/s);
	xl[2] = 0.0; xu[2] = 2.*M_PI;
	xl[3] = -1.; xu[3] = 1.;
	
	// Actuall integration, require the Xi-square to be close to 1,  (0.5, 1.5) 
	gsl_monte_vegas_state * sv = gsl_monte_vegas_alloc(4);
	do{ 
		gsl_monte_vegas_integrate(&G, xl, xu, 4, 10000, r, sv, &result, &error);
	}while(std::abs(gsl_monte_vegas_chisq(sv)-1.0)>0.5); 
	gsl_monte_vegas_free(sv);
	gsl_rng_free(r);
	delete params;
	return result/c256pi4/(s-M2);
}

void Xsection_2to3::sample_dXdPS(double * arg, std::vector< std::vector<double> > & final_states){
	// for 2->3, dXdPS is a 5-dimensional distribution,
	// In center of mass frame:
	// there is an overall azimuthal symmetry which allows a flat sampling 
	// the rese 4 variables are sampled from Affine-invariant MCMC procedure,
	// since the distribution scale of each variable could vary a lot.
	// returns heavy quark 4-momentum and radiated gluon 4-momentum
	double s = arg[0], Temp = arg[1], dt = arg[2];
	double * p = new double[4]; //s, T, dt, M
	double sqrts = std::sqrt(s);
	double M2 = M1*M1;
	p[0] = s; p[1] = Temp; p[2] = M1; p[3] = dt; // dt in CoM frame
	size_t n_dims = 4;
	double * guessl = new double[n_dims];
	double * guessh = new double[n_dims];
	double scale1 = 0.5*sqrts*(1.0 - M2/s);
	double scale2 = sqrts-M1;
	guessl[0] = scale1; guessl[1] = -scale1; guessl[2] = M_PI; guessl[3] = -1.0;
	guessh[0] = scale1 + ( scale2 - scale1 )*0.1; guessh[1] = -scale1*0.9; guessh[2] = 2.0*M_PI; guessh[3] = -0.5;
	double * vec4 = sampler.sample(dXdPS, n_dims, p, guessl, guessh);
	double k = 0.5*(vec4[0]+vec4[1]), p4 = 0.5*(vec4[0]-vec4[1]), phi4k = vec4[2], cos4 = vec4[3];
	double cos_star = ((s-M2)-2.*sqrts*(p4+k))/(2.*p4*k) +1.;
	double sin_star = std::sqrt(1. - cos_star*cos_star), sin4 = std::sqrt(1. - cos4*cos4);
	double cos_4k = std::cos(phi4k), sin_4k = std::sin(phi4k);
	// k-vec	
	double kxp = k*(sin_star*cos_4k*cos4 + sin4*cos_star), 
		   kyp = k*sin_star*sin_4k,
		   kz = k*(-sin_star*cos_4k*sin4 + cos4*cos_star);
	// HQ-vec
	double HQxp = -kxp - p4*sin4,
		   HQyp = -kyp,
		   HQz = -kz - p4*cos4,
		   EQ = std::sqrt(HQxp*HQxp+HQyp*HQyp+HQz*HQz+M2);
	// --- randomize the azimuthal angle phi4----
	double phi4 = dist_phi4(gen);
	double cos_phi4 = std::cos(phi4), sin_phi4 = std::sin(phi4);
	double kx = kxp*cos_phi4 + kyp*sin_phi4, ky = -kxp*sin_phi4 + kyp*cos_phi4;
	double HQx = HQxp*cos_phi4 + HQyp*sin_phi4, HQy = -HQxp*sin_phi4 + HQyp*cos_phi4;
	final_states.resize(2);
	final_states[0].resize(4); final_states[1].resize(4);
	final_states[0][0] = EQ; final_states[0][1] = HQx; 
	final_states[0][2] = HQy; final_states[0][3] = HQz;

	final_states[1][0] = k; final_states[1][1] = kx; 
	final_states[1][2] = ky; final_states[1][3] = kz;
	delete[] p;
}


//============Derived 3->2 Xsection class===================================
f_3to2::f_3to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double *, double), double M1_, std::string name_)
:	Xsection(dXdPS_, approx_X_, M1_, name_), rd(), gen(rd()), dist_phi4(0.0, 2.0*M_PI), Ns1k(10), s1kL(M1_*M1_*1.01), s1kH(M1_*M1_*100.), ds1k((s1kH-s1kL)/(Ns1k-1.)), Ns2k(10), s2kL(0.0), s2kH(100.), ds2k((s2kH-s2kL)/(Ns2k-1.))
{
	Xtab.resize(Ns*2);
	for (auto&& dim1 : Xtab) {
		dim1.resize(NT);
		for (auto && dim2 : dim1){
			dim2.resize(Ns1k);
			for (auto && dim3 : dim2){
				dim3.resize(Ns2k);
			}
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
		auto code = [this](size_t NTstart_, size_t dNT_) { this->tabulate_s_Temp(NTstart_, dNT_); };
		threads.push_back( std::thread(code, Nstart, dN) );
	}
	
	for (std::thread& t : threads)	t.join();

	std::ofstream file(name_);
	double * arg = new double[4];
	for (size_t i=0; i<2*Ns; i++) {
		if (i<Ns) arg[0] = sL + i*ds1; else arg[0] = sM + (i-Ns)*ds2;
		for (size_t j=0; j<NT; j++) {
			arg[1] = TL + j*dT;
			for (size_t k=0; k<Ns1k; k++) {
				arg[2] = s1kL + k*ds1k;
				for (size_t t=0; t<Ns2k; t++) {
					arg[3] = s2kL + t*ds2k;
					file << Xtab[i][j][k][t] << " ";
				}
			}
		}
	}
}

void Xf_3to2::tabulate_s_Temp(size_t T_start, size_t dnT){
	double * arg = new double[4];// s, T, s1k, s2k
	for (size_t i=0; i<2*Ns; i++) {
		if (i<Ns) arg[0] = sL + i*ds1; else arg[0] = sM + (i-Ns)*ds2;
		for (size_t j=0; j<NT; j++) {
			arg[1] = TL + j*dT;
			for (size_t k=0; k<Ns1k; k++) {
				arg[2] = s1kL + k*ds1k;
				for (size_t t=0; t<Ns2k; t++) {
					arg[3] = s2kL + t*ds2k;
					Xtab[i][j][k][t] = calculate(arg);
				}
			}
		}
	}
}

double f_3to2::interpX(double * arg){
	double s = arg[0], Temp = arg[1], s1k = arg[2], s2k = arg[3];
	if (Temp < TL) Temp = TL;
	if (Temp >= TH) Temp = TH-dT;
	if (dt < dtL) dt = TL;
	if (dt >= dtH) dt = dtH-ddt;
	if (s < sL) s = sL;
	if (s >= sH) s = sH-ds2;
	double xT, rT, xs, rs, ds, smin, xs1k, rs1k, xs2k, rs2k;
	size_t iT, is, Noffsets, is1k, is2k;
	xT = (Temp-TL)/dT;	iT = floor(xT); rT = xT - iT;
	xs1k = (s1k-s1kL)/ddt;	is1k = floor(xs1k); rs1k = xs1k - is1k;
	xs2k = (s2k-s2kL)/ddt;	is2k = floor(xs2k); rs2k = xs2k - is2k;
	if (s < sM) {ds = ds1; smin=sL; Noffsets=0;}
	else {ds = ds2; smin=sM; Noffsets=Ns;}
	xs = (s - smin)/ds; is = floor(xs); rs = xs - is; is += Noffsets;

	return interpolate4d(Xtab, is, iT, is1k, is2k, rs, rT, rs1k, rs2k);
}

double f_3to2::calculate(double * arg){

}

void f_3to2::sample_dXdPS(double * arg, std::vector< std::vector<double> > & final_states){

}

