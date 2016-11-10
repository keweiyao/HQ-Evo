#include <iostream>
#include "Xsection.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

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
	tabulate_s_Temp();
}

void Xsection_2to2::tabulate_s_Temp(void){
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

