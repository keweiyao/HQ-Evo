#ifndef XSECTION_H
#define XSECTION_H

#include <cstdlib>
#include <vector>
#include "matrix_elements.h"

//=============================================================================
class Xsection_2to2{
private:
	void tabulate_s_Temp(void);
	double (*dXdt)(double t, void * params);
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
    Xsection_2to2(double (*dXdt_)(double t, void * params), double (*approx_X_)(double s, double Temp, double M), double M1);
    double calculate(double s, double Temp);
	double interpX(double s, double Temp);
	double sample_dXdt(double s, double Temp);
};

#endif
