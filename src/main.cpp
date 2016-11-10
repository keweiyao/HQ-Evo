#include <iostream>
#include <fstream>
#include "sample_methods.h"
#include "matrix_elements.h"
#include "rates.h"

int main(){
	//rejection_1d sampler1d;
	//double * p = new double[2]; p[0] = 0.1; p[1] = 0.1;
	//for(int i=0; i< 1000000; i++) std::cout << sampler1d.sample(test, 0., 1., p) << std::endl;

	//--------test MC sampler------------
	double * p = new double[3]; //s, T, M
	double M = 1.3;
	double sqrts = M*5.;
	double s = sqrts*sqrts;
	p[0] = s; p[1] = 0.2;  p[2] = M;
	size_t n_dims = 4;
	double * guess = new double[n_dims];
	guess[0] = sqrts/2. - M/2.;
	guess[1] = guess[0];
	guess[2] = M_PI;
	guess[3] = 0.2;
	AiMS sampler;
	for (size_t i = 0; i<1000; i++)
	sampler.sample(M2_Qq2Qqg, n_dims, p, guess);	

	//Xsection_2to2 XQg2Qg(&dX_Qg2Qg_dt, &approx_XQg2Qg, 1.3);
	//Xsection_2to2 XQq2Qq(&dX_Qq2Qq_dt, &approx_XQq2Qq, 1.3);
	//rates<Xsection_2to2> RQg2Qg(&XQg2Qg, 8*2);
	//rates<Xsection_2to2> RQq2Qq(&XQq2Qq, 3*4);

	return 0;
}
