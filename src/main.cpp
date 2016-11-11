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
	p[0] = s; p[1] = 0.3;  p[2] = M;
	size_t n_dims = 4;
	double * guess = new double[n_dims];
	guess[0] = (sqrts/2. - M/2.+M*M/sqrts/2.)*2;
	guess[1] = (sqrts/2. - M/2.+M*M/sqrts/2.)*0.1;
	guess[2] = M_PI;
	guess[3] = 0.1;
	AiMS sampler;
	for (size_t i = 0; i<1000; i++)
	sampler.sample(M2_Qq2Qqg, n_dims, p, guess);	
	/*
	// Xsection: Q+q->Q+q
	Xsection_2to2 XQq2Qq(&dX_Qq2Qq_dPS, &approx_XQq2Qq, 1.3, "./tables/X-Qq-Qq.dat");
	// Xsection: Q+g->Q+g
	Xsection_2to2 XQg2Qg(&dX_Qg2Qg_dPS, &approx_XQg2Qg, 1.3, "./tables/X-Qg-Qg.dat");
	// Xsection: Q+q->Q+q+g
	Xsection_2to3 XQq2Qqg(&M2_Qq2Qqg, &approx_XQq2Qqg, 1.3, "./tables/X-Qq-Qqg.dat");
	// Xsection: Q+g->Q+g+g
	Xsection_2to3 XQg2Qgg(&M2_Qg2Qgg, &approx_XQg2Qgg, 1.3, "./tables/X-Qg-Qgg.dat");

	// Rate: Q+q->Q+q
	rates<Xsection_2to2> RQq2Qq(&XQq2Qq, 3*4, "./tables/R-Qq-Qq.dat");
	// Rate: Q+g->Q+g
	rates<Xsection_2to2> RQg2Qg(&XQg2Qg, 8*2, "./tables/R-Qg-Qg.dat");
	// Rate: Q+q->Q+q+g
	rates<Xsection_2to3> RQq2Qqg(&XQq2Qqg, 3*4, "./tables/R-Qq-Qqg.dat");
	// Rate: Q+g->Q+g+g
	rates<Xsection_2to3> RQg2Qgg(&XQg2Qgg, 8*2, "./tables/R-Qg-Qgg.dat");
*/
	return 0;
}
