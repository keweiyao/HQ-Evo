#include <iostream>
#include <fstream>
#include "sample_methods.h"
#include "matrix_elements.h"
#include "rates.h"

int main(){
	double M = 1.3;
	// Xsection: Q+q->Q+q
	//Xsection_2to2 XQq2Qq(&dX_Qq2Qq_dPS, &approx_XQq2Qq, M, "./tables/X-Qq-Qq.dat");
	// Xsection: Q+g->Q+g
	//Xsection_2to2 XQg2Qg(&dX_Qg2Qg_dPS, &approx_XQg2Qg, M, "./tables/X-Qg-Qg.dat");
	// Xsection: Q+q->Q+q+g
	Xsection_2to3 XQq2Qqg(&M2_Qq2Qqg, &approx_XQq2Qqg, M, "./tables/X-Qq-Qqg.dat");
	double * result = new double[8];
	double s = M*M*400., T = 0.3;
	std::ofstream f("HQ-k.dat")	;
	for (int i=0; i<40000; i++){	
		XQq2Qqg.sample_dXdPS(s, T, result);
		for (int j=0; j<8; j++){
			f << result[j] << " ";
		}
		f << std::endl;
	}
	delete result;
	// Xsection: Q+g->Q+g+g
	//Xsection_2to3 XQg2Qgg(&M2_Qg2Qgg, &approx_XQg2Qgg, M, "./tables/X-Qg-Qgg.dat");

	// Rate: Q+q->Q+q
	//rates<Xsection_2to2> RQq2Qq(&XQq2Qq, 3*4, "./tables/R-Qq-Qq.dat");
	// Rate: Q+g->Q+g
	//rates<Xsection_2to2> RQg2Qg(&XQg2Qg, 8*2, "./tables/R-Qg-Qg.dat");
	// Rate: Q+q->Q+q+g
	//rates<Xsection_2to3> RQq2Qqg(&XQq2Qqg, 3*4, "./tables/R-Qq-Qqg.dat");
	// Rate: Q+g->Q+g+g
	//rates<Xsection_2to3> RQg2Qgg(&XQg2Qgg, 8*2, "./tables/R-Qg-Qgg.dat");

	return 0;
}
