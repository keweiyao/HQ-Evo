#include <iostream>
#include <fstream>
#include "sample_methods.h"
#include "matrix_elements.h"
#include "rates.h"
#include "utility.h"
#include <vector>
#include <cmath>

int main(){
	double M = 1.3;

	// Xsection: Q+q->Q+q
	Xsection_2to2 XQq2Qq(&dX_Qq2Qq_dPS, &approx_XQq2Qq, M, "./tables/X-Qq-Qq.dat");

	// Xsection: Q+g->Q+g
	//Xsection_2to2 XQg2Qg(&dX_Qg2Qg_dPS, &approx_XQg2Qg, M, "./tables/X-Qg-Qg.dat");

	// Xsection: Q+q->Q+q+g
	//Xsection_2to3 XQq2Qqg(&M2_Qq2Qqg, &approx_XQq2Qqg, M, "./tables/X-Qq-Qqg.dat");

	// Xsection: Q+g->Q+g+g
	//Xsection_2to3 XQg2Qgg(&M2_Qg2Qgg, &approx_XQg2Qgg, M, "./tables/X-Qg-Qgg.dat");

	// Rate: Q+q->Q+q
	rates<Xsection_2to2> RQq2Qq(&XQq2Qq, 3*4, "./tables/R-Qq-Qq.dat");

	// Rate: Q+g->Q+g
	//rates<Xsection_2to2> RQg2Qg(&XQg2Qg, 8*2, "./tables/R-Qg-Qg.dat");
	
	// Rate: Q+q->Q+q+g
	//rates<Xsection_2to3> RQq2Qqg(&XQq2Qqg, 3*4, "./tables/R-Qq-Qqg.dat");
	
	// Rate: Q+g->Q+g+g
	//rates<Xsection_2to3> RQg2Qgg(&XQg2Qgg, 8*2, "./tables/R-Qg-Qgg.dat");
	

	// Here is a sample code for a medium cell
	double temp=0.3;
	std::vector<double> u(4), v(3), p_lab(4), p_cell(4), p_cell_z(4), p_com(4);
	v[0] = 0.0; v[1] = 0.2; v[2] = 0.2;
	
	p_lab[0] = 10.0; p_lab[1] = 0.0; p_lab[2] = std::sqrt(10.0*10.0-M*M); p_lab[3] = 0.0;
	std::cout << "Heavy quark momentum in Lab frame" << std::endl;
	std::cout << p_lab[0] << "\t" << p_lab[1] << "\t" << p_lab[2] << "\t" << p_lab[3] << std::endl;
	
	// boost to cell frame 
	boost_by3(p_lab, p_cell, v);
	std::cout << "Heavy quark momentum in Cell frame" << std::endl;
	std::cout << p_cell[0] << "\t" << p_cell[1] << "\t" << p_cell[2] << "\t" << p_cell[3] << std::endl;
	
	// within z axis;
	// sample E2, s; given E1, temp 
	double s, E2;
	RQq2Qq.sample_initial(p_cell[0], temp, E2, s);
	std::cout << "E2 and s " << std::endl << E2 << "\t" << s << std::endl;
	double E1 = p_cell[0], p1 = std::sqrt(E1*E1-M*M);
	double costheta2 = (M*M + 2.*p_cell[0]*E2 - s)/2./p1/E2;
	double sintheta2 = std::sqrt(1. - costheta2*costheta2);
	std::cout << "cos2 "<< std::endl << costheta2 << std::endl;
	double phi = 2*M_PI*(std::rand()*1./RAND_MAX);
	double cosphi = std::cos(phi), sinphi = std::sin(phi);
	std::vector<double> p2(4);
	p2[0] = E2; p2[1] = E2*sintheta2*cosphi; p2[2] = E2*sintheta2*sinphi; p2[3] = E2*costheta2;
	std::cout << "p2"<< std::endl;
	std::cout << p2[0] << "\t" << p2[1] << "\t" << p2[2] << "\t" << p2[3] << std::endl;
	
	// rotate to back to Cell frame;
	double alpha = std::atan2(p_cell[2], p_cell[1])+M_PI/2.;
	double beta = std::atan2(std::sqrt(p_cell[1]*p_cell[1]+p_cell[2]*p_cell[2]), p_cell[3]);
	double gamma = 0.0;
	rotate_Euler(p_cell_z, p_cell, -gamma, -beta, -alpha);

	std::cout << p_cell[0] << "\t" << p_cell[1] << "\t" << p_cell[2] << "\t" << p_cell[3] << std::endl;

	
	return 0;
}
