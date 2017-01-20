#ifndef CONSTANTS_H
#define CONSTATNS_H

#include <cmath>
#include <vector>
#include <boost/multi_array.hpp>
//=============constants=======================================================
const double c4d9 = 4./9.;
const double c1d9 = 1./9.;
const double c16pi = 16.*M_PI;
const double c48pi = 48.*M_PI;
const double c16pi2 = 16.*M_PI*M_PI;
const double c64d9pi2 = 64./9.*M_PI*M_PI;
const double c256pi4 = 256.*std::pow(M_PI, 4);
const int Nc = 3, nf = 3;
const double pf_g = 4.*M_PI/3.*(Nc + nf/2.); // prefractor for gluon self energy^2 
const double pf_q = M_PI/2.*(Nc*Nc - 1)/2./Nc; // prefractor for quark self energy^2 
const double alpha0 = 4.*M_PI/(11. - 2./3.*nf); // alpha_s(Q2 = e*Lambda2)
const double Lambda2 = 0.2*0.2; // [GeV^2] Lambda QCD squared
const double Q2cut_l = -Lambda2*exp(alpha0), Q2cut_h = 0.; // [GeV^2] ranges within which alphas > 1 and will be cut


double interpolate2d(	boost::multi_array<double, 2> * A, 
					 	const int& ni, const int& nj, 
					 	const double& ri, const double& rj);
double interpolate3d(	boost::multi_array<double, 3> * A, 
						const int& ni, const int& nj, const int& nk,
						const double& ri, const double& rj, const double& rk);
double interpolate4d(	boost::multi_array<double, 4> * A, 
						const int& ni, const int& nj, const int& nk, const int& nt, 
						const double& ri, const double& rj, const double& rk, const double& rt);

#endif
