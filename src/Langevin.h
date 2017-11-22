#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <vector>

double kperp_coeff(double p, double M, double T);
double kpara_coeff(double p, double M, double T);
void Langevin_step(	double E0, double M, double T, 
					double delta_t_lrf, std::vector<double> & pnew);
void initialize_transport_coeff(double _A, double _B);
#endif
