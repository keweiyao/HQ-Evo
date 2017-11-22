#include "Langevin.h"
#include <vector>
#include <random>
#include <chrono>
#include <iostream>

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> white_noise(0.0, 1.0);
double A, B;
double const tiny = 1e-10;

double kperp_coeff(double E, double M, double T){
	return std::pow(T,3)*( A + B/(E*T) );
}

double kpara_coeff(double E, double M, double T){
	return std::pow(T,3)*( A + B/(E*T) );
}

void initialize_transport_coeff(double _A, double _B){
	A = _A; B = _B;
	std::cout << "A = " << A << ", B = " << B << std::endl;
};

void Langevin_step(	double E0, double M, double T, 
					double delta_t_lrf, std::vector<double> & pnew){
	pnew.resize(4);
	double pz0 = std::sqrt(E0*E0 - M*M + tiny); // in case M*M-M*M = 0-
	// step-1
	double kperp = kperp_coeff(E0, M, T),
		   kpara = kpara_coeff(E0, M, T);
	double drag = kpara/(2.*E0*T) - std::pow((std::sqrt(kpara)-std::sqrt(kperp))/pz0, 2);
		   
    double white_noise_holder[3];
    for (size_t i=0; i<3; ++i) white_noise_holder[i] = white_noise(generator);

	double perp_scale = std::sqrt(kperp*delta_t_lrf);
	double para_scale = std::sqrt(kpara*delta_t_lrf);

    pnew[1] = perp_scale * white_noise_holder[0];
    pnew[2] = perp_scale * white_noise_holder[1];
    pnew[3] = pz0 * (1. - drag * delta_t_lrf) + para_scale * white_noise_holder[2];
    pnew[0] = std::sqrt(M*M + pnew[1]*pnew[1] + pnew[2]*pnew[2] + pnew[3]*pnew[3]);

	// step-2
	kperp = kperp_coeff(pnew[0], M, T);
	kpara = kpara_coeff(pnew[0], M, T);
	drag = kpara/(2.*E0*T) - std::pow((std::sqrt(kpara)-std::sqrt(kperp))/pz0, 2);

	perp_scale = std::sqrt(kperp*delta_t_lrf);
	para_scale = std::sqrt(kpara*delta_t_lrf);

    pnew[1] = perp_scale * white_noise_holder[0];
    pnew[2] = perp_scale * white_noise_holder[1];
    pnew[3] = pz0 * (1. - drag * delta_t_lrf) + para_scale * white_noise_holder[2];
    pnew[0] = std::sqrt(M*M + pnew[1]*pnew[1] + pnew[2]*pnew[2] + pnew[3]*pnew[3]);
	return;
}

