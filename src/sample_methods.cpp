#include "sample_methods.h"
#include <random>
#include <fstream>
#include <algorithm>

void rejection_1d::build_interval(double xL, double xH, double fxL, double fxH){
	double xM = 0.5*(xL+xH);
	double fxM = f(xM, params);
	double mid_h = 0.5*(fxL+fxH);
	if (fxM < mid_h*0.1){
		build_interval(xL, xM, fxL, fxM);
		build_interval(xM, xH, fxM, fxH);
	}
	else{
		rectangle A; A.xL = xL; A.dx = xH-xL; A.fL = fxL; A.df = fxH-fxL; A.w = A.dx*mid_h;
		intervals.push_back(A);
		total_weight += A.w;
		return;
	}
}

double rejection_1d::sample(double (*f_) (double x, void * params), double xlo_, double xhi_, void * params_){
	f = f_;
	xlo = xlo_;
	xhi = xhi_;
	params = params_;
	intervals.clear();
	total_weight = 0.0;
			
	double fL = f(xlo, params), fH = f(xhi, params);
	build_interval(xlo, xhi, fL, fH);
	double cumulate = 0.0;
	for (auto&& ele : intervals){
		cumulate += ele.w/total_weight;
		ele.w = cumulate;
	}
	double r1, r2, xl, dx, fl, df, fh, lambda, xtry;
	size_t i;
	do{
		r1 = std::rand()*1./RAND_MAX;
		i=0;
		while (intervals[i].w < r1) i++;
		xl = intervals[i].xL; dx = intervals[i].dx;
		fl = intervals[i].fL; df = intervals[i].df;
		fh = fl+df;
		r2 = std::rand()*1./RAND_MAX;
		lambda = ( std::sqrt(fl*fl*(1.-r2) + fh*fh*r2) -fl )/df;
		xtry = xl + lambda*dx;
		}while (f(xtry, params)/(fl + df*lambda)*RAND_MAX < std::rand());
	return xtry;
}

double rejection_1d::plain_sample(double (*f_) (double x, void * params), double xlo_, double xhi_, void * params_){
	f = f_;
	xlo = xlo_;
	xhi = xhi_;
	params = params_;
	double dx = (xhi-xlo), xtry;
	double h = std::max(f(xlo, params), f(xhi, params));
	do{
		xtry = xlo + dx*std::rand()*1./RAND_MAX;
	}while (f(xtry, params)/h < std::rand()*1./RAND_MAX);
	return xtry;
}

// ----------Affine-invariant metropolis sample-------------------
AiMS::AiMS(void)
:	a(0.3), rd(), gen(rd()), sqrtZ(std::sqrt(1./a), std::sqrt(a)),
	reject(0.0, 1.0)
{
}
void AiMS::initialize(void){
    std::uniform_real_distribution<double> init_dis(-1, 1);
	for (auto&& w : walkers){
		for (size_t j=0; j < n_dims; ++j){
			w[j] = guess[j]*(1.0+0.1*init_dis(gen));
		}
	}
}
void AiMS::update(void){
	size_t ri;
	double sqz, z, Pnow, Ptry, Paccept;
	double * xtry = new double[n_dims];
    for (size_t i=0; i<Nwalker; ++i){
		do{ 
			ri = rd() % Nwalker;
		}while(i==ri);
		sqz = sqrtZ(gen);
		z = sqz*sqz;
		for (size_t j=0; j < n_dims; ++j) xtry[j] = walkers[ri][j] + z*(walkers[i][j] - walkers[ri][j]);
		Pnow = f(walkers[i], n_dims, params);
		Ptry = f(xtry, n_dims, params);
		Paccept = Ptry/Pnow*std::pow(z, n_dims-1);
		if (Paccept >= 1.0){
			for (size_t j=0; j < n_dims; ++j) buff_walkers[i][j] = xtry[j];
		}
		else if (Paccept >= reject(gen)){
			for (size_t j=0; j < n_dims; ++j) buff_walkers[i][j] = xtry[j];
		}
		else{
			for (size_t j=0; j < n_dims; ++j) buff_walkers[i][j] = walkers[i][j];
		}
	}
	for (size_t i=0; i<Nwalker; ++i){
		for (size_t j=0; j < n_dims; ++j) walkers[i][j] = buff_walkers[i][j];
	}
	delete[] xtry;
}
double AiMS::sample(double (*f_) (double*, size_t, void*), size_t n_dims_, void * params_, double * guess_){
	walkers.clear();
	f = f_; n_dims = n_dims_; params = params_; guess = guess_;
	Nwalker = n_dims*3;
	walkers.resize(Nwalker);
	buff_walkers.resize(Nwalker);
	for (auto&& w : walkers) w = new double[n_dims];
	for (auto&& w : buff_walkers) w = new double[n_dims];

	initialize();
	for (size_t i = 0; i<Nwalker*200; i++){
		update();
	}
	std::ofstream initf("samples.dat", std::ofstream::out | std::ofstream::app);
	for (size_t i = 0; i<100; i++){
		update();
		for (auto&& w : walkers){
			double k = w[0], p4 = w[1], phi4k = w[2], cos4 = w[3];
			double * p_ = static_cast<double*>(params);
			double s = p_[0], M = p_[2];
			double sqrts = std::sqrt(s);
			double cos_star = ((s-M*M)-2.*sqrts*(p4+k))/(2.*p4*k) +1.;
			double sin_star = std::sqrt(1. - cos_star*cos_star), sin4 = std::sqrt(1. - cos4*cos4);
			double cos_4k = std::cos(phi4k), sin_4k = std::sin(phi4k);
			double kx = k*(sin_star*cos_4k*cos4 - sin4*cos_star), ky = sin_star*sin_4k,
		   			kz = k*(sin_star*cos_4k*sin4 + cos4*cos_star);
			double Qx = -p4*sin4-kx, Qy = -ky, Qz = -kz-p4*cos4;
			double EQ = std::sqrt(Qx*Qx+Qy*Qy+Qz*Qz+M*M);
			initf << std::atan2(ky, kx) << " "
				  << 0.5*std::log((EQ+Qz)/(EQ-Qz)) <<  " "
				  << 0.5*std::log((k+kz)/(k-kz)) << std::endl;
		}
	}
	for (auto&& w : walkers) delete[] w;
	for (auto&& w : buff_walkers) delete[] w;

	return 1.0;
}














