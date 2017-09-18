#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip> // required to use std::setw

#include <stdlib.h>
#include <fstream>
#include <thread>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <boost/filesystem.hpp>

#include "Langevin.h"
#include "qhat.h"
#include "utility.h"

// two step Langevin: step 1, giving (temp, drag, kpara, kperp), return the first step (xi_Tx, xi_Ty, xi_L)
// E is the HQ energy, M is HQ mass, temp is the local temperature
// drag, kperp, kpara have the unit of [GeV^3]
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> distribution(0.0, 1.0);

double randUniform()
{
    return (double)rand()/RAND_MAX;
}

void Langevin_pre(double p_length, double M, double temp, double drag, double kperp, double kpara, double deltat_lrf, std::vector<double> &pre_result)
{
         // sample normal distribution with mean=0.0, width=1.0

        // save the values of rho as those will carry on to the post_point scheme second step
        std::vector<double> rho_xi(3, 0.0);
        for (size_t i=0; i< 3; ++i)
                rho_xi[i] = distribution(generator);


        std::vector<double> xi(3, 0.0);
        xi[0] = std::sqrt(kperp/deltat_lrf) * rho_xi[0];
        xi[1] = std::sqrt(kperp/deltat_lrf) * rho_xi[1];
        xi[2] = std::sqrt(kpara/deltat_lrf) * rho_xi[2];

        
        std::vector<double> new_p(4, 0.0);
        new_p[1] = xi[0] * deltat_lrf;
        new_p[2] = xi[1] * deltat_lrf;
        new_p[3] = p_length + (-drag * p_length + xi[2]) *deltat_lrf;
        new_p[0] = std::sqrt(M*M + new_p[1]*new_p[1] + new_p[2]*new_p[2] + new_p[3]*new_p[3]);

        pre_result={new_p[0], new_p[1], new_p[2], new_p[3], rho_xi[0], rho_xi[1], rho_xi[2]};


}


// for the two step post-point scheme, between step 1 (pre) and step 2(post), there are a few steps missing 
// the pre-step1 returns new_p, with new_p, we can update kpara and kperp 

void Langevin_post(double p_length, double M, double temp, double drag, double kperp, double kpara, double deltat_lrf, std::vector<double> const& pre_result, std::vector<double> & post_result)
{
// use the pre-point stored rho

        std::vector<double> rho={pre_result[4], pre_result[5], pre_result[6]};
        std::vector<double> xi(3, 0.0);
        xi[0] = std::sqrt(kperp/deltat_lrf) * rho[0];
        xi[1] = std::sqrt(kperp/deltat_lrf) * rho[1];
        xi[2] = std::sqrt(kpara/deltat_lrf) * rho[2];

		post_result.resize(4);
        post_result[1] = xi[0] * deltat_lrf;
        post_result[2] = xi[1] * deltat_lrf;
        post_result[3] = p_length + (-drag * p_length + xi[2]) * deltat_lrf;
        post_result[0] = std::sqrt(M*M + post_result[1]*post_result[1] + post_result[2]*post_result[2] + post_result[3]*post_result[3]);

}




//>>>>>>>>>>>>>>>>>>>>>> radiative energy loss from Higher Twist formula >>>>>>>>>>>>>>>>>>>>>>>
RadiationHT::RadiationHT(int flavor, std::string filename, bool plain, bool refresh)
:   HQflavor_(flavor),
    Ntime_(2), Ntemp_(11), NE_(10),
    timeL_(2), timeH_(4), tempL_(0.15), tempH_(0.45), EL_(10), EH_(100),
    dtime_((timeH_-timeL_)/(Ntime_-1.)), dtemp_((tempH_-tempL_)/(Ntemp_-1.)), dE_((EH_-EL_)/(NE_-1.)),
    RadTab_(boost::extents[Ntime_][Ntemp_][NE_]), Max_dNg_(boost::extents[Ntime_][Ntemp_][NE_])
{
    if (flavor == 4) HQmass_ = cMass;
    else if (flavor == 5) HQmass_ = bMass;
    else
    {
        std::cerr << "Flavor " << flavor << " is currently not implemented!" << std::endl;
        exit(EXIT_FAILURE);
    }

    bool file_exist = boost::filesystem::exists(filename);
    if (!file_exist || (file_exist && refresh))
    {
        std::cout << "Generating radiation table using HT formula... " << std::endl;
        std::cout << "threads " << std::thread::hardware_concurrency() << std::endl;
        std::vector<std::thread> threads;

        size_t Ncores = std::min(size_t(std::thread::hardware_concurrency()), NE_);
        size_t call_per_core = size_t(NE_ * 1./Ncores);
        size_t call_last_core = NE_ - call_per_core * (Ncores -1);

        auto code = [this](size_t Nstart_, size_t dN_) {this->tabulate(Nstart_, dN_);};

        for (size_t i=0; i<Ncores; ++i)
        {
            size_t Nstart = i*call_per_core;
            size_t dN = (i==Ncores-1) ? call_last_core : call_per_core;
            threads.push_back(std::thread(code, Nstart, dN));
        }

        for (std::thread& t: threads) t.join();
        saveToFile(filename, plain);
    } else {
        readFromFile(filename, plain);
    }
}



void RadiationHT::saveToFile(std::string filename, bool plain)
// save radiation table to plain .dat (plain) or hdf5 (!plain) file
{
    if (plain) {
        boost::filesystem::ofstream fRad{filename};
        fRad.precision(6);
        fRad << "#   time      temp      HQenergy    Ng/qhat   Max_dNgdxdydt/qhat" << std::endl;
        
        double time, temp, energy;
        for (size_t i=0; i<Ntime_; i++) {
          time = timeL_ + i*dtime_;
          for (size_t j=0; j<Ntemp_; j++) {
            temp = tempL_ + j*dtemp_;
            for (size_t k=0; k<NE_; k++) {
              energy = EL_ + k*dE_;
              fRad << std::setw(8) << time
                   << std::setw(8) << temp
                   << std::setw(8) << energy
                   << std::setw(15) << RadTab_[i][j][k]
                   << std::setw(15) << Max_dNg_[i][j][k] << std::endl;
            }
          }
        }

        std::cout << "Writing to file successfully :)" << std::endl;
    }
}

void RadiationHT::readFromFile(std::string filename, bool plain)
{
    if (plain){
        boost::filesystem::ifstream fRad{filename};
        
        double time, temp, energy;
        std::string header;
        std::getline(fRad, header);

        for (size_t i=0; i<Ntime_; i++)
          for (size_t j=0; j<Ntemp_; j++)
            for (size_t k=0; k<NE_; k++) 
              fRad >> time >> temp >> energy >> RadTab_[i][j][k] >> Max_dNg_[i][j][k];
    }
}



double RadiationHT::AlphaS(double kT, double temp)
// strong coupling constant, kT is the gluon transverse momentum, temp is medium temperature
{
    if (kT < M_PI*temp) kT = M_PI*temp;
    int nflavor;
    double lambdas;
    if (kT < cMass) {
        nflavor = 3;
        lambdas = 0.2;
    } else if (kT < bMass){
        nflavor = 4;
        lambdas = 0.172508;
    } else {
        nflavor = 5;
        lambdas = 0.130719;
    }

    double result = 4.*M_PI/((11.-2.*nflavor/3.)*(2.*log(kT/lambdas)));
    return result;
}


double RadiationHT::SplittingP(double x)
{
    return (2. - 2.*x + x*x) * CF / x;
}


double RadiationHT::TauF(double x, double y, double HQenergy)
{
    double k0_gluon = x*HQenergy;
    double kT_gluon = y*k0_gluon;
    double result= 2.*k0_gluon * (1. - x)/(pow(kT_gluon,2) + pow(x*HQmass_, 2));
    return result;
}


double RadiationHT::dNgOverDxdydt(double time, double temp, double HQenergy, double x, double y, double & max_dNg)
// HT formula (as a function of delta_time, temperature, HQenergy, x, y)
// arxiv: 1065.06477 Eqn.14
// returns actually dNg_over_dxdydt/qhat, I figured this is the best without all those hastles
{
    double cutoff = M_PI*temp;
    double k0_gluon = x*HQenergy;
    double kT_gluon = y*k0_gluon;
    double tauf = TauF(x, y, HQenergy);

    if (k0_gluon < cutoff)
        return 1e-12;
    if (tauf < 1./cutoff)
        return 1e-12;

    double result = 4./M_PI*Nc*AlphaS(kT_gluon, temp) * SplittingP(x) * pow(y, 5)
                    * pow(sin(time/(2*tauf*hbarC)), 2)
                    * pow(HQenergy*HQenergy/(y*y*HQenergy*HQenergy + HQmass_*HQmass_), 4)
                    / (k0_gluon * k0_gluon * hbarC);

    if (result > max_dNg)
        max_dNg = result;

    return result;
}


double wrapper_dNgdxdydt(double* args, size_t dim, void* params)
{
    (void) dim;
    struct gsl_NgIntegral_params *p = (struct gsl_NgIntegral_params*) params;
    double x = args[0];
    double y = args[1];

    double result = p->ptr->dNgOverDxdydt(p->time, p->temp, p->HQenergy, x, y, p->max_dNgdxdydt);
    return result;
}


void RadiationHT::calculate(double time, double temp, double HQenergy, double& result, double& max_dNg)
{
    double xl[2] = {0, 0};
    double xu[2] = {1, 1};

    struct gsl_NgIntegral_params ps;
    ps.time = time; 
    ps.temp = temp;
    ps.HQenergy = HQenergy;
    ps.max_dNgdxdydt = 0; // initialize the maxdNg before every integration
    ps.ptr = this;

    double res, err;
    const gsl_rng_type *T;
    gsl_rng* r;
    gsl_monte_function G;

    G.f = &wrapper_dNgdxdydt;
    G.dim = 2;
    G.params = &ps;

    size_t calls = 500000;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    {
        gsl_monte_vegas_state * s = gsl_monte_vegas_alloc(2);
        gsl_monte_vegas_integrate(&G, xl, xu, 2, 100000, r, s, &res, &err);
        int ctl = 0;

        do {
            gsl_monte_vegas_integrate(&G, xl, xu, 2, calls/5, r, s, &res, &err);
            ctl ++;
        } while(fabs(gsl_monte_vegas_chisq(s)-1.0)>0.1 && ctl <=20);

        gsl_monte_vegas_free(s);
    }

    gsl_rng_free(r);

    result = res;
    max_dNg = ps.max_dNgdxdydt;
}



void RadiationHT::tabulate(size_t NEstart, size_t dNE)
{
    double time, temp, energy;
    double result, max_dNg;
    for (size_t i=0; i<Ntime_; ++i) {
      time = timeL_ + i*dtime_;
      for (size_t j=0; j<Ntemp_; ++j) {
        temp = tempL_ + j*dtemp_;
        for (size_t k = NEstart; k<(NEstart + dNE); ++k) {
          energy = EL_ + k*dE_;
          calculate(time, temp, energy, result, max_dNg);
          RadTab_[i][j][k] = result;
          Max_dNg_[i][j][k] = max_dNg;
        }
      }
    }
}


double RadiationHT::interpR(double time, double temp, double HQenergy)
{
    if (time > timeH_) {
        time  = timeH_;
        std::cerr << "Attention! Tabulate in time is not big enough!" << std::endl;
    }

    if (temp > tempH_) {
        temp = tempH_;
        std::cerr << "Attention! Tabulate in temperature is not big enough!" << std::endl;
    }

    if (HQenergy > EH_) {
        HQenergy  = EH_;
        std::cerr << "Attention! Tabulate in energy is not big enough!" << std::endl;
    }

    double xTime, xTemp, xE, rTime, rTemp, rE;
    size_t iTime, iTemp, iE;

    xTime = (time-timeL_)/dtime_; iTime=floor(xTime+1e-6); rTime = xTime-iTime;
    xTemp = (temp-tempL_)/dtemp_; iTemp=floor(xTemp+1e-6); rTemp = xTemp-iTemp;
    xE  = (HQenergy - EL_)/dE_;   iE = floor(xE + 1e-6);   rE = xE-iE;

    return interpolate3d(&RadTab_, iTime, iTemp, iE, rTime, rTemp, rE);
}


double RadiationHT::interpMax(double time, double temp, double HQenergy)
{
    if (time > timeH_) {
        time  = timeH_;
        std::cerr << "Attention! Tabulate in time is not big enough!" << std::endl;
    }

    if (temp > tempH_) {
        temp = tempH_;
        std::cerr << "Attention! Tabulate in temperature is not big enough!" << std::endl;
    }

    if (HQenergy > EH_) {
        HQenergy  = EH_;
        std::cerr << "Attention! Tabulate in energy is not big enough!" << std::endl;
    }

    double xTime, xTemp, xE, rTime, rTemp, rE;
    size_t iTime, iTemp, iE;

    xTime = (time-timeL_)/dtime_; iTime=floor(xTime+1e-6); rTime = xTime-iTime;
    xTemp = (temp-tempL_)/dtemp_; iTemp=floor(xTemp+1e-6); rTemp = xTemp-iTemp;
    xE  = (HQenergy - EL_)/dE_;   iE = floor(xE + 1e-6);   rE = xE-iE;

    return interpolate3d(&Max_dNg_, iTime, iTemp, iE, rTime, rTemp, rE);
}


double RadiationHT::getNg(double time, double temp, double HQenergy, double qhat)
//return the average number of radiated gluon from a hard parton with (E,T) in unit time
{
    double delta_ng = interpR(time, temp, HQenergy);
    return delta_ng * qhat;
}


double RadiationHT::getMaxdNg(double time, double temp, double HQenergy)
{
    double max_ = interpMax(time, temp, HQenergy);
    return max_;
}

bool RadiationHT::emitGluon(double time, double temp, double HQenergy, double qhat, double deltat)
// check with a gluon can be emitted in this timestep
{
    double delta_ng = getNg(time, temp, HQenergy, qhat) * deltat;
    if (delta_ng > 1)
        std::cerr << "Error! Gluon emission probability exceeds 1! " << std::endl;

    double xLow = M_PI*temp/HQenergy;
    double xHigh = 1.0;
    double xRange = xHigh - xLow;

    if ((xRange > eSMALL) && (randUniform() < delta_ng - eSMALL) && (2.*HQenergy*(HQenergy-M_PI*temp) > HQmass_*HQmass_))
        return true;
    else
        return false;
}


bool RadiationHT::sampleGluon(double time, double temp, double HQenergy, double qhat, double deltat, std::vector<double> & gluon)
// given delta_time, temp, HQenergy, qhat, deltat, sample emitted gluon momentum
// the gluon momentum should be rotated into HQ 
{
    gluon.resize(3);
    gluon[0] = 0.; gluon[1] = 0.; gluon[2] = 0.;

    double cutoff = M_PI*temp;
    double xLow = cutoff/HQenergy;
    double xHigh = 1.;
    double xRange = xHigh - xLow;
    double maxF = getMaxdNg(time, temp, HQenergy);

    double x, y, tauf, dNg_dxdydt, dum;
    int count = 0;
    
    do {
        x = xLow + randUniform()*xRange;
        y = randUniform();
        tauf = TauF(x, y, HQenergy);
        dNg_dxdydt = dNgOverDxdydt(time, temp, HQenergy, x, y, dum);
        count ++;

        if (count > 1e5)
        {
            std::cout << "Given up loop at point 1 ... " << std::endl;
            return false;
        }
    } while((tauf < 1./cutoff) || (maxF * randUniform() > dNg_dxdydt));

// <<<< successfully sampled Gluon <<<<<
    double kT_gluon = x*y*HQenergy;
    double theta = 2.*M_PI * randUniform();
    gluon[0] = kT_gluon*cos(theta);
    gluon[1] = kT_gluon*sin(theta);
    gluon[2] = x * HQenergy * sqrt(1.-y*y);


// >>>>> remove the momentum broadenning from gluons
// in SCao's Langevin, using senario narrow 4
    int count2 = 0;
    double width_gluon = sqrt(Nc * qhat * tauf / 2.);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::normal_distribution<double> distribution(0, width_gluon);
    std::vector<double> dk_gluon = {0., 0., 0.};
    double dkT_;
    do {
        dk_gluon[0] = distribution(generator);
        dk_gluon[1] = distribution(generator);
        dkT_ = sqrt(dk_gluon[0] * dk_gluon[0] + dk_gluon[1]*dk_gluon[1]);
        count2 ++;

        if (count2 > 1e6)
        {
            std::cout << "Given up loop at point 2 ... " << dk_gluon[0] << " " << dk_gluon[1] << std::endl;
            gluon[0] = 0.; gluon[1] = 0.; gluon[2] = 0.;
            return false;
        }
    } while((dkT_ > kT_gluon) || (dkT_*dkT_ - 2.*dk_gluon[0]*gluon[0] - 2*dk_gluon[1]*gluon[1])>0 );


    gluon[0] = gluon[0] - dk_gluon[0];
    gluon[1] = gluon[1] - dk_gluon[1];

    return true;
}
