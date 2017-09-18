#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <vector>
#include <boost/multi_array.hpp>
#include <string>
#include <random>

#include "qhat.h"

void Langevin_pre(double p_length, double M, double temp, double drag, double kperp, double kpara, double delta_lrf, std::vector<double> & pre_result);

void Langevin_post(double p_length, double M, double temp, double drag, double kperp, double kpara, double delta_lrf, std::vector<double> const& pre_result, std::vector<double> & post_result);


// this is the contribution from radiative energy loss 
class RadiationHT
{
private:
    int HQflavor_;
    double HQmass_;
    size_t Ntime_, Ntemp_, NE_;
    double timeL_, timeH_, tempL_, tempH_, EL_, EH_;
    double dtime_, dtemp_, dE_;

    boost::multi_array<double, 3> RadTab_, Max_dNg_;

    void saveToFile(std::string filename, bool plain);
    void readFromFile(std::string filename, bool plain);

    double AlphaS(double kT, double temp); // a fancy running coupling (with flavor dependence)
    double SplittingP(double x); // splitting function for quarks
    double TauF(double x, double y, double HQenergy); // formation time for emitted gluon
    void tabulate(size_t NEstart, size_t dNE);
    double interpR(double time, double temp, double HQenergy);
    double interpMax(double time, double temp, double HQenergy);
    double getMaxdNg(double time, double temp, double HQenergy);
    void calculate(double time, double temp, double HQenergy, double& result, double& max_dNg);
 
public:
    RadiationHT(double mass, std::string filename, bool plain, bool refresh);
    double dNgOverDxdydt(double time, double temp, double HQenergy, double x, double y, double& max_dNg);
    double getNg(double time, double temp, double HQenergy, double qhat);
    bool emitGluon(double time, double temp, double HQenergy, double qhat, double deltat);
    bool sampleGluon(double time, double temp, double HQenergy, double qhat, double deltat, std::vector<double> & gluon);
};

struct gsl_NgIntegral_params
{
    double time;
    double temp;
    double HQenergy;
    double max_dNgdxdydt;
    RadiationHT * ptr;
};



#endif
