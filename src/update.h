#ifndef UPDATE_H
#define UPDATE_H

#include "Xsection.h"
#include "rates.h"
#include <random>

class update{
private:
	std::random_device rd;
    std::mt19937 gen;
	std::uniform_real_distribution<double> dist01;
	Xsection_2to2 XQq2Qq;
	Xsection_2to2 XQg2Qg;
	Xsection_2to3 XQq2Qqg;
	Xsection_2to3 XQg2Qgg;
	rates<Xsection_2to2> RQq2Qq;
	rates<Xsection_2to2> RQg2Qg;
	rates<Xsection_2to3> RQq2Qqg;
	rates<Xsection_2to3> RQg2Qgg;
	
public:
	update(double M);
	int sample_channel(double E1, double T, double dt);
	
};

#endif


