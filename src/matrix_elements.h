#ifndef MATRIX_ELEMENTS_H
#define MATRIX_ELEMENTS_H
#include <cstdlib>

//======running coupling=======================================================
double alpha_s(double Q2);
//=============Baisc function for Q+q --> Q+q==================================
double M2_Qq2Qq(double t, void * params);
double dX_Qq2Qq_dt(double t, void * params);
double approx_XQq2Qq(double s, double Temp, double M);
//=============Baisc function for Q+g --> Q+g==================================
double M2_Qg2Qg(double t, void * params);
double dX_Qg2Qg_dt(double t, void * params);
double approx_XQg2Qg(double s, double Temp, double M);
//=============Baisc function for Q+q --> Q+q+g==================================
double M2_Qq2Qqg(double * x_, size_t n_dims_, void * params_);
//=============Baisc function for Q+g --> Q+g+g==================================
double M2_Qg2Qgg(double * x_, size_t n_dims_, void * params_);


#endif
