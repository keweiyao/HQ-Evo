from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "../src/matrix_elements.h":
	cdef double dX_Qq2Qq_dPS(double * PS, size_t n_dims, void * params)
	cdef double approx_XQq2Qq(double s, double Temp, double M)
	cdef double dX_Qg2Qg_dPS(double * PS, size_t n_dims, void * params)
	cdef double approx_XQg2Qg(double s, double Temp, double M)
	cdef double M2_Qq2Qqg(double * x_, size_t n_dims_, void * params_)
	cdef double approx_XQq2Qqg(double s, double Temp, double M)
	cdef double M2_Qg2Qgg(double * x_, size_t n_dims_, void * params_)
	cdef double approx_XQg2Qgg(double s, double Temp, double M)

cdef extern from "../src/Xsection.h":
	cdef cppclass Xsection_2to2:
		Xsection_2to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double, double, double), double M1_, string name_)
		double calculate(double s, double Temp)
		void sample_dXdPS(double s, double Temp, vector[ vector[double] ] & final_states)
		double get_M1()
		double interpX(double s, double Temp)
	
	cdef cppclass Xsection_2to3:
		Xsection_2to3(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double, double, double), double M1_, string name_)
		double calculate(double s, double Temp)
		void sample_dXdPS(double s, double Temp, vector[ vector[double] ] & final_states)
		double get_M1()
		double interpX(double s, double Temp)

cdef extern from "../src/rates.h":
	cdef cppclass rates[T]:
		rates(T * Xprocess_, int degeneracy_, string name_)
		double calculate(double E1, double Temp)
		double interpR(double E1, double Temp)
		void sample_initial(double E1, double Temp, double &E2, double &s)

cdef Xsection_2to2 * x1 = new Xsection_2to2(&dX_Qq2Qq_dPS, &approx_XQq2Qq, 1.3, "./tables/X-Qq-Qq.dat")
cdef Xsection_2to2 * x2 = new Xsection_2to2(&dX_Qg2Qg_dPS, &approx_XQg2Qg, 1.3, "./tables/X-Qg-Qg.dat")
cdef Xsection_2to3 * x3 = new Xsection_2to3(&M2_Qq2Qqg, &approx_XQq2Qqg, 1.3, "./tables/X-Qq-Qqg.dat")
cdef Xsection_2to3 * x4 = new Xsection_2to3(&M2_Qg2Qgg, &approx_XQg2Qgg, 1.3, "./tables/X-Qg-Qgg.dat")

cdef rates[Xsection_2to2] * r1 = new rates[Xsection_2to2](x1, 3*4, "./tables/R-Qq-Qq.dat")
cdef rates[Xsection_2to2] * r2 = new rates[Xsection_2to2](x2, 8*2, "./tables/R-Qg-Qg.dat")
cdef rates[Xsection_2to3] * r3 = new rates[Xsection_2to3](x3, 3*4, "./tables/R-Qq-Qqg.dat")
cdef rates[Xsection_2to3] * r4 = new rates[Xsection_2to3](x4, 8*2, "./tables/R-Qg-Qgg.dat")
