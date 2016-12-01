from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdlib cimport malloc, free
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport fmin
import os

#------------------Import C++ fucntions and class for Xsection and rates------------------
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
		double interpX(double s, double Temp)
	
	cdef cppclass Xsection_2to3:
		Xsection_2to3(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double, double, double), double M1_, string name_)
		double calculate(double s, double Temp)
		void sample_dXdPS(double s, double Temp, vector[ vector[double] ] & final_states)
		double interpX(double s, double Temp)

cdef extern from "../src/rates.h":
	cdef cppclass rates[T]:
		rates(T * Xprocess_, int degeneracy_, string name_)
		double calculate(double E1, double Temp)
		double interpR(double E1, double Temp)
		void sample_initial(double E1, double Temp, double &E2, double &s)

#-------------------Wrap Xsection and rate class---------------------------
cdef class pyX2to2:
	cdef Xsection_2to2 * cX2to2
	def __cinit__(self, channel, double mass, string filename):
		if channel == 'Qq->Qq':
			self.cX2to2 = new Xsection_2to2(&dX_Qq2Qq_dPS, &approx_XQq2Qq, mass, filename)
		elif channel == 'Qg->Qg':
			self.cX2to2 = new Xsection_2to2(&dX_Qg2Qg_dPS, &approx_XQg2Qg, mass, filename)
		else:
			raise ValueError("channel %s not implemented"%channel)
	cpdef double calculate(self, double s, double Temp):
		return self.cX2to2.calculate(s, Temp)
	cpdef sample_dXdPS(self, double s, double Temp):
		cdef vector[ vector[double] ] final_states
		self.cX2to2.sample_dXdPS(s, Temp, final_states)
		return final_states
	cpdef double interpX(self, double s, double Temp):
		return self.cX2to2.interpX(s, Temp)

cdef class pyX2to3:
	cdef Xsection_2to3 * cX2to3 
	def __cinit__(self, channel, double mass, string filename):
		if channel == 'Qq->Qqg':
			self.cX2to3 = new Xsection_2to3(&M2_Qq2Qqg, &approx_XQq2Qqg, mass, filename)
		elif channel == 'Qg->Qgg':
			self.cX2to3 = new Xsection_2to3(&M2_Qg2Qgg, &approx_XQg2Qgg, mass, filename)
		else:
			raise ValueError("channel %s not implemented"%channel)
	cpdef double calculate(self, double s, double Temp):
		return self.cX2to3.calculate(s, Temp)
	cpdef sample_dXdPS(self, double s, double Temp):
		cdef vector[ vector[double] ] final_states
		self.cX2to3.sample_dXdPS(s, Temp, final_states)
		return final_states
	cpdef double interpX(self, double s, double Temp):
		return self.cX2to3.interpX(s, Temp)

cdef class pyR2to2:
	cdef rates[Xsection_2to2] * cR2to2
	def __cinit__(self, pyX2to2 x2to2, int degeneracy, string filename):
		self.cR2to2 = new rates[Xsection_2to2](x2to2.cX2to2, degeneracy, filename)
	cpdef double calculate(self, double E1, double Temp):
		return self.cR2to2.calculate(E1, Temp)
	cpdef double interpR(self, double E1, double Temp):
		return self.cR2to2.interpR(E1, Temp)
	cpdef sample_initial(self, double E1, double Temp):
		cdef double E2 = 0.0, s = 0.0
		self.cR2to2.sample_initial(E1, Temp, E2, s)
		return E2, s

cdef class pyR2to3:
	cdef rates[Xsection_2to3] * cR2to3
	def __cinit__(self, pyX2to3 x2to3, int degeneracy, string filename):
		self.cR2to3 = new rates[Xsection_2to3](x2to3.cX2to3, degeneracy, filename)
	cpdef double calculate(self, double E1, double Temp):
		return self.cR2to3.calculate(E1, Temp)
	cpdef double interpR(self, double E1, double Temp):
		return self.cR2to3.interpR(E1, Temp)
	cpdef sample_initial(self, double E1, double Temp):
		cdef double E2 = 0.0, s = 0.0
		self.cR2to3.sample_initial(E1, Temp, E2, s)
		return E2, s
	
	
#-------------Heavy quark evolution class------------------------

cdef class HqEvo:
	cdef bool elastic, inelastic
	cdef object X22list, X23list, R22list, R23list
	cdef size_t Nchannels, Nelastic, Ninelastic
	def __cinit__(self, mass=1.3, elastic=True, inelastic=False, table_folder='./tables'):
		self.elastic=elastic
		self.inelastic=inelastic
		self.X22list = []
		self.X23list = []
		self.R22list = []
		self.R23list = []
		self.Nchannels = 0
		if not os.path.exists(table_folder):
			os.makedirs(table_folder)
		if self.elastic:
			xQq2Qq = pyX2to2('Qq->Qq', mass, "%s/XQq2Qq.dat"%table_folder)
			xQg2Qg = pyX2to2('Qg->Qg', mass, "%s/XQg2Qg.dat"%table_folder)
			rQq2Qq = pyR2to2(xQq2Qq, 3*4, "%s/RQq2Qq.dat"%table_folder)
			rQg2Qg = pyR2to2(xQg2Qg, 8*2, "%s/RQg2Qg.dat"%table_folder)
			self.X22list = [xQq2Qq, xQg2Qg]
			self.R22list = [rQq2Qq, rQg2Qg]
			self.Nelastic = len(self.X22list)
			self.Nchannels += self.Nelastic
		if self.inelastic:
			xQq2Qqg = pyX2to3('Qq->Qqg', mass, "%s/XQq2Qqg.dat"%table_folder)
			xQg2Qgg = pyX2to3('Qg->Qgg', mass, "%s/XQg2Qgg.dat"%table_folder)
			rQq2Qqg = pyR2to3(xQq2Qqg, 3*4, "%s/RQq2Qqg.dat"%table_folder)
			rQg2Qgg = pyR2to3(xQg2Qgg, 8*2, "%s/RQg2Qgg.dat"%table_folder)
			self.X23list = [xQq2Qqg, xQg2Qgg]
			self.R23list = [rQq2Qqg, rQg2Qgg]
			self.Ninelastic = len(self.X23list)
			self.Nchannels += self.Ninelastic 

	cpdef sample_channel(self, double E1, double T):
		cdef double r, psum = 0.0, dt, Pmax = 0.1
		cdef int i=0, channel_index
		cdef double * p = <double*>malloc(self.Nchannels*sizeof(double))
		if self.elastic:
			for Rchannel in self.R22list:
				psum += Rchannel.interpR(E1, T)
				p[i] = psum
				i += 1
		if self.inelastic:
			for Rchannel in self.R23list:
				psum += Rchannel.interpR(E1, T)
				p[i] = psum
				i += 1

		dt = Pmax/psum
		r = (<double>rand())/dt/RAND_MAX
		if r >= p[self.Nchannels-1]:
			return -1, dt
		for i in range(self.Nchannels):
			if r < p[i]:
				channel_index = i
				break
		return channel_index, dt

	cpdef sample_initial(self, int channel, double E1, double T):
		cdef double E2=0.0, s=0.0
		if channel < 0:
			raise ValueError("channel must be > 0, channel < 0 correspond to no scattering")
		elif channel < self.Nelastic:
			#print "An elastic scattering"
			return self.R22list[channel].sample_initial(E1, T)
		elif channel-self.Nelastic < self.Ninelastic:
			#print "An inelastic scattering"
			return self.R23list[channel-self.Nelastic].sample_initial(E1, T)
		else:
			raise ValueError("channel %d not implememented"%channel)

	cpdef sample_final(self, int channel, double s, double T):
		cdef vector[ vector[double] ] results 
		if channel < 0:
			raise ValueError("channel must be > 0, channel < 0 correspond to no scattering")
		elif channel < self.Nelastic:
			#print "An elastic scattering"
			results = self.X22list[channel].sample_dXdPS(s, T)
			return results
		elif channel-self.Nelastic < self.Ninelastic:
			#print "An inelastic scattering"
			results = self.X23list[channel-self.Nelastic].sample_dXdPS(s, T)
			return results
		else:
			raise ValueError("channel %d not implememented"%channel)
		
			


		
			

							
		










