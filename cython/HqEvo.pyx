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
	cdef double approx_XQq2Qq(double * arg, double M)
	cdef double dX_Qg2Qg_dPS(double * PS, size_t n_dims, void * params)
	cdef double approx_XQg2Qg(double * arg, double M)

	cdef double M2_Qq2Qqg(double * x_, size_t n_dims_, void * params_)
	cdef double approx_XQq2Qqg(double * arg, double M)
	cdef double M2_Qg2Qgg(double * x_, size_t n_dims_, void * params_)
	cdef double approx_XQg2Qgg(double * arg, double M)

	cdef double Ker_Qqg2Qq(double * x_, size_t n_dims_, void * params_)
	cdef double approx_XQqg2Qq(double * arg, double M)
	cdef double Ker_Qgg2Qg(double * x_, size_t n_dims_, void * params_)
	cdef double approx_XQgg2Qg(double * arg, double M)
	

cdef extern from "../src/Xsection.h":
	cdef cppclass Xsection_2to2:
		Xsection_2to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double * arg, double), double M1_, string name_)
		double calculate(double * arg)
		void sample_dXdPS(double * arg, vector[ vector[double] ] & final_states)
		double interpX(double * arg)
	
	cdef cppclass Xsection_2to3:
		Xsection_2to3(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double * arg, double), double M1_, string name_)
		double calculate(double * arg)
		void sample_dXdPS(double * arg, vector[ vector[double] ] & final_states)
		double interpX(double * arg)

	cdef cppclass f_3to2:
		f_3to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double * arg, double), double M1_, string name_)
		double calculate(double * arg)
		void sample_dXdPS(double * arg, vector[ vector[double] ] & final_states)
		double interpX(double * arg)

cdef extern from "../src/rates.h":
	cdef cppclass rates_2to2:
		rates_2to2(Xsection_2to2 * Xprocess_, int degeneracy_, string name_)
		double calculate(double * arg)
		double interpR(double * arg)
		vector[ vector[double] ] sample_initial(double * arg_in)

	cdef cppclass rates_2to3:
		rates_2to3(Xsection_2to3 * Xprocess_, int degeneracy_, string name_)
		double calculate(double * arg)
		double interpR(double * arg)
		vector[ vector[double] ] sample_initial(double * arg_in)

	cdef cppclass rates_3to2:
		rates_3to2(f_3to2 * Xprocess_, int degeneracy_, string name_)
		double calculate(double * arg)
		double interpR(double * arg)
		vector[ vector[double] ] sample_initial(double * arg_in)

#-------------------Wrap Xsection class---------------------------
cdef class pyX2to2:
	cdef Xsection_2to2 * cX2to2
	def __cinit__(self, channel, double mass, string filename):
		if channel == 'Qq->Qq':
			self.cX2to2 = new Xsection_2to2(&dX_Qq2Qq_dPS, &approx_XQq2Qq, mass, filename)
		elif channel == 'Qg->Qg':
			self.cX2to2 = new Xsection_2to2(&dX_Qg2Qg_dPS, &approx_XQg2Qg, mass, filename)
		else:
			raise ValueError("channel %s not implemented"%channel)
	cpdef double calculate(self, double & s, double & Temp):
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = s; arg[1] = Temp
		return self.cX2to2.calculate(arg)
	cpdef sample_dXdPS(self, double & s, double & Temp):
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = s; arg[1] = Temp
		cdef vector[ vector[double] ] final_states
		self.cX2to2.sample_dXdPS(arg, final_states)
		return final_states
	cpdef double interpX(self, double & s, double & Temp):
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = s; arg[1] = Temp
		return self.cX2to2.interpX(arg)

cdef class pyX2to3:
	cdef Xsection_2to3 * cX2to3 
	def __cinit__(self, channel, double mass, string filename):
		if channel == 'Qq->Qqg':
			self.cX2to3 = new Xsection_2to3(&M2_Qq2Qqg, &approx_XQq2Qqg, mass, filename)
		elif channel == 'Qg->Qgg':
			self.cX2to3 = new Xsection_2to3(&M2_Qg2Qgg, &approx_XQg2Qgg, mass, filename)
		else:
			raise ValueError("channel %s not implemented"%channel)
	cpdef double calculate(self, double & s, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = dt
		return self.cX2to3.calculate(arg)
	cpdef sample_dXdPS(self, double & s, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = dt
		cdef vector[ vector[double] ] final_states
		self.cX2to3.sample_dXdPS(arg, final_states)
		return final_states
	cpdef double interpX(self, double & s, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = dt
		return self.cX2to3.interpX(arg)

cdef class pyf3to2:
	cdef f_3to2 * cf3to2 
	def __cinit__(self, channel, double mass, string filename):
		if channel == 'Qqg->Qq':
			self.cf3to2  = new f_3to2(&Ker_Qqg2Qq, &approx_XQqg2Qq, mass, filename)
		elif channel == 'Qgg->Qg':
			self.cf3to2  = new f_3to2(&Ker_Qgg2Qg, &approx_XQgg2Qg, mass, filename)
		else:
			raise ValueError("channel %s not implemented"%channel)
	cpdef double calculate(self, double & s, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = dt
		return self.cf3to2.calculate(arg)
	cpdef sample_dXdPS(self, double & s, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = dt
		cdef vector[ vector[double] ] final_states
		self.cf3to2.sample_dXdPS(arg, final_states)
		return final_states
	cpdef double interpX(self, double & s, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = dt
		return self.cf3to2.interpX(arg)



#------------------Wrapper for Rates class-------------------------------
cdef class pyR2to2:
	cdef rates_2to2 * cR2to2
	def __cinit__(self, pyX2to2 x2to2, int degeneracy, string filename):
		self.cR2to2 = new rates_2to2(x2to2.cX2to2, degeneracy, filename)
	cpdef double calculate(self, double & E1, double & Temp):
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = E1; arg[1] = Temp
		return self.cR2to2.calculate(arg)
	cpdef double interpR(self, double & E1, double & Temp):
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = E1; arg[1] = Temp
		return self.cR2to2.interpR(arg)
	cpdef sample_initial(self, double & E1, double & Temp):
		cdef double * arg_in = <double*>malloc(2*sizeof(double))
		arg_in[0] = E1; arg_in[1] = Temp
		return self.cR2to2.sample_initial(arg_in)


cdef class pyR2to3:
	cdef rates_2to3 * cR2to3
	def __cinit__(self, pyX2to3 x2to3, int degeneracy, string filename):
		self.cR2to3 = new rates_2to3(x2to3.cX2to3, degeneracy, filename)
	cpdef double calculate(self, double & E1, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = Temp; arg[2] = dt
		return self.cR2to3.calculate(arg)
	cpdef double interpR(self, double & E1, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = Temp; arg[2] = dt
		return self.cR2to3.interpR(arg)
	cpdef sample_initial(self, double & E1, double & Temp, double & dt):
		cdef double * arg_in = <double*>malloc(3*sizeof(double))
		arg_in[0] = E1; arg_in[1] = Temp; arg_in[2] = dt
		return self.cR2to3.sample_initial(arg_in)

cdef class pyR3to2:
	cdef rates_3to2 * cR3to2
	def __cinit__(self, pyf3to2 f3to2, int degeneracy, string filename):
		self.cR3to2 = new rates_3to2(f3to2.cf3to2, degeneracy, filename)
	cpdef double calculate(self, double & E1, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = Temp; arg[2] = dt
		return self.cR3to2.calculate(arg)
	cpdef double interpR(self, double & E1, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = Temp; arg[2] = dt
		return self.cR3to2.interpR(arg)
	cpdef sample_initial(self, double & E1, double & Temp, double & dt):
		cdef double * arg_in = <double*>malloc(3*sizeof(double))
		arg_in[0] = E1; arg_in[1] = Temp; arg_in[2] = dt
		return self.cR3to2.sample_initial(arg_in)
	
	
#-------------Heavy quark evolution class------------------------

cdef class HqEvo:
	cdef bool elastic, inelastic, detailed_balance
	cdef object X22list, X23list, R22list, R23list
	cdef size_t Nchannels, Nelastic, Ninelastic
	def __cinit__(self, mass=1.3, elastic=True, inelastic=False, detailed_balance=False, table_folder='./tables'):
		self.elastic=elastic
		self.inelastic=inelastic
		self.detailed_balance=detailed_balance
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
			rQq2Qq = pyR2to2(xQq2Qq, 12, "%s/RQq2Qq.dat"%table_folder)
			rQg2Qg = pyR2to2(xQg2Qg, 16, "%s/RQg2Qg.dat"%table_folder)
			self.X22list = [xQq2Qq, xQg2Qg]
			self.R22list = [rQq2Qq, rQg2Qg]
			self.Nelastic = len(self.X22list)
			self.Nchannels += self.Nelastic
		else:
			self.Nelastic = 0
		
		if self.inelastic:
			xQq2Qqg = pyX2to3('Qq->Qqg', mass, "%s/XQq2Qqg.dat"%table_folder)
			xQg2Qgg = pyX2to3('Qg->Qgg', mass, "%s/XQg2Qgg.dat"%table_folder)
			rQq2Qqg = pyR2to3(xQq2Qqg, 12, "%s/RQq2Qqg.dat"%table_folder)
			rQg2Qgg = pyR2to3(xQg2Qgg, 16/2., "%s/RQg2Qgg.dat"%table_folder)
			self.X23list = [xQq2Qqg, xQg2Qgg]
			self.R23list = [rQq2Qqg, rQg2Qgg]
			self.Ninelastic = len(self.X23list)
			self.Nchannels += self.Ninelastic
		else:
			self.Ninelastic = 0

		if self.detailed_balance:
			xQqg2Qq = pyf3to2('Qqg->Qq', mass, "%s/XQqg2Qq.dat"%table_folder)
			xQgg2Qg = pyf3to2('Qgg->Qg', mass, "%s/XQgg2Qg.dat"%table_folder)
			rQqg2Qq = pyR3to2(xQqg2Qq, 12*16, "%s/RQqg2Qq.dat"%table_folder)
			rQgg2Qg = pyR3to2(xQgg2Qg, 16*16/2., "%s/RQgg2Qg.dat"%table_folder)

	cpdef sample_channel(self, double & E1, double & T, double & dt_from_last, double & Nc):
		cdef double r, psum = 0.0, dt, Pmax = 0.1
		cdef size_t i=0, channel_index=-1
		cdef double * p = <double*>malloc(self.Nchannels*sizeof(double))
		if self.elastic:
			for Rchannel in self.R22list:
				psum += Rchannel.interpR(E1, T)
				p[i] = psum
				i += 1
		
		if self.inelastic:
			for Rchannel in self.R23list:
				psum += Nc*Rchannel.interpR(E1, T, dt_from_last) # only use dt_from_last is not fully consistent
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

	cpdef sample_initial(self, int channel, double & E1, double & T, double & t_elapse_cell):
		if channel < 0:
			raise ValueError("channel must be > 0, channel < 0 correspond to no scattering")
		elif channel < self.Nelastic:
			#print "An elastic scattering"
			return self.R22list[channel].sample_initial(E1, T)
		elif channel-self.Nelastic < self.Ninelastic:
			#print "An inelastic scattering"
			return self.R23list[channel-self.Nelastic].sample_initial(E1, T, t_elapse_cell)
		else:
			raise ValueError("channel %d not implememented"%channel)

	cpdef sample_final(self, int channel, double & s, double & T, double & t_elapse_com):
		if channel < 0:
			raise ValueError("channel must be > 0, channel < 0 correspond to no scattering")
		elif channel < self.Nelastic:
			#print "An elastic scattering"
			return self.X22list[channel].sample_dXdPS(s, T)
		elif channel-self.Nelastic < self.Ninelastic:
			#print "An inelastic scattering"
			return self.X23list[channel-self.Nelastic].sample_dXdPS(s, T, t_elapse_com)
		else:
			raise ValueError("channel %d not implememented"%channel)
		
			


		
			

							
		










