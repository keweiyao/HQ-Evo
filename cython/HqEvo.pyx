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
		Xsection_2to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double * arg, double), double M1_, string name_, bool refresh)
		double calculate(double * arg)
		void sample_dXdPS(double * arg, vector[ vector[double] ] & FS)
		double interpX(double * arg)
	
	cdef cppclass Xsection_2to3:
		Xsection_2to3(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double * arg, double), double M1_, string name_, bool refresh)
		double calculate(double * arg)
		void sample_dXdPS(double * arg, vector[ vector[double] ] & FS)
		double interpX(double * arg)

	cdef cppclass f_3to2:
		f_3to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double * arg, double), double M1_, string name_, bool refresh)
		double calculate(double * arg)
		void sample_dXdPS(double * arg, vector[ vector[double] ] & FS)
		double interpX(double * arg)

cdef extern from "../src/rates.h":
	cdef cppclass rates_2to2:
		rates_2to2(Xsection_2to2 * Xprocess_, int degeneracy_, string name_, bool refresh)
		double calculate(double * arg)
		double interpR(double * arg)
		void sample_initial(double * arg, vector[ vector[double] ] & IS)

	cdef cppclass rates_2to3:
		rates_2to3(Xsection_2to3 * Xprocess_, int degeneracy_, string name_, bool refresh)
		double calculate(double * arg)
		double interpR(double * arg)
		void sample_initial(double * arg, vector[ vector[double] ] & IS)

	cdef cppclass rates_3to2:
		rates_3to2(f_3to2 * Xprocess_, int degeneracy_, string name_, bool refresh)
		double calculate(double * arg)
		double interpR(double * arg)
		void sample_initial(double * arg, vector[ vector[double] ] & IS)

#-------------------Wrap Xsection class---------------------------
cdef class pyX2to2:
	cdef Xsection_2to2 * cX2to2
	def __cinit__(self, channel, double mass, string filename, bool refresh):
		if channel == 'Qq->Qq':
			self.cX2to2 = new Xsection_2to2(&dX_Qq2Qq_dPS, &approx_XQq2Qq, mass, filename, refresh)
		elif channel == 'Qg->Qg':
			self.cX2to2 = new Xsection_2to2(&dX_Qg2Qg_dPS, &approx_XQg2Qg, mass, filename, refresh)
		else:
			raise ValueError("channel %s not implemented"%channel)
	cpdef double calculate(self, double s, double Temp):
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = s; arg[1] = Temp
		cdef double result = self.cX2to2.calculate(arg)
		free(arg)
		return result
	cpdef sample_dXdPS(self, double s, double Temp):
		cdef vector[ vector[double] ] FS
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = s; arg[1] = Temp
		self.cX2to2.sample_dXdPS(arg, FS)
		free(arg)
		return FS
	cpdef double interpX(self, double s, double Temp):
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = s; arg[1] = Temp
		cdef double result = self.cX2to2.interpX(arg)
		free(arg)
		return result

cdef class pyX2to3:
	cdef Xsection_2to3 * cX2to3 
	def __cinit__(self, channel, double mass, string filename, bool refresh):
		if channel == 'Qq->Qqg':
			self.cX2to3 = new Xsection_2to3(&M2_Qq2Qqg, &approx_XQq2Qqg, mass, filename, refresh)
		elif channel == 'Qg->Qgg':
			self.cX2to3 = new Xsection_2to3(&M2_Qg2Qgg, &approx_XQg2Qgg, mass, filename, refresh)
		else:
			raise ValueError("channel %s not implemented"%channel)
	cpdef double calculate(self, double & s, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = dt
		cdef double result = self.cX2to3.calculate(arg)
		free(arg)
		return result
	cpdef sample_dXdPS(self, double & s, double & Temp, double & dt):
		cdef vector[ vector[double] ] FS
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = dt
		self.cX2to3.sample_dXdPS(arg, FS)
		free(arg)
		return FS
	cpdef double interpX(self, double & s, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = dt
		cdef double result = self.cX2to3.interpX(arg)
		free(arg)
		return result

cdef class pyf3to2:
	cdef f_3to2 * cf3to2 
	def __cinit__(self, channel, double mass, string filename, bool refresh):
		if channel == 'Qqg->Qq':
			self.cf3to2  = new f_3to2(&Ker_Qqg2Qq, &approx_XQqg2Qq, mass, filename, refresh)
		elif channel == 'Qgg->Qg':
			self.cf3to2  = new f_3to2(&Ker_Qgg2Qg, &approx_XQgg2Qg, mass, filename, refresh)
		else:
			raise ValueError("channel %s not implemented"%channel)
	cpdef double calculate(self, double & s, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = dt
		cdef double result = self.cf3to2.calculate(arg)
		free(arg)
		return result
	cpdef sample_dXdPS(self, double & s, double & Temp, double & a1, double & a2):
		cdef vector[ vector[double] ] FS
		cdef double * arg = <double*>malloc(4*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = a1; arg[3] = a2;
		self.cf3to2.sample_dXdPS(arg, FS)
		free(arg)
		return FS
	cpdef double interpX(self, double & s, double & Temp, double & a1, double & a2, double & dt):
		cdef double * arg = <double*>malloc(5*sizeof(double))
		arg[0] = s; arg[1] = Temp; arg[2] = a1; arg[3] = a2; arg[4] = dt
		cdef double result = self.cf3to2.interpX(arg)
		free(arg)
		return result



#------------------Wrapper for Rates class-------------------------------
cdef class pyR2to2:
	cdef rates_2to2 * cR2to2
	def __cinit__(self, pyX2to2 x2to2, int degeneracy, string filename, bool refresh):
		self.cR2to2 = new rates_2to2(x2to2.cX2to2, degeneracy, filename, refresh)
	cpdef double calculate(self, double & E1, double & Temp):
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = E1; arg[1] = Temp
		cdef double result = self.cR2to2.calculate(arg)
		free(arg)
		return result
	cpdef double interpR(self, double & E1, double & Temp):
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = E1; arg[1] = Temp
		cdef double result =self.cR2to2.interpR(arg)
		free(arg)
		return result
	cpdef sample_initial(self, double & E1, double & Temp):
		cdef vector[ vector[double] ] IS
		cdef double * arg = <double*>malloc(2*sizeof(double))
		arg[0] = E1; arg[1] = Temp
		self.cR2to2.sample_initial(arg, IS)
		free(arg)
		return IS

cdef class pyR2to3:
	cdef rates_2to3 * cR2to3
	def __cinit__(self, pyX2to3 x2to3, int degeneracy, string filename, bool refresh):
		self.cR2to3 = new rates_2to3(x2to3.cX2to3, degeneracy, filename, refresh)
	cpdef double calculate(self, double E1, double Temp, double dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = Temp; arg[2] = dt
		cdef double result = self.cR2to3.calculate(arg)
		free(arg)
		return result
	cpdef double interpR(self, double & E1, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = Temp; arg[2] = dt
		cdef double result =self.cR2to3.interpR(arg)
		free(arg)
		return result
	cpdef sample_initial(self, double & E1, double & Temp, double & dt):
		cdef vector[ vector[double] ] IS
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = Temp; arg[2] = dt
		self.cR2to3.sample_initial(arg, IS)
		free(arg)
		return IS

cdef class pyR3to2:
	cdef rates_3to2 * cR3to2
	def __cinit__(self, pyf3to2 f3to2, int degeneracy, string filename, bool refresh):
		self.cR3to2 = new rates_3to2(f3to2.cf3to2, degeneracy, filename, refresh)
	cpdef double calculate(self, double & E1, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = Temp; arg[2] = dt
		cdef double result = self.cR3to2.calculate(arg)
		free(arg)
		return result
	cpdef double interpR(self, double & E1, double & Temp, double & dt):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = Temp; arg[2] = dt
		cdef double result =self.cR3to2.interpR(arg)
		free(arg)
		return result
	cpdef sample_initial(self, double & E1, double & Temp, double & dt):
		cdef vector[ vector[double] ] IS
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = Temp; arg[2] = dt
		self.cR3to2.sample_initial(arg, IS)
		free(arg)
		return IS
	
	
#-------------Heavy quark evolution class------------------------

cdef class HqEvo:
	cdef bool elastic, inelastic, detailed_balance
	cdef object X22list, X23list, X32list, R22list, R23list, R32list
	cdef size_t Nchannels, N22, N23, N32
	def __cinit__(self, mass=1.3, elastic=True, inelastic=False, detailed_balance=False, table_folder='./tables', refresh_table=False):
		self.elastic=elastic
		self.inelastic=inelastic
		self.detailed_balance=detailed_balance
		self.X22list = []
		self.X23list = []
		self.X32list = []
		self.R22list = []
		self.R23list = []
		self.R32list = []
		self.Nchannels = 0
		if not os.path.exists(table_folder):
			os.makedirs(table_folder)

		self.N22 = 0
		if self.elastic:
			xQq2Qq = pyX2to2('Qq->Qq', mass, "%s/XQq2Qq.hdf5"%table_folder, refresh_table)
			xQg2Qg = pyX2to2('Qg->Qg', mass, "%s/XQg2Qg.hdf5"%table_folder, refresh_table)
			rQq2Qq = pyR2to2(xQq2Qq, 12., "%s/RQq2Qq.hdf5"%table_folder, refresh_table)
			rQg2Qg = pyR2to2(xQg2Qg, 16., "%s/RQg2Qg.hdf5"%table_folder, refresh_table)
			self.X22list = [xQq2Qq, xQg2Qg]
			self.R22list = [rQq2Qq, rQg2Qg]
			self.N22 += len(self.X22list)
			self.Nchannels += self.N22
			
		self.N23 = 0
		if self.inelastic:
			xQq2Qqg = pyX2to3('Qq->Qqg', mass, "%s/XQq2Qqg.hdf5"%table_folder, refresh_table)
			xQg2Qgg = pyX2to3('Qg->Qgg', mass, "%s/XQg2Qgg.hdf5"%table_folder, refresh_table)
			rQq2Qqg = pyR2to3(xQq2Qqg, 12., "%s/RQq2Qqg.hdf5"%table_folder, refresh_table)
			rQg2Qgg = pyR2to3(xQg2Qgg, 16./2., "%s/RQg2Qgg.hdf5"%table_folder, refresh_table)
			self.X23list = [xQq2Qqg, xQg2Qgg]
			self.R23list = [rQq2Qqg, rQg2Qgg]
			self.N23 += len(self.X23list)
			self.Nchannels += self.N23

		self.N32 = 0
		if self.detailed_balance:
			xQqg2Qq = pyf3to2('Qqg->Qq', mass, "%s/XQqg2Qq.hdf5"%table_folder, refresh_table)
			xQgg2Qg = pyf3to2('Qgg->Qg', mass, "%s/XQgg2Qg.hdf5"%table_folder, refresh_table)
			rQqg2Qq = pyR3to2(xQqg2Qq, 12.*16., "%s/RQqg2Qq.hdf5"%table_folder, refresh_table)
			rQgg2Qg = pyR3to2(xQgg2Qg, 16.*16./2., "%s/RQgg2Qg.hdf5"%table_folder, refresh_table)
			self.X32list = [xQqg2Qq, xQgg2Qg]
			self.R32list = [rQqg2Qq, rQgg2Qg]
			self.N32 += len(self.X32list)
			self.Nchannels += self.N32
		print "Number of Channels", self.Nchannels

	cpdef sample_channel(self, double & E1, double & T, double & dt23, double & dt32):
		cdef double r, psum = 0.0, dt, Pmax = 0.1
		cdef size_t i=0
		cdef int channel_index = -1
		cdef double p[6]
		if self.elastic:
			for Rchannel in self.R22list:
				psum += Rchannel.interpR(E1, T)
				p[i] = psum
				i += 1
		if self.inelastic:
			for Rchannel in self.R23list:
				psum += Rchannel.interpR(E1, T, dt23) 
				# only use dt_from_last is not fully consistent
				p[i] = psum
				i += 1
		if self.detailed_balance:
			for Rchannel in self.R32list:
				psum += Rchannel.interpR(E1, T, dt32) 
				# only use dt_from_last is not fully consistent
				p[i] = psum
				i += 1
		#print E1, T, dt23, dt32
		#print p[0], p[1]-p[0], p[2]-p[1], p[3]-p[2], p[4]-p[3], p[5]-p[4]
		dt = Pmax/psum
		r = (<double>rand())/dt/RAND_MAX
		if r >= p[self.Nchannels-1]:
			return -1, dt
		for i in range(self.Nchannels):
			if r < p[i]:
				channel_index = i
				break
		return channel_index, dt

	cpdef sample_initial(self, int channel, double & E1, double & T, double & mean_dt23, double & mean_dt32):
		if channel < 0:
			raise ValueError("channel must be > 0, channel < 0 correspond to no scattering")
		elif channel < self.N22:
			return self.R22list[channel].sample_initial(E1, T)
		elif channel-self.N22 < self.N23:
			return self.R23list[channel-self.N22].sample_initial(E1, T, mean_dt23)
		elif channel-self.N22-self.N23 < self.N32:
			return self.R32list[channel-self.N22-self.N23].sample_initial(E1, T, mean_dt32)
		else:
			raise ValueError("channel %d not implememented"%channel)

	cpdef sample_final(self, int channel, double & s, double & T, double & mean_dt23, double & a1, double & a2):
		if channel < 0:
			raise ValueError("channel must be > 0, channel < 0 correspond to no scattering")
		elif channel < self.N22:
			return self.X22list[channel].sample_dXdPS(s, T)
		elif channel-self.N22 < self.N23:
			return self.X23list[channel-self.N22].sample_dXdPS(s, T, mean_dt23)
		elif channel-self.N22-self.N23 < self.N32:
			return self.X32list[channel-self.N22-self.N23].sample_dXdPS(s, T, a1, a2)
		else:
			raise ValueError("channel %d not implememented"%channel)
		
			


		
			

							
		










