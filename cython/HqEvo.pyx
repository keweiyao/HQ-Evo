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
		void sample_dXdPS(double * arg, vector[ vector[double] ] & FS)
	
	cdef cppclass Xsection_2to3:
		Xsection_2to3(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double * arg, double), double M1_, string name_, bool refresh)
		void sample_dXdPS(double * arg, vector[ vector[double] ] & FS)

	cdef cppclass f_3to2:
		f_3to2(double (*dXdPS_)(double *, size_t, void *), double (*approx_X_)(double * arg, double), double M1_, string name_, bool refresh)
		void sample_dXdPS(double * arg, vector[ vector[double] ] & FS)

cdef extern from "../src/rates.h":
	cdef cppclass rates_2to2:
		rates_2to2(Xsection_2to2 * Xprocess_, int degeneracy_, double eta_2_, string name_, bool refresh)
		double interpR(double * arg)
		void sample_initial(double * arg, vector[ vector[double] ] & IS)

	cdef cppclass rates_2to3:
		rates_2to3(Xsection_2to3 * Xprocess_, int degeneracy_, double eta_2_, string name_, bool refresh)
		double interpR(double * arg)
		void sample_initial(double * arg, vector[ vector[double] ] & IS)

	cdef cppclass rates_3to2:
		rates_3to2(f_3to2 * Xprocess_, int degeneracy_, double eta_2_, double eta_k_, string name_, bool refresh)
		double interpR(double * arg)
		void sample_initial(double * arg, vector[ vector[double] ] & IS)


#-------------Heavy quark evolution class------------------------
cdef double dtmin = 0.15
cdef class HqEvo(object):
	cdef bool elastic, inelastic, detailed_balance
	cdef Xsection_2to2 * x_Qq_Qq
	cdef Xsection_2to2 * x_Qg_Qg
	cdef Xsection_2to3 * x_Qq_Qqg
	cdef Xsection_2to3 * x_Qg_Qgg
	cdef f_3to2 * x_Qqg_Qq
	cdef f_3to2 * x_Qgg_Qg
	cdef rates_2to2 * r_Qq_Qq
	cdef rates_2to2 * r_Qg_Qg
	cdef rates_2to3 * r_Qq_Qqg
	cdef rates_2to3 * r_Qg_Qgg
	cdef rates_3to2 * r_Qqg_Qq
	cdef rates_3to2 * r_Qgg_Qg
	cdef size_t Nchannels
	cdef public vector[ vector[double] ] IS, FS
	def __cinit__(self, mass=1.3, elastic=True, inelastic=False, detailed_balance=False, table_folder='./tables', refresh_table=False):
		self.elastic=elastic
		self.inelastic=inelastic
		self.detailed_balance=detailed_balance
		self.Nchannels = 0

		if not os.path.exists(table_folder):
			os.makedirs(table_folder)


		if self.elastic:
			self.x_Qq_Qq = new Xsection_2to2(&dX_Qq2Qq_dPS, &approx_XQq2Qq, mass, "%s/XQq2Qq.hdf5"%table_folder, refresh_table)
			self.x_Qg_Qg = new Xsection_2to2(&dX_Qg2Qg_dPS, &approx_XQg2Qg, mass, "%s/XQg2Qg.hdf5"%table_folder, refresh_table)
			self.r_Qq_Qq = new rates_2to2(self.x_Qq_Qq, 36, 0., "%s/RQq2Qq.hdf5"%table_folder, refresh_table)
			self.r_Qg_Qg = new rates_2to2(self.x_Qg_Qg, 16, 0., "%s/RQg2Qg.hdf5"%table_folder, refresh_table)
			self.Nchannels += 2
			
		if self.inelastic:
			self.x_Qq_Qqg = new Xsection_2to3(&M2_Qq2Qqg, &approx_XQq2Qqg, mass, "%s/XQq2Qqg.hdf5"%table_folder, refresh_table)
			self.x_Qg_Qgg = new Xsection_2to3(&M2_Qg2Qgg, &approx_XQg2Qgg, mass, "%s/XQg2Qgg.hdf5"%table_folder, refresh_table)
			self.r_Qq_Qqg = new rates_2to3(self.x_Qq_Qqg, 36, 0., "%s/RQq2Qqg.hdf5"%table_folder, refresh_table)
			self.r_Qg_Qgg = new rates_2to3(self.x_Qg_Qgg, 16/2, 0., "%s/RQg2Qgg.hdf5"%table_folder, refresh_table)
			self.Nchannels += 2

		if self.detailed_balance:
			self.x_Qqg_Qq = new f_3to2(&Ker_Qqg2Qq, &approx_XQq2Qqg, mass, "%s/XQqg2Qq.hdf5"%table_folder, refresh_table)
			self.x_Qgg_Qg = new f_3to2(&Ker_Qgg2Qg, &approx_XQg2Qgg, mass, "%s/XQgg2Qg.hdf5"%table_folder, refresh_table)
			self.r_Qqg_Qq = new rates_3to2(self.x_Qqg_Qq, 36*16, 0., 0., "%s/RQqg2Qq.hdf5"%table_folder, refresh_table)
			self.r_Qgg_Qg = new rates_3to2(self.x_Qgg_Qg, 16*16/2, 0., 0., "%s/RQgg2Qg.hdf5"%table_folder, refresh_table)
			self.Nchannels += 2
			
		print "Number of Channels", self.Nchannels

	cpdef sample_channel(self, double E1, double T, double Tc, double dt23, double dt32):	
		cdef double r, psum = 0.0, dt, Pmax = 0.1, R1, R2, Relastic
		cdef size_t i=0
		cdef int channel_index = -1
		cdef double p[6]
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = T;
		if self.elastic:
			arg[2] = 0.
			psum += self.r_Qq_Qq.interpR(arg)
			p[i] = psum; i += 1
			psum += self.r_Qg_Qg.interpR(arg)
			p[i] = psum; i += 1
			Relastic = psum
		if self.inelastic:
			arg[2] = max(dt23, 0.1)
			R1 = self.r_Qq_Qqg.interpR(arg)
			R2 = self.r_Qg_Qgg.interpR(arg)
			psum += R1
			p[i] = psum; i += 1
			psum += R2
			p[i] = psum; i += 1
		if self.detailed_balance:
			arg[2] = max(dt32, 0.1)
			R1 = self.r_Qqg_Qq.interpR(arg)
			R2 = self.r_Qgg_Qg.interpR(arg)
			psum += R1
			p[i] = psum; i += 1
			psum += R2
			p[i] = psum; i += 1
		free(arg)
		dt = Pmax/psum
		r = (<double>rand())/dt/RAND_MAX
		if r >= p[self.Nchannels-1] or T < Tc:
			return -1, dt
		for i in range(self.Nchannels):
			if r < p[i]:
				channel_index = i
				break
		return channel_index, dt

	cpdef sample_initial(self, int channel, double E1, double T, double dt23, double dt32):
		cdef double * arg = <double*>malloc(3*sizeof(double))
		arg[0] = E1; arg[1] = T; arg[2] = 0.;
		if channel < 0:
			raise ValueError("channel must be > 0, channel < 0 correspond to no scattering")
		elif channel == 0:
			self.r_Qq_Qq.sample_initial(arg, self.IS)
		elif channel == 1:
			self.r_Qg_Qg.sample_initial(arg, self.IS)
		elif channel == 2:
			arg[2] = max(dt23, dtmin)
			self.r_Qq_Qqg.sample_initial(arg, self.IS)
		elif channel == 3:
			arg[2] = max(dt23, dtmin)
			self.r_Qg_Qgg.sample_initial(arg, self.IS)
		elif channel == 4:
			arg[2] = max(dt32, dtmin)
			self.r_Qqg_Qq.sample_initial(arg, self.IS)
		elif channel == 5:
			arg[2] = max(dt32, dtmin)
			self.r_Qgg_Qg.sample_initial(arg, self.IS)
		else:
			raise ValueError("channel %d not implememented"%channel)
		free(arg)

	cpdef sample_final(self, int channel, double s, double T, double dt23, double a1, double a2):
		cdef double * arg = <double*>malloc(4*sizeof(double))
		arg[0] = s; arg[1] = T; arg[2] = 0.; arg[3] = 0.
		if channel < 0:
			raise ValueError("channel must be > 0, channel < 0 correspond to no scattering")
		elif channel == 0:
			self.x_Qq_Qq.sample_dXdPS(arg, self.FS)
		elif channel == 1:
			self.x_Qg_Qg.sample_dXdPS(arg, self.FS)
		elif channel == 2:
			arg[2] = max(dt23, dtmin)
			self.x_Qq_Qqg.sample_dXdPS(arg, self.FS)
		elif channel == 3:
			arg[2] = max(dt23, dtmin)
			self.x_Qg_Qgg.sample_dXdPS(arg, self.FS)
		elif channel == 4:
			arg[2] = a1; arg[3] = a2
			self.x_Qqg_Qq.sample_dXdPS(arg, self.FS)
		elif channel == 5:
			arg[2] = a1; arg[3] = a2
			self.x_Qgg_Qg.sample_dXdPS(arg, self.FS)
		else:
			raise ValueError("channel %d not implememented"%channel)
		free(arg)


		
			

							
		










