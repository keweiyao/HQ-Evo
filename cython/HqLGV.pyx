# cython: c_string_type=str, c_string_encoding=ascii
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdlib cimport malloc, free
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport *
import cython
import os

cdef double GeV_to_Invfm = 5.068

#------------- Import c++ function for Langevin evolution
cdef extern from "../src/Langevin.h":
	cdef double kperp(double p, double M, double T)
	cdef double kpara(double p, double M, double T)
	cdef void Langevin_step(double pz0, double M, double T, 
							double delta_t_lrf, vector[double] & pnew)
	cdef void initialize_transport_coeff(double A, double B)



#------------ Heavy quark Langevin transport evolution class -------------
cdef class HqLGV:
	cdef double mass
	cdef public vector[double] pnew
		
	def __cinit__(self, options):
		self.mass = options['mass']
		cdef double A = options['transport']['A']
		cdef double B = options['transport']['B']
		initialize_transport_coeff(A, B)

	# Giving [E0, 0, 0, pz0] return [E', px', py', pz']
	# delta_t_lrf [GeV^-1]
	cpdef update_by_Langevin(self, double E0, double T, double delta_t_lrf):
		Langevin_step(E0, self.mass, T, delta_t_lrf, self.pnew)
		return self.pnew
