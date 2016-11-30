from libcpp.vector cimport vector

#------------Import c++ utility function for rotation and transformation-----------
cdef extern from "../src/utility.h":
	cdef double product4(const vector[double] & A, const vector[double] & B)
	cdef void print4vec(const vector[double] & A)
	cdef void rotate_Euler(const vector[double] & A, vector[double] & Ap, double alpha, double beta, double gamma)
	cdef void rotate_axis(const vector[double] & A, vector[double] & Ap, double alpha, unsigned int dir)
	cdef void boost_by3(const vector[double] & A, vector[double] & Ap, const vector[double] v)
	cdef void boost_by4(const vector[double] & A, vector[double] & Ap, const vector[double] u)
	cdef void boost_axis(const vector[double] & A, vector[double] & Ap, const double vd, unsigned int dir)
	cdef void go_to_CoM(const vector[double] & Pcom,
			   const vector[double] & A, const vector[double] & B,
			   vector[double] & Ap, vector[double] & Bp);

#-------------C/Python Wrapper functions--------------------------------------------
cpdef double dot4(vector[double] & A, vector[double] & B):
	return product4(A, B)

cpdef print4(const vector[double] & A):
	print4vec(A)

cpdef rotate_ByEuler(vector[double] & A, double alpha, double beta, double gamma):
	cdef vector[double] Ap
	Ap.resize(4)
	rotate_Euler(A, Ap, alpha, beta, gamma)
	return Ap

cpdef rotate_ByAxis(vector[double] & A, double alpha, unsigned int dir):
	if not (dir in [1,2,3]):
		raise ValueError("Direction can only be 1(x), 2(y), 3(z)")
	cdef vector[double] Ap
	Ap.resize(4)
	rotate_axis(A, Ap, alpha, dir)
	return Ap

cpdef boost4_By3(vector[double] A, vector[double] v):
	cdef vector[double] Ap
	Ap.resize(4)
	boost_by3(A, Ap, v)
	return Ap

cpdef boost4_By4(vector[double] A, vector[double] u):
	cdef vector[double] Ap
	Ap.resize(4)
	boost_by4(A, Ap, u)
	return Ap

cpdef boost4_ByAxis(vector[double] A, double vd, unsigned int dir):
	if not (dir in [1,2,3]):
		raise ValueError("Direction can only be 1(x), 2(y), 3(z)")
	cdef vector[double] Ap
	Ap.resize(4)
	boost_axis(A, Ap, vd, dir)
	return Ap

cpdef boost4_ByCoM(vector[double] A, vector[double] B):
	cdef vector[double] Ap, Bp, Pcom
	Ap.resize(4)
	Bp.resize(4)
	Pcom.resize(4)
	cdef i
	for i in range(4):
		Pcom[i] = A[i] + B[i]
	go_to_CoM(Pcom, A, B, Ap, Bp)
	return Pcom, Ap, Bp


