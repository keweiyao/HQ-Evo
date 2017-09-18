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

#------------------Import C++ fucntions and class for Xsection and rates------------------
cdef extern from "../src/matrix_elements.h":
    cdef void initialize_Debye_mass(const unsigned int mDtype, const double mDTc,
                           const double mDslope, const double mDcurv, 
                           const double Tc)

    cdef double dqhat_Qq2Qq_dPS(double* PS, size_t ndims, void* params) 
    cdef double dqhat_Qg2Qg_dPS(double* PS, size_t ndims, void* params) 

cdef extern from "../src/qhat_Xsection.h":
    cdef cppclass QhatXsection_2to2:
        QhatXsection_2to2(double (*dXdPS_)(double*, size_t, void *), double M1_, string name_, bool refresh) 
        double calculate(double *args) 
        double interpX(double *args) 

cdef extern from "../src/qhat.h":
    cdef cppclass Qhat_2to2:
        Qhat_2to2(QhatXsection_2to2 *Xprocess_, int degeneracy_, double eta_2_, string name_, bool refresh) 
        double calculate(double *args) 
        double interpQ(double *args) 


cdef extern from "../src/Langevin.h":
    cdef cppclass RadiationHT:
        RadiationHT(double mass, string filename, bool plain, bool refresh)
        double getNg(double time, double temp, double HQenergy, double qhat)
        bool emitGluon(double time, double temp, double HQenergy, double qhat, double deltat);
        bool sampleGluon(double time, double temp, double HQenergy, double qhat, double deltat, vector[double] & gluon)



#------------- Import c++ function for Langevin evolution
cdef extern from "../src/Langevin.h":
    cdef void Langevin_pre(double p_length, double mass, double temp, double drag, double kpara, double kperp, double delta_lrf, vector[double] & pre_result) 
    cdef void Langevin_post(double p_length, double mass, double temp, double drag, double kpara, double kperp, double delta_lrf, vector[double] & pre_result, vector[double] & post_result) 



#------------ Heavy quark Langevin transport evolution class -------------
cdef class HqLGV:
    cdef bool collisional, EinR, radiative
    cdef QhatXsection_2to2 * qhatX_Qq2Qq
    cdef QhatXsection_2to2 * qhatX_Qg2Qg
    cdef Qhat_2to2 * qhat_Qq2Qq
    cdef Qhat_2to2 * qhat_Qg2Qg 
    cdef RadiationHT * radHT
    cdef size_t Nf
    cdef double mass, deltat_lrf, Kfactor
    cdef public vector[double] pre_result, post_result
        
    def __cinit__(self, options, table_folder='./tables', refresh_table=False):
        self.mass = options['mass']
        self.Nf = options['Nf']
        self.elastic = options['transport']['elastic']
        self.EinR = options['transport']['Einstein']
        self.radiative = options['transport']['radiative']
        self.deltat_lrf = options['transport']['dt_lrf']
        self.Kfactor = options['Kfactor']

        # set mD
        cdef double Tc = options['Tc']
        cdef unsigned int mD_type = options['mD']['mD-model']
        cdef double mDTc = options['mD']['mTc']
        cdef double mDslope = options['mD']['slope']
        cdef double mDcurv = options['mD']['curv']
        print "mD:", mD_type, mDTc, mDslope, mDcurv, Tc
        initialize_Debye_mass(mD_type, mDTc, mDslope, mDcurv, Tc)

        # initialize diffusion coefficient (with unit GeV*fm-2)
        self._drag = 0.
        self._kperp = 0.
        self._kpara = 0.

        if not os.path.exists(table_folder):
            os.makedirs(table_folder)

# this is the qhat calculated from pQCD calculation
        if self.elastic:
            self.qhatX_Qq2Qq = new QhatXsection_2to2(&dqhat_Qq2Qq_dPS, self.mass, "%s/QhatX_Qq2Qq.hdf5"%table_folder, refresh_table)
            self.qhatX_Qg2Qg = new QhatXsection_2to2(&dqhat_Qg2Qg_dPS, self.mass, "%s/QhatX_Qg2Qg.hdf5"%table_folder, refresh_table)
            self.qhat_Qq2Qq = new Qhat_2to2(self.qhatX_Qq2Qq, 12*self.Nf, 0., "%s/Qhat_Qq2Qq.hdf5"%table_folder, refresh_table)
            self.qhat_Qg2Qg = new Qhat_2to2(self.qhatX_Qg2Qg, 16, 0., "%s/Qhat_Qg2Qg.hdf5"%table_folder, refresh_table)
    
# this is emitted gluon from HT calculation
        if self.radiative:
            self.radHT = new RadiationHT(self.mass, "%s/dNg_dxdydt_HT.dat"%table_folder, True, refresh_table)


    cpdef getDiffusion(self, double temp, double HQenergy, int case):
        #input: medium temp, HQenergy, case (0 for Langevin-pre, 1 for Langevin-post)
        #return: assign value for self._kperp, self._kpara, self._drag(case 0), unit [GeV.fm-2]
        cdef double drag_Qq, drag_Qg
        cdef double kperp_Qq, kperp_Qg, kpara_Qq, kpara_Qg

        cdef double p_length = sqrt(HQenergy*HQenergy - self.mass*self.mass)
        cdef double *arg = <double*> malloc(3*sizeof(double))
        arg[0] = HQenergy; arg[1] = temp;
        arg[2] = 1
        kperp_Qq = self.qhat_Qq2Qq.interpQ(arg)
        kperp_Qg = self.qhat_Qg2Qg.interpQ(arg)
        arg[2] = 2
        kpara_Qq = self.qhat_Qq2Qq.interpQ(arg)
        kperp_Qg = self.qhat_Qg2Qg.interpQ(arg)

        self._kperp = (kperp_Qq + kperp_Qg) * GeV_to_Invfm
        self._kpara = (kpara_Qq + kpara_Qg) * GeV_to_invfm  # problematic in this definition!!
        
        if not case: # Langevin_pre
            if self.EinR:
                self._drag = self._kperp / (2.*temp*HQenergy)
            else:
                arg[2] = 0
                drag_Qq = self.qhat_Qq2Qq.interpQ(arg) 
                drag_Qg = self.qhat_Qg2Qg.interpQ(arg)
                self._drag = (drag_Qq + drag_Qg) / p_length * GeV_toInvfm

        free(arg)
        
        
    cpdef update_by_Langevin(self, double E1, double temp):
        #input: incoming HQenergy (E1 in cell frame)
        #return: new HQ momentum (p0, px, py, pz) in (0, 0, p_length) frame
        cdef double p_length = sqrt(E1*E1 - self.mass*self.mass)
        getDiffusion(E1, temp, 0)
        Langevin_pre(p_length, self.mass, temp, self._drag, self._kperp, self._kpara, self.deltat_lrf, self.pre_result)

        cdef double new_energy = self.pre_result[0]
        getDiffusion(new_energy, temp, 1)
        Langevin_post(p_length, self.mass, temp, self._drag, self._kperp, self._kpara, self.deltat_lrf, self.pre_result)
