from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdlib cimport malloc, free
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport *
import os

#------------------Import C++ fucntions and class for Xsection and rates------------------
cdef extern from "../src/matrix_elements.h":
        cdef double approx_XQq2Qq(double* args, double M)
        cdef double approx_XQg2Qg(double* args, double M)

        cdef double dqhat_Qq2Qq_dPS(double* PS, size_t ndims, void* params)
        cdef double dqhat_Qg2Qg_dPS(double* PS, size_t ndims, void* params)



cdef extern from "../src/qhat_Xsection.h":
        cdef cppclass QhatXsection_2to2:
                QhatXsection_2to2(double (*dXdPS_)(double*, size_t, void *), double (*approx_X_)(double*, double), double M1_, string name_, bool refresh)
                double calculate(double *args)
                double interpX(double *args)

cdef extern from "../src/qhat.h":
        cdef cppclass Qhat_2to2:
                Qhat_2to2(QhatXsection_2to2 *Xprocess_, int degeneracy_, double eta_2_, string name_, bool refresh)
                double calculate(double *args)
                double interpQ(double *args)
                


# ----- Wrapper for QhatXsection_2to2 and Qhat class ------------
cdef class pyQhatX2to2:
        cdef QhatXsection_2to2 *cQhatX2to2
        def __cinit__(self, channel, double mass, string filename, bool refresh):
                if channel == 'Qq->Qq':
                        self.cQhatX2to2 = new QhatXsection_2to2(&dqhat_Qq2Qq_dPS, &approx_XQq2Qq, mass, filename, refresh)
                elif channel == 'Qg->Qg':
                        self.cQhatX2to2 = new QhatXsection_2to2(&dqhat_Qg2Qg_dPS, &approx_XQg2Qg, mass, filename, refresh)
                else:
                        raise ValueError("channel %s is not implemented"%channel)

        cpdef double calculate(self, double s, double Temp, int index):
                cdef double * args = <double*>malloc(3*sizeof(double))
                args[0] = s; args[1] = Temp; args[2] = index
                cdef double result = self.cQhatX2to2.calculate(args)
                free(args)
                return result

        cpdef double interpX(self, double s, double Temp, int index):
                cdef double * args = <double*>malloc(3*sizeof(double))
                args[0] = s; args[1] = Temp; args[2] = index
                cdef double result = self.cQhatX2to2.interpX(args)
                free(args)
                return result

cdef class pyQhat2to2:
        cdef Qhat_2to2 *cQhat2to2
        def __cinit__(self, pyQhatX2to2 x2to2, int degeneracy, double eta2, string filename, bool refresh):
                self.cQhat2to2 = new Qhat_2to2(x2to2.cQhatX2to2, degeneracy, eta2, filename, refresh)

        cpdef double calculate(self, double &E1, double &Temp, int &iweight, int &qidx):
                cdef double * args = <double*>malloc(4*sizeof(double))
                args[0] = E1; args[1] = Temp; args[2] = iweight; args[3] = qidx
                cdef double result = self.cQhat2to2.calculate(args)
                free(args)
                return result

        cpdef double interpQ(self, double &E1, double &Temp, int &iweight, int &qidx):
                cdef double * args = <double*>malloc(4*sizeof(double))
                args[0] = E1; args[1] = Temp; args[2] = iweight; args[3] = qidx
                cdef double result = self.cQhat2to2.interpQ(args)
                free(args)
                return result

#------------- Import c++ function for Langevin evolution
cdef extern from "../src/Langevin.h":
        cdef vector[double] Langevin_pre(double E1, double mass, double temp, double drag, double kpara, double kperp, double delta_lrf)
        cdef vector[double] Langevin_post(double E1, double mass, double temp, double drag, double kpara, double kperp, double delta_lrf, vector[double] & pre_result)



#------------ Heavy quark Langevin transport evolution class -------------
cdef class HqLGV:
        cdef bool elastic, EinR
        cdef object qhat_Qq2Qq, qhat_Qg2Qg  # no need to qhatX_Qq2Qq object, since we are not actually publicly using them....
        cdef size_t Nchannels, N22
        cdef double mass, deltat_lrf
        def __cinit__(self, mass =1.3, elastic = True, EinR = False, deltat_lrf=0.01, table_folder='./tables', refresh_table=False):
                self.elastic = elastic
                self.EinR = EinR
                self.mass = mass
                self.deltat_lrf = deltat_lrf
                if not os.path.exists(table_folder):
                        os.makedirs(table_folder)

                if self.elastic:
                        qhatX_Qq2Qq = pyQhatX2to2('Qq->Qq', mass, "%s/QhatX_Qq2Qq.hdf5"%table_folder, refresh_table)
                        qhatX_Qg2Qg = pyQhatX2to2('Qg->Qg', mass, "%s/QhatX_Qg2Qg.hdf5"%table_folder, refresh_table)
                        self.qhat_Qq2Qq = pyQhat2to2(qhatX_Qq2Qq, 36., 0., "%s/Qhat_Qq2Qq.hdf5"%table_folder, refresh_table)
                        self.qhat_Qg2Qg = pyQhat2to2(qhatX_Qg2Qg, 16., 0., "%s/Qhat_Qg2Qg.hdf5"%table_folder, refresh_table)
        
        # giving the incoming heavy quark energy E1 in cell frame, return new_p(p0, px, py, pz) in (0, 0, p_length) frame
        cpdef update_by_Langevin(self, double E1, double temp):
                cdef double GeV_to_Invfm = 5.068
                cdef double drag_Qq, drag_Qg, drag, kperp_Qq, kperp_Qg, kperp, kpara_Qq, kpara_Qg
                cdef double p_length = sqrt(E1*E1 - self.mass*self.mass)
                drag_Qq = self.qhat_Qq2Qq.interpQ(E1, temp, 0, 0)
                drag_Qg = self.qhat_Qg2Qg.interpQ(E1, temp, 0, 0)
                kperp_Qq = self.qhat_Qq2Qq.interpQ(E1, temp, 0, 1)
                kperp_Qg = self.qhat_Qg2Qg.interpQ(E1, temp, 0, 1)
                kpara_Qq = self.qhat_Qq2Qq.interpQ(E1, temp, 0, 2)
                kpara_Qg = self.qhat_Qg2Qg.interpQ(E1, temp, 0, 2)

                drag = (drag_Qq + drag_Qg) / p_length * GeV_to_Invfm
                kperp = (kperp_Qq + kperp_Qg) * GeV_to_Invfm
                kpara = (kpara_Qq + kpara_Qg) * GeV_to_Invfm

                if self.EinR:
                        drag = kperp / (2.*temp*E1);

                cdef vector[double] pre_result = Langevin_pre(E1, self.mass, temp, drag, kperp, kpara, self.deltat_lrf)
                cdef double new_energy = sqrt(self.mass**2 + pre_result[0]**2 + pre_result[1]**2 + pre_result[2]**2)

                kperp_Qq = self.qhat_Qq2Qq.interpQ(new_energy, temp, 0, 1)
                kperp_Qg = self.qhat_Qg2Qg.interpQ(new_energy, temp, 0, 1)
                kpara_Qq = self.qhat_Qq2Qq.interpQ(new_energy, temp, 0, 2)
                kpara_Qg = self.qhat_Qg2Qg.interpQ(new_energy, temp, 0, 2)

                kperp = (kperp_Qq + kperp_Qg) * GeV_to_Invfm
                kpara = (kpara_Qq + kpara_Qg) * GeV_to_Invfm

                cdef vector[double] new_p = Langevin_post(E1, self.mass, temp, drag, kperp, kpara, self.deltat_lrf, pre_result)
                return new_p

                
               
