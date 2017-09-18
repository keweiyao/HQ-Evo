#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <boost/multi_array.hpp>
#include <H5Cpp.h>

#include "matrix_elements.h"
#include "Xsection.h"
#include "rates.h"
//#include "qhat_Xsection.h"
#include "qhat.h"
#include "Langevin.h"


using std::vector;

int main()
{
    RadiationHT rad = RadiationHT(4, "dNg_over_dxdydt.dat", true, false);

    double time=2.0;
    double temp = 0.18;
    double HQenergy = 50.;
    double deltat_lrf = 0.2;
    double qhat = M_PI * pow(temp , 3);

    double result, max_dNg;
    rad.calculate(time, temp, HQenergy, result, max_dNg);
    std::cout << "gluon emission probability: " << result * qhat * deltat_lrf << " " << max_dNg << std::endl; 
    bool emit, success;
    int count1 = 0, count2 = 0, trial = 1e4;
    std::vector<double> gluon;
    for (int i=0; i<trial; ++i)
    {
        emit = rad.emitGluon(time, temp, HQenergy, qhat, deltat_lrf);
        if (emit)
        {
            count1 ++;
            success = rad.sampleGluon(time, temp, HQenergy, qhat, deltat_lrf, gluon);
            if (success) count2 ++;
        }
    }

    std::cout << "in reality: " << 1.*count1/trial << " " << 1.*count2/trial << std::endl;

/*
        double M = 1.3;

        double temp = 0.3;

        initialize_Debye_mass(0, 2.2, 1.0, 0.0, 0.154);
        
        bool refresh = true;
        Xsection_2to2 xQq2Qq(&dX_Qq2Qq_dPS, M, "XQq2Qq.hdf5", refresh);
        Xsection_2to2 xQg2Qg(&dX_Qg2Qg_dPS, M, "XQg2Qg.hdf5", refresh);
        rates_2to2 rQq2Qq(&xQq2Qq, 36, 0., "rQq2Qq.hdf5", refresh);
        rates_2to2 rQg2Qg(&xQg2Qg, 16, 0., "rQg2Qg.hdf5", refresh);

        QhatXsection_2to2 qhat_xQq2Qq(&dqhat_Qq2Qq_dPS, M, "qhat_XQq2Qq.hdf5", refresh);
        QhatXsection_2to2 qhat_xQg2Qg(&dqhat_Qg2Qg_dPS, M, "qhat_XQg2Qg.hdf5", refresh);
        Qhat_2to2 qhatQq2Qq(&qhat_xQq2Qq, 36, 0., "qhat_Qq2Qq.hdf5", refresh);
        Qhat_2to2 qhatQg2Qg(&qhat_xQg2Qg, 16, 0., "qhat_Qg2Qg.hdf5", refresh);
*/


/*

        std::cout << "read in table successfully :) " << std::endl;

        particle testHQ;
    
        double E1 = 10.;
    
        testHQ.p={E1, 0., 0, std::sqrt(E1*E1 - M*M)};
        testHQ.x={0., 0., 0.};

        double drag_Qq, drag_Qg, drag;
        double kperp_Qq, kperp_Qg, kperp;
        double kpara_Qq, kpara_Qg, kpara;

        int npart = 20000;
        vector<particle> events(npart, testHQ);

        int ntime = 10000;
        double deltat = 0.01;

        int NTau = int(ntime*deltat);


        typedef boost::multi_array<double, 3> array3D;
        array3D part_info(boost::extents[NTau][npart][4]);

*/
        // test 
        /*
        particle HQ;
        HQ.p = testHQ.p;
        HQ.x = testHQ.x;

        double* args = new double[4];
        args[0] = HQ.p[0];
        args[1] = temp;
        args[2] = 0.0;
 
        args[3] = 1;  drag_Qq = qhatQq2Qq.interpQ(args); drag_Qg = qhatQg2Qg.interpQ(args);
        args[3] = 2;  kperp_Qq = qhatQq2Qq.interpQ(args); kperp_Qg = qhatQg2Qg.interpQ(args);
        args[3] = 3;  kpara_Qq = qhatQq2Qq.interpQ(args); kpara_Qg = qhatQg2Qg.interpQ(args);
        delete [] args;
        drag = (drag_Qq + drag_Qg) * 5.068;
        kperp = (kperp_Qq + kperp_Qg) * 5.068;
        kpara = ((kpara_Qq - drag_Qq*drag_Qq) + (kpara_Qg - drag_Qg*drag_Qg) ) * 5.068;

        std::cout << "main: " << HQ.p[0] <<" " <<  drag << " " << kperp << " " << kpara << std::endl;


        //update_by_Langevin_test(HQ, &qhatQq2Qq, &qhatQg2Qg, temp, deltat, true);
        update_by_Langevin(HQ, temp, drag, kpara, kperp, deltat, false);
        std::cout << "main: " << HQ.p[0] <<" " <<  HQ.p[1]<< " " << HQ.p[2] << " " << HQ.p[3] << std::endl;
        */

/*
        for (size_t ipart = 0; ipart < npart; ++ipart)
        {
                particle HQ;
                HQ.p = testHQ.p;
                HQ.x = testHQ.x;

                for (size_t itime=0; itime<ntime; ++itime)
                {
                        update_by_Langevin_test(HQ, &qhatQq2Qq, &qhatQg2Qg, temp, deltat, false);
                        if (itime % 100 == 0)
                        {
                                int idx_t = itime/100;
                                part_info[idx_t][ipart][0] = HQ.p[0];
                                part_info[idx_t][ipart][1] = HQ.p[1];
                                part_info[idx_t][ipart][2] = HQ.p[2];
                                part_info[idx_t][ipart][3] = HQ.p[3];
                        }
                }

        }
      
*/
 
/*
        for (int ipart = 0; ipart < npart; ++ipart)
        {

        particle HQ;
        HQ.p = testHQ.p;
        HQ.x = testHQ.x;

        for (int itime = 0; itime < ntime; ++itime)
        {
                double *args = new double[4];
                args[0] = HQ.p[0];
                args[1] = temp;
                args[2] = 0.0;
                args[3] = 1;  drag_Qq = qhatQq2Qq.interpQ(args); drag_Qg = qhatQg2Qg.interpQ(args);
                args[3] = 2;  kperp_Qq = qhatQq2Qq.interpQ(args); kperp_Qg = qhatQg2Qg.interpQ(args);
                args[3] = 3;  kpara_Qq = qhatQq2Qq.interpQ(args); kpara_Qg = qhatQg2Qg.interpQ(args);

                drag = (drag_Qq + drag_Qg)*5.068;
                kperp = (kperp_Qq + kperp_Qg) * 5.068;
                kpara = ((kpara_Qq - drag_Qq*drag_Qq) + (kpara_Qg - drag_Qg*drag_Qg)) * 5.068;

                delete [] args;

                update_by_Langevin(HQ, temp, drag, kpara, kperp, deltat, false);

                if (itime % 100 == 0)
                {
                    int idx_t = itime/100;
                    part_info[idx_t][ipart][0] = HQ.p[0];
                    part_info[idx_t][ipart][1] = HQ.p[1];
                    part_info[idx_t][ipart][2] = HQ.p[2];
                    part_info[idx_t][ipart][3] = HQ.p[3];
                }
        }
        }
*/


/*

        H5::H5File* file = new H5::H5File("Event_static_EinRFalse_preIto.hdf5", H5F_ACC_TRUNC);
        std::string datasetname("Langevin-events");
        const size_t rank=3;
        hsize_t dims[rank] = {NTau, npart, 4};
        H5::DSetCreatPropList proplist{};
        proplist.setChunk(rank, dims);

        H5::DataSpace dataspace(rank, dims);
        auto datatype(H5::PredType::NATIVE_DOUBLE);
        H5::DataSet* dataset = new H5::DataSet( file->createDataSet(datasetname, datatype, dataspace, proplist) );
        dataset->write(part_info.data(), datatype);

        delete dataset;
        delete file;
*/

        return 0;
}
