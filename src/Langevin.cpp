#include "Langevin.h"
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include "qhat.h"

//giving vec(xi) in (px, py, pz) coordinate, Rotate vec(xi) to (px, py, py) --> (0, 0, p) coordinate
std::vector<double> rotate_toZ(std::vector<double> const& xi, double px, double py, double pz)
{
        double p_perp = std::sqrt(px*px + py*py);
        double p_length = std::sqrt(p_perp*p_perp + pz*pz);
        double cos_theta = pz / p_length;
        double sin_theta = p_perp / p_length;
        double cos_phi = px / p_perp;
        double sin_phi = py / p_perp;

        std::vector<double> Rxi;
        Rxi.resize(3);

        Rxi[0] = cos_theta*cos_phi*xi[0] + cos_theta*sin_phi*xi[1] - sin_theta*xi[2];
        Rxi[1] = - sin_phi*xi[0] + cos_phi * xi[1];
        Rxi[2] = sin_theta*cos_phi*xi[0] + sin_theta*sin_phi*xi[1] + cos_theta*xi[2];
        return Rxi;
}



//giving vec(xi) in (0, 0, p) coordinate, Rotate vec(xi) back to (px, py, pz) coordinate: (0, 0, p) --> (px, py, pz)
std::vector<double> rotate_fromZ(std::vector<double> const& xi, double px, double py, double pz)
{
        double p_perp = std::sqrt(px*px + py*py);

        if (p_perp < 1e-9) // in z direction already, no need for rotate
        {
                std::vector<double> Rxi(3);
                Rxi[0] = xi[0];
                Rxi[1] = xi[1];
                Rxi[2] = xi[2];
                return Rxi;
        }
        else // rotate
        {
                double p_length = std::sqrt(p_perp*p_perp + pz*pz);
                double cos_theta = pz / p_length;
                double sin_theta = p_perp / p_length;
                double cos_phi = px / p_perp;
                double sin_phi = py / p_perp;

                std::vector<double> Rxi(3);
                Rxi[0] = cos_theta*cos_phi*xi[0] - sin_phi*xi[1] + sin_theta*cos_phi*xi[2];
                Rxi[1] = cos_theta*sin_phi*xi[0] + cos_phi*xi[1] + sin_theta*sin_phi*xi[2];
                Rxi[2] = - sin_theta*xi[0] + cos_theta*xi[2];

                //std::cout << "Rxi: " << Rxi[0] << " " << Rxi[1] << " " << Rxi[2] << std::endl;
                return Rxi;
        }
}

// two step Langevin: step 1, giving (temp, drag, kpara, kperp), return the first step (xi_Tx, xi_Ty, xi_L)
// E is the HQ energy, M is HQ mass, temp is the local temperature
// drag, kperp, kpara have the unit of [GeV^3]
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<double> distribution(0.0, 1.0);
        
int Langevin_pre(double E, double M, double temp, double drag, double kperp, double kpara, double deltat_lrf, std::vector<double> &pre_result)
{
        double p_length = std::sqrt(E*E - M*M);

         // sample normal distribution with mean=0.0, width=1.0

        // save the values of rho as those will carry on to the post_point scheme second step
        std::vector<double> rho_xi(3, 0.0);
        for (size_t i=0; i< 3; ++i)
                rho_xi[i] = distribution(generator);


        std::vector<double> xi(3, 0.0);
        xi[0] = std::sqrt(kperp/deltat_lrf) * rho_xi[0];
        xi[1] = std::sqrt(kperp/deltat_lrf) * rho_xi[1];
        xi[2] = std::sqrt(kpara/deltat_lrf) * rho_xi[2];

        
        std::vector<double> new_p(4, 0.0);
        new_p[1] = xi[0] * deltat_lrf;
        new_p[2] = xi[1] * deltat_lrf;
        new_p[3] = p_length + (-drag * p_length + xi[2]) *deltat_lrf;
        new_p[0] = std::sqrt(M*M + new_p[0]*new_p[0] + new_p[1]*new_p[1] + new_p[2]*new_p[2]);

        pre_result={new_p[0], new_p[1], new_p[2], new_p[3], rho_xi[0], rho_xi[1], rho_xi[2]};

        // what returns is the new_p in pZ frame, and the recorded gaussian noise
        // if there is no second step, it will be pre-point scheme Langevin
	
//	std::cout << "Langevin_pre: " << pre_result[0] << " " << pre_result[1] << " " << pre_result[2] << " " << pre_result[3] <<" " << pre_result[4] << " " << pre_result[5] << " " << pre_result[6] << std::endl;
	return 0;

}


// for the two step post-point scheme, between step 1 (pre) and step 2(post), there are a few steps missing 
// the pre-step1 returns new_p, with new_p, we can update kpara and kperp 

int Langevin_post(double E, double M, double temp, double drag, double kperp, double kpara, double deltat_lrf, std::vector<double> const& pre_result, std::vector<double> & post_result)
{
        // use the pre-point stored rho
//        std::cout << "Langevin_post: " << pre_result[0] << " " << pre_result[1] << " " << pre_result[2] << " " << pre_result[3] <<" " << pre_result[4] << " " << pre_result[5] << " " << pre_result[6] << std::endl;


        std::vector<double> rho={pre_result[4], pre_result[5], pre_result[6]};

        double p_length = std::sqrt(E*E - M*M);
        
        std::vector<double> xi(3, 0.0);
        xi[0] = std::sqrt(kperp/deltat_lrf) * rho[0];
        xi[1] = std::sqrt(kperp/deltat_lrf) * rho[1];
        xi[2] = std::sqrt(kpara/deltat_lrf) * rho[2];

    	post_result.resize(4);
        post_result[1] = xi[0] * deltat_lrf;
        post_result[2] = xi[1] * deltat_lrf;
        post_result[3] = p_length + (-drag * p_length + xi[2]) * deltat_lrf;
        post_result[0] = std::sqrt(M*M + post_result[1]*post_result[1] + post_result[2]*post_result[2] + post_result[3]*post_result[3]);

//	std::cout << "Langevin_post: " << post_result[0] << " " << post_result[1] << " " << post_result[2] << " " << post_result[3] << std::endl;
        return 0;
}


int update_by_Langevin(particle &HQ, Qhat_2to2* qhatQq2Qq, Qhat_2to2* qhatQg2Qg, double temp, double deltat_lrf, bool EinR=false)
{
        double GeV_to_Invfm = 5.068;
        double drag_Qq, drag_Qg, drag;
        double kperp_Qq, kperp_Qg, kperp;
        double kpara_Qq, kpara_Qg, kpara;

        double p_perp = std::sqrt(HQ.p[1]*HQ.p[1] + HQ.p[2]*HQ.p[2]);
        double p_length = std::sqrt(p_perp*p_perp + HQ.p[3]*HQ.p[3]);
        double M = std::sqrt(HQ.p[0] * HQ.p[0] - p_length*p_length);

        // step 1: pre-point shceme, estimate delta_p
        double* args = new double[4];
        args[0] = HQ.p[0];
        args[1] = temp;
        args[2] = 0.0;
        args[3] = 0; drag_Qq = qhatQq2Qq->interpQ(args); drag_Qg = qhatQg2Qg->interpQ(args);
        args[3] = 1; kperp_Qq = qhatQq2Qq->interpQ(args); kperp_Qg = qhatQg2Qg->interpQ(args);
        args[3] = 2; kpara_Qq = qhatQq2Qq->interpQ(args); kpara_Qg = qhatQg2Qg->interpQ(args);
        delete [] args;

        drag = (drag_Qq + drag_Qg) / p_length * GeV_to_Invfm;
        kperp = (kperp_Qq + kperp_Qg) * GeV_to_Invfm;
        kpara = (kpara_Qq  + kpara_Qg) * GeV_to_Invfm;
        if (EinR) drag = kperp/(2*temp*HQ.p[0]);

	std::vector<double> pre_result;
	std::vector<double> post_result;
        Langevin_pre(HQ.p[0], M, temp, drag, kperp, kpara, deltat_lrf, pre_result);
        std::vector<double> new_p = {pre_result[0], pre_result[1], pre_result[2], pre_result[3]};



        double new_Energy = new_p[0];

        double* args_ = new double[4];
        args_[0] = new_Energy;
        args_[1] = temp;
        args_[2] = 0.0;
        args_[3] = 1; kperp_Qq = qhatQq2Qq->interpQ(args); kperp_Qg = qhatQg2Qg->interpQ(args);
        args_[3] = 2; kpara_Qq = qhatQq2Qq->interpQ(args); kpara_Qg = qhatQg2Qg->interpQ(args);
        delete [] args_;

        kperp = (kperp_Qq + kperp_Qg) * GeV_to_Invfm;
        kpara = (kpara_Qq + kpara_Qg) * GeV_to_Invfm;
        Langevin_post(HQ.p[0], M, temp, drag, kperp, kpara, deltat_lrf, pre_result, post_result);

        std::vector<double> new3_p = {post_result[1], post_result[2], post_result[3]};
       
        // rotate back from Z direction to lrf frame

        std::vector<double> p_lrf(3, 0.0);
        p_lrf = rotate_fromZ(new3_p, HQ.p[1], HQ.p[2], HQ.p[3]);

        for (size_t i=0; i<3; ++i)
        {
                HQ.x[i] += HQ.p[i+1]/HQ.p[0] * deltat_lrf;
                HQ.p[i+1] = p_lrf[i];
        }

        HQ.p[0] = std::sqrt(M*M + HQ.p[1]*HQ.p[1] + HQ.p[2]*HQ.p[2] + HQ.p[3]*HQ.p[3]);
        return 0;


}
