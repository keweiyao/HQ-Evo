#ifndef UTILITY_H
#define UTILITY_H
#include <cmath>
#include <vector>
#include <iostream>

double product4(std::vector<double> const& A, std::vector<double> const& B);
void print4vec(std::vector<double> const& A);
 
//=======================Rotation==============================================================
// Rotation can take care of 
//---------------------General Euler angles----------------------------------------------------
// This rotation operation takes 3 euler angles and rotate (new frame relative to the old frame)
// And returns the A vector by its components in the new frame.
// the new frame is achieved by the old frame from Z(1_=alpha), X(2_=beta) and Z(3_=gamma)
// R = Z3T*X2T*Z1T
void rotate_Euler(std::vector<double> const& A, std::vector<double> & Ap, double alpha, double beta, double gamma);

// Rotation around ith axis, return vector components in the new frame, passive
// dir = 1(x), 2(y), 3(z)
void rotate_axis(std::vector<double> const& A, std::vector<double> &Ap, double alpha, unsigned int dir);

//=======================Boost==============================================================

//---------------------General Boost (beta_x, beta_y, beta_z)-----------------------------------------
// This boost operation takes 3 velocity (vx, vy, vz, of the new frame relative to the old frame)
// And returns the A vector by its components in the new frame.
void boost_by3(std::vector<double> const& A, std::vector<double> &Ap, std::vector<double> const& v);

// This boost operation takes 4 velocity (u0, u1, u2, u3 of the new frame relative to the old frame)
// And returns the A vector by its components in the new frame.
void boost_by4(std::vector<double> const& A, std::vector<double> &Ap, std::vector<double> const& u);
// Boost along ith axis, return vector components in the new frame, passive
// dir = 1(x), 2(y), 3(z)
void boost_axis(std::vector<double> const& A, std::vector<double> &Ap, double const vd, unsigned int dir);

// boost two 4-vectors to their center of mass frame:
void go_to_CoM(std::vector<double> const& Pcom,
			   std::vector<double> const& A, std::vector<double> const& B,
			   std::vector<double> & Ap, std::vector<double> & Bp);


#endif
