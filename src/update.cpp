#include "update.h"
#include <iostream>


update::update(double M)
:	rd(),
	gen(rd()),
	dist01(0., 1.),
	XQq2Qq(&dX_Qq2Qq_dPS, &approx_XQq2Qq, M, "./tables/X-Qq-Qq.dat"), // channel -- 1
	XQg2Qg(&dX_Qg2Qg_dPS, &approx_XQg2Qg, M, "./tables/X-Qg-Qg.dat"), // channel -- 2
	XQq2Qqg(&M2_Qq2Qqg, &approx_XQq2Qqg, M, "./tables/X-Qq-Qqg.dat"), // channel -- 3
	XQg2Qgg(&M2_Qg2Qgg, &approx_XQg2Qgg, M, "./tables/X-Qg-Qgg.dat"), // channel -- 4
	RQq2Qq(&XQq2Qq, 3*4, "./tables/R-Qq-Qq.dat"),
	RQg2Qg(&XQg2Qg, 8*2, "./tables/R-Qg-Qg.dat"),
	RQq2Qqg(&XQq2Qqg, 3*4, "./tables/R-Qq-Qqg.dat"),
	RQg2Qgg(&XQg2Qgg, 8*2, "./tables/R-Qg-Qgg.dat")
{
	std::cout << "hello" << std::endl;
}

int update::sample_channel(double E1, double T, double dt){

}
