/*
 * LBM_HaloDataOrganizer.cpp
 *
 *  Created on: Jul 21, 2017
 *      Author: stu
 */
#include "LBM_HaloDataOrganizer.h"


LBM_HaloDataOrganizer::LBM_HaloDataOrganizer()
{

}

LBM_HaloDataOrganizer::~LBM_HaloDataOrganizer()
{

}


void LBM_HaloDataOrganizer::insert_ngb(const int ngbNum,const int numData, const int numSpd)
{
	HaloData[ngbNum]=LBM_HaloData(numData,numSpd);
}

void LBM_HaloDataOrganizer::initialize_ngb_pointers(const int ngbNum, int * nd_num, int * spd, float * data)
{
	HaloData[ngbNum].set_local_nn(nd_num);
	HaloData[ngbNum].set_spd(spd);
	HaloData[ngbNum].set_data(data);
}
