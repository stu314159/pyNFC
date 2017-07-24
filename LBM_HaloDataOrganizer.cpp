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


void LBM_HaloDataOrganizer::insert(const int ngbNum,const int numData, const int numSpd)
{
	HaloData[ngbNum]=LBM_HaloData(numData,numSpd);
}
