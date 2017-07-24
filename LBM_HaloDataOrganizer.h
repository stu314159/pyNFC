/*
 * LBM_HaloDataOrganizer.h
 *
 *  Created on: Jul 21, 2017
 *      Author: stu
 */

#ifndef LBM_HALODATAORGANIZER_H_
#define LBM_HALODATAORGANIZER_H_

#include "LBM_HaloData.h"
#include <map>

class LBM_HaloDataOrganizer
{
public:
	LBM_HaloDataOrganizer();
	~LBM_HaloDataOrganizer();
	void insert(const int ngbNum,const int numData, const int numSpd);


private:
	std::map<int,LBM_HaloData> HaloData;

};




#endif /* LBM_HALODATAORGANIZER_H_ */
