/*
 * LBM_HaloDataOrganizer.h
 *
 *  Created on: Jul 21, 2017
 *      Author: stu
 */

#ifndef LBM_HALODATAORGANIZER_H_
#define LBM_HALODATAORGANIZER_H_

#include "LBM_DataHandler.h"
#include <map>

class LBM_HaloDataOrganizer
{
public:
LBM_HaloDataOrganizer();
~LBM_HaloDataOrganizer();


private:
std::map<int,LBM_DataHandler> HaloData;

};




#endif /* LBM_HALODATAORGANIZER_H_ */
