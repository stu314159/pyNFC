#include "LBM_HaloData.h"
#include <cstdlib>

LBM_HaloData::LBM_HaloData(const int numData, const int numSpd) :
local_nn(NULL), spd(NULL), data_buf(NULL), numData(numData),numSpd(numSpd)
{


}

LBM_HaloData::~LBM_HaloData()
{


}

void LBM_HaloData::extractHaloData(const float * f)
{
	for(int hIdx=0;hIdx<numData;hIdx++)
	{
		data_buf[hIdx] = f[local_nn[hIdx]*numSpd+spd[hIdx]];
	}


}

void LBM_HaloData::distributeHaloData(float * f)
{
	for(int hIdx=0;hIdx<numData;hIdx++)
	{
		f[local_nn[hIdx]*numSpd+spd[hIdx]] = data_buf[hIdx];
	}

}
