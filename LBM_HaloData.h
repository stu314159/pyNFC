#ifndef LBM_HALODATA_H_
#define LBM_HALODATA_H_

class LBM_HaloData
{
public:

LBM_HaloData(const int numData, const int numSpd);
~LBM_HaloData();
void set_numData(const int nd){numData = nd;};
void set_local_nn(int* nn){local_nn = nn;};
void set_spd(int * s){spd = s;};
void set_data(float * db){data_buf = db;};
void extractHaloData(const float * f);
void distributeHaloData(float * f);

private:

// these will hold the pointers to the python HDO object
int * local_nn; // pointer to array of local node numbers
int * spd; // pointer to array of speeds
float * data_buf; // pointer to array of data

int numData; // number of data items.
const int numSpd; // number of speeds for the associated LBM lattice (needed for dist)
};


#endif
