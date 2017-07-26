#ifndef PYLBM_INTERFACE_H_
#define PYLBM_INTERFACE_H_

#include "Lattice.h"
#include "D3Q15Lattice.h"
#include "D3Q19Lattice.h"
#include "D3Q27Lattice.h"
#include "LBM_DataHandler.h"
#include "LBM_HaloData.h"
#include "LBM_HaloDataOrganizer.h"

#include <boost/python.hpp>
#include <cstdlib>
#include <omp.h>

class PyLBM_Interface
{
public:
	PyLBM_Interface(const int numSpd);
	~PyLBM_Interface();
	void set_fIn(const float* fIn, const int nd,LBM_DataHandler& f);
	//void get_fOut(boost::python::object obj);
	void set_fEven(boost::python::object obj);
	void set_fOdd(boost::python::object obj);
	void set_ux(boost::python::object obj);
	void set_uy(boost::python::object obj);
	void set_uz(boost::python::object obj);
	void set_rho(boost::python::object obj);
	void compute_local_data(const bool isEven);
	void set_adjacency(boost::python::object obj);
	void set_boundaryNL(boost::python::object obj);
	void set_bnlSZ(int sz);
	void set_inlSZ(int sz);
	void set_interiorNL(boost::python::object obj);
	void set_inl(boost::python::object obj);
	void set_onl(boost::python::object obj);
	void set_snl(boost::python::object obj);
	void set_ndT(boost::python::object obj);
	void set_ndType(const int nt,LBM_DataHandler& f);
	void set_Ubc(const float u);
	void set_rhoBC(const float rho);
	void set_omega(const float o);
	void set_totalNodes(const int tn);
	void process_nodeList(const bool isEven,const int nodeList);
	void computeFout(LBM_DataHandler& f);
	void streamData(float * fOut,const int nd,LBM_DataHandler& f);
	int get_numSpd();
	int get_ndType();
	void registerNeighbor(const int ngbNumber,const int numData);
	void getHaloInPointers(boost::python::object lnd_num,
			boost::python::object spd, boost::python::object data, int ngb);
	void getHaloOutPointers(boost::python::object lnd_num,
			boost::python::object spd, boost::python::object data, int ngb);
	void extract_halo_data(bool isEven);
	void insert_boundary_data(bool isEven);
	LBM_DataHandler fData;
	Lattice * myLattice;

	LBM_HaloDataOrganizer myHalo_in;
	LBM_HaloDataOrganizer myHalo_out;


private:
	float * fIn;
	float * fOut;
	float * fEven;
	float * fOdd;
	int * adjacency;
	int * boundary_nl;
	int bnl_sz; //boundary node list size
	int * interior_nl;
	int inl_sz; // interior node list size
	int * inl;
	int * onl;
	int * snl;
	int * ndT;
	int numSpd;
	int totalNodes;

	float * ux;
	float * uy;
	float * uz;
	float * rho;


};


#endif
