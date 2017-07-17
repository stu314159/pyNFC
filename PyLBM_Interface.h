#ifndef PYLBM_INTERFACE_H_
#define PYLBM_INTERFACE_H_

#include "Lattice.h"
#include "D3Q15Lattice.h"
#include "D3Q19Lattice.h"
#include "D3Q27Lattice.h"
#include "LBM_DataHandler.h"

#include <boost/python.hpp>

class PyLBM_Interface
{
public:
  PyLBM_Interface(const int numSpd);
  ~PyLBM_Interface();
  void set_fIn(boost::python::object obj);
  void get_fOut(boost::python::object obj);
  void set_fEven(boost::python::object obj);
  void set_fOdd(boost::python::object obj);
  void set_adjacency(boost::python::object obj);
  void set_boundaryNL(boost::python::object obj);
  void set_bnlSZ(int sz);
  void set_inlSZ(int sz);
  void set_interiorNL(boost::python::object obj);
  void set_inl(boost::python::object obj);
  void set_onl(boost::python::object obj);
  void set_snl(boost::python::object obj);
  void set_ndType(const int nt);
  void set_Ubc(const float u);
  void set_rhoBC(const float rho);
  void set_omega(const float o);
  void set_totalNodes(const int tn);
  void computeFout();
  int get_numSpd();
  int get_ndType();
  LBM_DataHandler fData;
  Lattice * myLattice;


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
  int numSpd;
  int totalNodes;
   

};


#endif
