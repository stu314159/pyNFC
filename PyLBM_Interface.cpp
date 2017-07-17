#include "PyLBM_Interface.h"


PyLBM_Interface::PyLBM_Interface(const int numSpd) :
numSpd(numSpd),fData(numSpd),fIn(NULL),fOut(NULL)
{
  // create the lattice object
  switch (numSpd)
  {
  case(15):
	myLattice = new D3Q15Lattice(0,0,0); break;
  case(19):
	myLattice = new D3Q19Lattice(0,0,0); break;
  case(27):
	myLattice = new D3Q27Lattice(0,0,0);

  }

}

PyLBM_Interface::~PyLBM_Interface()
{
  
  
}

void PyLBM_Interface::set_fIn(const float * fIn, const int nd)
{
  for(int spd = 0; spd<numSpd; spd++)
  {
    fData.f[spd] = fIn[nd*numSpd+spd]; //keep an eye on how this data is organized
  }

}

void PyLBM_Interface::streamData(float * fOut, const int nd)
{
  int tgtNode;
  for(int spd = 0; spd<numSpd; spd++)
  {
    tgtNode = adjacency[nd*numSpd+spd];
    fOut[tgtNode*numSpd+spd] = fData.fOut[spd];
  }

}

//void PyLBM_Interface::set_fIn(boost::python::object obj)
//{
//  PyObject* pobj = obj.ptr();
//  Py_buffer pybuf;
//  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
//  void * buf = pybuf.buf;
//  fIn = (float *)buf;
//  PyBuffer_Release(&pybuf);
//
//  // initialize fData with incoming data
//  for(int spd=0;spd<numSpd;spd++)
//  {
//    fData.f[spd] = fIn[spd];
//  }
//}

void PyLBM_Interface::set_inl(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  inl = (int *)buf;
}

void PyLBM_Interface::set_fEven(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  fEven = (float *)buf;
}

void PyLBM_Interface::set_fOdd(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  fOdd = (float *)buf;
}

void PyLBM_Interface::set_adjacency(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  adjacency = (int *)buf;
}

void PyLBM_Interface::set_boundaryNL(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  boundary_nl = (int *)buf;
}

void PyLBM_Interface::set_interiorNL(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  interior_nl = (int *)buf;
}

void PyLBM_Interface::set_onl(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  onl = (int *)buf;
}

void PyLBM_Interface::set_snl(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  snl = (int *)buf;
}

//void PyLBM_Interface::get_fOut(boost::python::object obj)
//{
//  PyObject* pobj = obj.ptr();
//  Py_buffer pybuf;
//  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
//  void * buf = pybuf.buf;
//  fOut = (float *)buf;
//  PyBuffer_Release(&pybuf);
//
//  //load fOut with fData.fOut
//  for(int spd=0;spd<numSpd;spd++)
//  {
//    fOut[spd]=fData.fOut[spd];
//  }
//
//}

int PyLBM_Interface::get_numSpd()
{
  return numSpd;
}

int PyLBM_Interface::get_ndType()
{
  return fData.nodeType;
}

void PyLBM_Interface::computeFout()
{
  myLattice->computeFout(fData);
}

void PyLBM_Interface::set_ndType(const int nd)
{
  fData.nodeType = 0;
  if(snl[nd]==1)
  {
    fData.nodeType=1;
  }else if(inl[nd]==1)
  {
    fData.nodeType=2;
  }else if(onl[nd]==1)
  {
    fData.nodeType=3;
  }
  
}

void PyLBM_Interface::set_Ubc(const float u)
{
  fData.u_bc = u;
}

void PyLBM_Interface::set_rhoBC(const float rho)
{
  fData.rho_bc = rho;
}

void PyLBM_Interface::set_omega(const float o)
{
  fData.omega = o;
}

void PyLBM_Interface::set_bnlSZ(int sz)
{
  bnl_sz = sz;
}

void PyLBM_Interface::set_inlSZ(int sz)
{
  inl_sz = sz;
}

void PyLBM_Interface::set_totalNodes(const int tn)
{
  totalNodes = tn;
}

void PyLBM_Interface::process_nodeList(const bool isEven,const int nodeListnum)
{
  if(isEven)
  {
   fIn = fEven; fOut = fOdd;
  }else
  {
   fIn = fOdd; fOut = fEven;
  }

  //node list 0 is the boundary node list; node list 1 is the interior node list
  int * ndList;
  int ndList_len;
  if(nodeListnum == 0)//boundary node list
  {
    ndList = boundary_nl; ndList_len = bnl_sz;
  }else if(nodeListnum == 1)
  {
    ndList = interior_nl; ndList_len = inl_sz;
  }

  for(int ndI=0; ndI<ndList_len;ndI++)
  {
     // get the node number (for the local partition)
     int nd = ndList[ndI];
     // set the node type in fData
     set_ndType(nd); 
     // get the incoming data
     set_fIn(fIn,nd);
     // compute fOut
     computeFout(); // passes fData to the appropriate lattice function and gets fOut
     // stream data to fOut array
     streamData(fOut,nd);

  }



}

using namespace boost::python;

BOOST_PYTHON_MODULE(LBM_Interface)
{
    class_<PyLBM_Interface>("PyLBM_Interface",init<int>())
        .def("get_numSpd",&PyLBM_Interface::get_numSpd)
        .def("computeFout",&PyLBM_Interface::computeFout)
        .def("set_ndType",&PyLBM_Interface::set_ndType)
        .def("get_ndType",&PyLBM_Interface::get_ndType)
        .def("set_Ubc",&PyLBM_Interface::set_Ubc)
        .def("set_rhoBC",&PyLBM_Interface::set_rhoBC)
        .def("set_omega",&PyLBM_Interface::set_omega)
        .def("set_inl",&PyLBM_Interface::set_inl)
        .def("set_onl",&PyLBM_Interface::set_onl)
        .def("set_snl",&PyLBM_Interface::set_snl)
        .def("set_fEven",&PyLBM_Interface::set_fEven)
        .def("set_fOdd",&PyLBM_Interface::set_fOdd)
        .def("set_adjacency",&PyLBM_Interface::set_adjacency)
        .def("set_boundaryNL",&PyLBM_Interface::set_boundaryNL)
        .def("set_interiorNL",&PyLBM_Interface::set_interiorNL)
        .def("set_bnlSZ",&PyLBM_Interface::set_bnlSZ)
        .def("set_inlSZ",&PyLBM_Interface::set_inlSZ)
        .def("set_totalNodes",&PyLBM_Interface::set_totalNodes)
        .def("process_nodeList",&PyLBM_Interface::process_nodeList)
     ;
}
