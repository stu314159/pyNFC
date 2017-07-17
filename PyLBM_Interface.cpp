#include "PyLBM_Interface.h"


PyLBM_Interface::PyLBM_Interface(const int numSpd) :
numSpd(numSpd),fData(numSpd),fIn(NULL),fOut(NULL),u_bc(0),rho_bc(0),omega(0),nd_type(0)
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

void PyLBM_Interface::set_fIn(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  fIn = (float *)buf;
  PyBuffer_Release(&pybuf);

  // initialize fData with incoming data
  for(int spd=0;spd<numSpd;spd++)
  {
    fData.f[spd] = fIn[spd];
  }
}


void PyLBM_Interface::get_fOut(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  fOut = (float *)buf;
  PyBuffer_Release(&pybuf);

  //load fOut with fData.fOut
  for(int spd=0;spd<numSpd;spd++)
  {
    fOut[spd]=fData.fOut[spd];
  }

}

int PyLBM_Interface::get_numSpd()
{
  return numSpd;
}

int PyLBM_Interface::get_ndType()
{
  return nd_type;
}

void PyLBM_Interface::computeFout()
{
  myLattice->computeFout(fData);
}

void PyLBM_Interface::set_ndType(const int nt)
{
  nd_type = nt;
}


using namespace boost::python;

BOOST_PYTHON_MODULE(LBM_Interface)
{
    class_<PyLBM_Interface>("PyLBM_Interface",init<int>())
        .def("set_fIn",&PyLBM_Interface::set_fIn)
        .def("get_numSpd",&PyLBM_Interface::get_numSpd)
        .def("get_fOut",&PyLBM_Interface::get_fOut)
        .def("computeFout",&PyLBM_Interface::computeFout)
        .def("set_ndType",&PyLBM_Interface::set_ndType)
        .def("get_ndType",&PyLBM_Interface::get_ndType)
     ;
}
