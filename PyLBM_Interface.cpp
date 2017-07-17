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

void PyLBM_Interface::set_inl(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  inl = (int *)buf;
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
        .def("set_Ubc",&PyLBM_Interface::set_Ubc)
        .def("set_rhoBC",&PyLBM_Interface::set_rhoBC)
        .def("set_omega",&PyLBM_Interface::set_omega)
        .def("set_inl",&PyLBM_Interface::set_inl)
        .def("set_onl",&PyLBM_Interface::set_onl)
        .def("set_snl",&PyLBM_Interface::set_snl)
     ;
}
