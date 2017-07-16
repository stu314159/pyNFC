#include "PyLBM_Interface.h"


PyLBM_Interface::PyLBM_Interface(const int numSpd) :
numSpd(numSpd),fIn(NULL),fOut(NULL),u_bc(0),rho_bc(0),omega(0),nd_type(0)
{
  

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
}

int PyLBM_Interface::get_numSpd()
{
  return numSpd;
}


using namespace boost::python;

BOOST_PYTHON_MODULE(LBM_Interface)
{
    class_<PyLBM_Interface>("PyLBM_Interface",init<int>())
        .def("set_fIn",&PyLBM_Interface::set_fIn)
        .def("get_numSpd",&PyLBM_Interface::get_numSpd)
     ;
}
