/*
 * LBM_DataHandler.cpp
 *
 *  Created on: Jul 10, 2017
 *      Author: stu
 */


#include "LBM_DataHandler.h"

LBM_DataHandler::LBM_DataHandler(const int numSpd) :
ux(0),uy(0),uz(0),rho(0),u_bc(0),rho_bc(0),nodeType(0),omega(0),
piFlat{0,0,0,0,0,0,0,0,0},numSpd(numSpd)
{
	f = new float[numSpd];
	fEq = new float[numSpd];
	fOut = new float[numSpd];


}

LBM_DataHandler::~LBM_DataHandler()
{
	delete [] f;
	delete [] fEq;
	delete [] fOut;

}

int LBM_DataHandler::get_numSpd()
{
  return numSpd;
}

void LBM_DataHandler::set_fIn(boost::python::object obj)
{
  PyObject* pobj = obj.ptr();
  Py_buffer pybuf;
  PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
  void * buf = pybuf.buf;
  f = (float *)buf;
  PyBuffer_Release(&pybuf);
}


void LBM_DataHandler::multFin(float mult)
{
   for(int i=0;i<numSpd;i++)
   {
      f[i]=mult*f[i];
   }
}



using namespace boost::python;

BOOST_PYTHON_MODULE(LDH)
{
    class_<LBM_DataHandler>("LBM_DataHandler",init<int>())
        .def("get_numSpd",&LBM_DataHandler::get_numSpd)
        .def("set_fIn",&LBM_DataHandler::set_fIn)
        .def("multFin",&LBM_DataHandler::multFin)
    ;
}
