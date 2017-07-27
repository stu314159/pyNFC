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

void PyLBM_Interface::registerNeighbor(const int ngbNum,const int numData)
{
	myHalo_in.insert_ngb(ngbNum,numData,numSpd);
	myHalo_out.insert_ngb(ngbNum,numData,numSpd);
}

void PyLBM_Interface::set_fIn(const float * fIn, const int nd, LBM_DataHandler& fData)
{
	for(int spd = 0; spd<numSpd; spd++)
	{
		fData.f[spd] = fIn[nd*numSpd+spd]; //keep an eye on how this data is organized
	}

}

void PyLBM_Interface::streamData(float * fOut, const int nd,LBM_DataHandler& fData)
{
	int tgtNode;
	for(int spd = 0; spd<numSpd; spd++)
	{
		tgtNode = adjacency[nd*numSpd+spd];
		fOut[tgtNode*numSpd+spd] = fData.fOut[spd];
	}

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

void PyLBM_Interface::set_ux(boost::python::object obj)
{
	PyObject* pobj = obj.ptr();
	Py_buffer pybuf;
	PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
	void * buf = pybuf.buf;
	ux = (float *)buf;
}

void PyLBM_Interface::set_uy(boost::python::object obj)
{
	PyObject* pobj = obj.ptr();
	Py_buffer pybuf;
	PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
	void * buf = pybuf.buf;
	uy = (float *)buf;
}

void PyLBM_Interface::set_uz(boost::python::object obj)
{
	PyObject* pobj = obj.ptr();
	Py_buffer pybuf;
	PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
	void * buf = pybuf.buf;
	uz = (float *)buf;
}

void PyLBM_Interface::set_rho(boost::python::object obj)
{
	PyObject* pobj = obj.ptr();
	Py_buffer pybuf;
	PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
	void * buf = pybuf.buf;
	rho = (float *)buf;
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


void PyLBM_Interface::set_ndT(boost::python::object obj)
{
	PyObject* pobj = obj.ptr();
	Py_buffer pybuf;
	PyObject_GetBuffer(pobj,&pybuf,PyBUF_SIMPLE);
	void * buf = pybuf.buf;
	ndT = (int *)buf;
}

void PyLBM_Interface::getHaloInPointers(boost::python::object nd,
		boost::python::object spd, boost::python::object data, int ngb)
{
	PyObject * pobj1 = nd.ptr();
	PyObject * pobj2 = spd.ptr();
	PyObject * pobj3 = data.ptr();

	Py_buffer pybuf1;
	Py_buffer pybuf2;
	Py_buffer pybuf3;

	PyObject_GetBuffer(pobj1,&pybuf1,PyBUF_SIMPLE);
	PyObject_GetBuffer(pobj2,&pybuf2,PyBUF_SIMPLE);
	PyObject_GetBuffer(pobj3,&pybuf3,PyBUF_SIMPLE);

	void * buf1 = pybuf1.buf;
	void * buf2 = pybuf2.buf;
	void * buf3 = pybuf3.buf;

	int * lnd_list = (int *) buf1;
	int * spd_list = (int *) buf2;
	float * dt = (float *) buf3;

	myHalo_in.initialize_ngb_pointers(ngb,lnd_list,spd_list,dt);


}

void PyLBM_Interface::getHaloOutPointers(boost::python::object nd,
		boost::python::object spd, boost::python::object data,int ngb)
{
	PyObject * pobj1 = nd.ptr();
		PyObject * pobj2 = spd.ptr();
		PyObject * pobj3 = data.ptr();

		Py_buffer pybuf1;
		Py_buffer pybuf2;
		Py_buffer pybuf3;

		PyObject_GetBuffer(pobj1,&pybuf1,PyBUF_SIMPLE);
		PyObject_GetBuffer(pobj2,&pybuf2,PyBUF_SIMPLE);
		PyObject_GetBuffer(pobj3,&pybuf3,PyBUF_SIMPLE);

		void * buf1 = pybuf1.buf;
		void * buf2 = pybuf2.buf;
		void * buf3 = pybuf3.buf;

		int * lnd_list = (int *) buf1;
		int * spd_list = (int *) buf2;
		float * dt = (float *) buf3;

		myHalo_out.initialize_ngb_pointers(ngb,lnd_list,spd_list,dt);

}


int PyLBM_Interface::get_numSpd()
{
	return numSpd;
}

int PyLBM_Interface::get_ndType()
{
	return fData.nodeType;
}

void PyLBM_Interface::computeFout(LBM_DataHandler & fData)
{
	myLattice->computeFout(fData);
}

void PyLBM_Interface::set_ndType(const int nd, LBM_DataHandler& fData)
{

	fData.nodeType = ndT[nd];

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

void PyLBM_Interface::set_dynamics(const int d)
{
	fData.dynamics = d;
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

void PyLBM_Interface::compute_local_data(const bool isEven)
{
	// get pointer to data
	float * fIn;
	if(isEven)
	{
		fIn = fEven;
	}else
	{
		fIn = fOdd;
	}
	// declare local iteration variables
	float ux_i, uy_i, uz_i,rho_i;
	float * f;
	int nd;
	//iterate through the boundary nodes
	for(int ndId = 0; ndId<bnl_sz; ndId++)
	{
		nd = boundary_nl[ndId]; //find local node number
		f = fIn+(nd*numSpd); //calculate f ptr
		//compute macroscopic data using myLattice object
		myLattice->computeMacroscopicData(rho_i,ux_i,uy_i,uz_i,f);
		// insert result into arrays
		ux[nd] = ux_i; uy[nd]=uy_i; uz[nd]=uz_i; rho[nd]=rho_i;
	}

	//iterate through the interior nodes
	for(int ndId = 0; ndId<inl_sz; ndId++)
	{
		nd = interior_nl[ndId]; //find local node number
		f = fIn+(nd*numSpd); //calculate f ptr
		//compute macroscopic data using myLattice object
		myLattice->computeMacroscopicData(rho_i,ux_i,uy_i,uz_i,f);
		// insert result into arrays
		ux[nd] = ux_i; uy[nd]=uy_i; uz[nd]=uz_i; rho[nd]=rho_i;
	}

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
#pragma omp parallel for
	for(int ndI=0; ndI<ndList_len;ndI++)
	{
		LBM_DataHandler fData_l(numSpd);
                fData_l.u_bc = fData.u_bc;
                fData_l.rho_bc = fData.rho_bc;
                fData_l.omega = fData.omega;
                // get the node number (for the local partition)
		int nd = ndList[ndI];
		// set the node type in fData
		set_ndType(nd,fData_l);
		// get the incoming data
		set_fIn(fIn,nd,fData_l);
		// compute fOut
		computeFout(fData_l); // passes fData to the appropriate lattice function and gets fOut
		// stream data to fOut array
		streamData(fOut,nd,fData_l);

	}
  
}


void PyLBM_Interface::extract_halo_data(bool isEven)
{
	float * f;
	if (isEven)
	{
		f = fOdd;
	}else
	{
		f = fEven;
	}
	myHalo_out.extractHaloData(f);


}

void PyLBM_Interface::insert_boundary_data(bool isEven)
{
	float * f;
	if (isEven)
	{
		f = fOdd;
	}else
	{
		f = fEven;
	}

	myHalo_in.distributeHaloData(f);


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
        		.def("set_ndT",&PyLBM_Interface::set_ndT)
        		.def("set_fEven",&PyLBM_Interface::set_fEven)
        		.def("set_fOdd",&PyLBM_Interface::set_fOdd)
        		.def("set_adjacency",&PyLBM_Interface::set_adjacency)
        		.def("set_boundaryNL",&PyLBM_Interface::set_boundaryNL)
        		.def("set_interiorNL",&PyLBM_Interface::set_interiorNL)
        		.def("set_bnlSZ",&PyLBM_Interface::set_bnlSZ)
        		.def("set_inlSZ",&PyLBM_Interface::set_inlSZ)
        		.def("set_totalNodes",&PyLBM_Interface::set_totalNodes)
        		.def("process_nodeList",&PyLBM_Interface::process_nodeList)
        		.def("set_ux",&PyLBM_Interface::set_ux)
        		.def("set_uy",&PyLBM_Interface::set_uy)
        		.def("set_uz",&PyLBM_Interface::set_uz)
        		.def("set_rho",&PyLBM_Interface::set_rho)
        		.def("compute_local_data",&PyLBM_Interface::compute_local_data)
        		.def("registerNeighbor",&PyLBM_Interface::registerNeighbor)
        		.def("getHaloOutPointers",&PyLBM_Interface::getHaloOutPointers)
        		.def("getHaloInPointers",&PyLBM_Interface::getHaloInPointers)
        		.def("extract_halo_data",&PyLBM_Interface::extract_halo_data)
        		.def("insert_boundary_data",&PyLBM_Interface::insert_boundary_data)
        		.def("set_dynamics",&PyLBM_Interface::set_dynamics)
        		;
}
