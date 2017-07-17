#include "Lattice.h"
#include <cstdlib>
#include <iostream> //<-- used for debugging

Lattice::Lattice(const int Nx, const int Ny, const int Nz):
Nx(Nx), Ny(Ny), Nz(Nz), numSpd(0),

ex(NULL), ey(NULL), ez(NULL),w(NULL), bbSpd(NULL), Qflat(NULL)
{

}

Lattice::~Lattice()
{


}

void Lattice::bounceBack(LBM_DataHandler& f)
{
	float fTmp[numSpd];
	for(int spd=0;spd<numSpd;spd++)
	{
		fTmp[spd]=f.f[bbSpd[spd]];
	}
	for(int spd=0;spd<numSpd;spd++)
	{
		f.f[spd]=fTmp[spd];
	}

}
void Lattice::computeMacroscopicData(float& rho, float& ux, float& uy, float& uz, const float* f)
{
	rho = 0.; ux = 0.; uy = 0.; uz = 0.;
	for(int spd = 0; spd<numSpd; spd++)
	{
		rho+=f[spd];
		ux+=ex[spd]*f[spd];
		uy+=ey[spd]*f[spd];
		uz+=ez[spd]*f[spd];
	}
	ux/=rho; uy/=rho; uz/=rho;

}
void Lattice::computeMacroscopicData(LBM_DataHandler& f)
{
	float rho, ux, uy, uz;
	rho = 0.; ux = 0.; uy = 0.; uz = 0.;
	for(int spd = 0; spd<numSpd; spd++)
	{

		rho+=f.f[spd];
		ux+=ex[spd]*f.f[spd];
		uy+=ey[spd]*f.f[spd];
		uz+=ez[spd]*f.f[spd];
	}
	ux/=rho; uy/=rho; uz/=rho;
	f.rho = rho; f.ux = ux; f.uy = uy; f.uz = uz;

}

void Lattice::computeEquilibrium(float * fEq,const float ux, const float uy, const float uz, const float rho)
{
	//this will assume that up-to-date data exists in rho, ux, uy, and uz class variables
	//and that ex,ey,ez,and w have been updated with the lattice-specific info
	float cu;

	for(int spd = 0; spd<numSpd;spd++)
	{
		fEq[spd] = 0;
		cu = 3.f*(ex[spd]*ux+ey[spd]*uy+ez[spd]*uz);
		fEq[spd]=w[spd]*rho*(1.f+cu+0.5f*(cu*cu)-3.f/2.f*(ux*ux+uy*uy+uz*uz));

	}
}

void Lattice::computeEquilibrium(LBM_DataHandler& f)
{

	computeEquilibrium(f.fEq,f.ux,f.uy,f.uz,f.rho);
}

void Lattice::compute_piFlat(LBM_DataHandler& f)
{

	//f.piFlat = {0,0,0,0,0,0,0,0,0}; // initialized in constructor
	for(int k=0; k<9;k++)
	{
		f.piFlat[k]=0.;
	}
	float fNeq;
	for(int spd = 0; spd<numSpd; spd++)
	{
		fNeq = f.f[spd] - f.fEq[spd];
		f.piFlat[0]+=ex[spd]*ex[spd]*fNeq;
		f.piFlat[1]+=ey[spd]*ex[spd]*fNeq;
		f.piFlat[2]+=ez[spd]*ex[spd]*fNeq;
		f.piFlat[3]+=ex[spd]*ey[spd]*fNeq;
		f.piFlat[4]+=ey[spd]*ey[spd]*fNeq;
		f.piFlat[5]+=ez[spd]*ey[spd]*fNeq;
		f.piFlat[6]+=ex[spd]*ez[spd]*fNeq;
		f.piFlat[7]+=ey[spd]*ez[spd]*fNeq;
		f.piFlat[8]+=ez[spd]*ez[spd]*fNeq;

	}
}

void Lattice::regularize(LBM_DataHandler& f)
{

	float wa;

	for(int spd = 0; spd<numSpd; spd++)
	{
		// get leading constant
		wa = (9.f/2.f)*w[spd];
		f.f[spd]=f.fEq[spd];
		// load chunk of Qflat:
		for(int k=0;k<9;k++)
		{
			f.f[spd]+= wa*Qflat[spd*9+k]*f.piFlat[k];
		}
	}
}

void Lattice::relax(LBM_DataHandler& f)
{

	for(int spd=0;spd<numSpd;spd++)
	{
		f.fOut[spd]=f.f[spd]-f.omega*(f.f[spd] - f.fEq[spd]);
	}

}

void Lattice::computeFout(LBM_DataHandler& f)
{
	// node type 1: just bounce back
	// node type 0, 2, and 3 continue with the following steps:

	// compute macroscopic velocity and pressure
	computeMacroscopicData(f);
	if(f.nodeType==1)
	{
		f.ux = 0.; f.uy = 0.; f.uz = 0; // solid nodes, zero velocity
	}

	// node type 2 and 3 apply macroscopic boundary conditions
	if(f.nodeType==2) //inlet node
	{
		set_inlet_bc_macro(f);
	}
	if(f.nodeType==3) // outlet node
	{
		set_outlet_bc_macro(f);
	}

	// compute equilibrium
	computeEquilibrium(f);

	// node type 2 and 3 apply microscopic boundary conditions and regularization
	if(f.nodeType==2)
	{
		set_inlet_bc_micro(f);
	}
	if(f.nodeType==3)
	{
		set_outlet_bc_micro(f);
	}

	// get (flattened) second-order moment of particle density distribution
	compute_piFlat(f);
	regularize(f);
	relax(f);

}
