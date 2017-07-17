/*
 * D3Q27Lattice.cpp
 *
 *  Created on: Jul 10, 2017
 *      Author: sblair
 */
#include "D3Q27Lattice.h"
#include <cstdlib>

D3Q27Lattice::D3Q27Lattice(const int Nx, const int Ny, const int Nz):
Lattice(Nx,Ny,Nz),
ex{0,-1,0,0,-1,-1,-1,-1,0,0,-1,-1,-1,-1,1,0,0,1,1,1,1,0,0,1,1,1,1},
ey{0,0,-1,0,-1,1,0,0,-1,-1,-1,-1,1,1,0,1,0,1,-1,0,0,1,1,1,1,-1,-1},
ez{0,0,0,-1,0,0,-1,1,-1,1,-1,1,-1,1,0,0,1,0,0,1,-1,1,-1,1,-1,1,-1},
w{8.f/27.f,2.f/27.f,2.f/27.f,2.f/27.f,1.f/54.f,1.f/54.f,1.f/54.f,1.f/54.f,1.f/54.f,
	1.f/54.f,1.f/216.f,1.f/216,1.f/216.f,1.f/216.f,2.f/27.f,2.f/27.f,
	2.f/27.f,1.f/54.f,1.f/54.f,1.f/54.f,1.f/54.f,1.f/54.f,
	1.f/54.f,1.f/216.f,1.f/216,1.f/216.f,1.f/216.f},
bbSpd{0,14,15,16,17,18,19,20,21,22,23,24,25,26,
		1,2,3,4,5,6,7,8,9,10,11,12,13},
Qflat{-1./3.,0,0,0,-1./3.,0,0,0,-1./3., // 0
	  2./3.,0,0,0,-1./3.,0,0,0,-1./3.,  // 1
	  -1./3.,0,0,0,2./3.,0,0,0,-1./3.,// 2
	  -1./3.,0,0,0,-1./3.,0,0,0,2./3.,// 3
	  2./3.,1,0,1,2./3.,0,0,0,-1./3.,// 4
	  2./3.,-1,0,-1,2./3.,0,0,0,-1./3.,// 5
	  2./3.,0,1,0,-1./3.,0,1,0,2./3.,// 6
	  2./3.,0,-1,0,-1./3.,0,-1,0,2./3.,// 7
	  -1./3.,0,0,0,2./3.,1,0,1,2./3.,//8
	  -1./3.,0,0,0,2./3.,-1,0,-1,2./3.,//9
	  2./3.,1,1,1,2./3.,1,1,1,2./3.,//10
	  2./3.,1,-1,1,2./3.,-1,-1,-1,2./3.,//11
	  2./3.,-1,1,-1,2./3.,-1,1,-1,2./3.,//12
	  2./3.,-1,-1,-1,2./3.,1,-1,1,2./3.,//13
	  2./3.,0,0,0,-1./3.,0,0,0,-1./3.,//14
	  -1./3.,0,0,0,2./3.,0,0,0,-1./3.,//15
	  -1./3.,0,0,0,-1./3.,0,0,0,2./3.,//16
	  2./3.,1,0,1,2./3.,0,0,0,-1./3.,//17
	  2./3.,-1,0,-1,2./3.,0,0,0,-1./3.,//18
	  2./3.,0,1,0,-1./3.,0,1,0,2./3.,//19
	  2./3.,0,-1,0,-1./3.,0,-1,0,2./3.,//20
	  -1./3.,0,0,0,2./3.,1,0,1,2./3.,//21
	  -1./3.,0,0,0,2./3.,-1,0,-1,2./3.,//22
	  2./3.,1,1,1,2./3.,1,1,1,2./3.,//23
	  2./3.,1,-1,1,2./3.,-1,-1,-1,2./3.,//24
	  2./3.,-1,1,-1,2./3.,-1,1,-1,2./3.,//25
	  2./3.,-1,-1,-1,2./3.,1,-1,1,2./3.}//26
{
	// direct base-class pointers to lattice variables
	setNumSpd(numSpd);
	setEx(ex);
	setEy(ey);
	setEz(ez);
	setW(w);
	setBBspd(bbSpd);
	setQflat(Qflat);

}

D3Q27Lattice::~D3Q27Lattice()
{

}

void D3Q27Lattice::set_inlet_bc_micro(LBM_DataHandler& f)
{
	int sp[9]={7,9,11,13,16,19,21,23,25};
	int bbSp[9]={20,22,24,26,3,6,8,10,12};
	int numBB = 9;
	for(int s=0;s<numBB;s++)
	{
		f.f[sp[s]]=f.fEq[sp[s]]+f.f[bbSp[s]]-f.fEq[bbSp[s]];
	}
}

void D3Q27Lattice::set_inlet_bc_macro(LBM_DataHandler& f)
{
	f.uz = f.u_bc;
	f.ux = 0; f.uy = 0.;
	f.rho = (1./(1. - f.uz))*(2.*(f.f[3]+f.f[6]+f.f[8]+
			f.f[10]+f.f[12]+f.f[20]+f.f[22]+f.f[24]+f.f[26])+
			(f.f[0]+f.f[1]+f.f[2]+f.f[4]+f.f[5]+
					f.f[14]+f.f[15]+f.f[17]+f.f[18]));
}

void D3Q27Lattice::set_outlet_bc_micro(LBM_DataHandler& f)
{
	int sp[9]={20,22,24,26,3,6,8,10,12};
	int bbSp[9]={7,9,11,13,16,19,21,23,25};
	int numBB = 9;
	for(int s=0;s<numBB;s++)
	{
		f.f[sp[s]]=f.fEq[sp[s]]+f.f[bbSp[s]]-f.fEq[bbSp[s]];
	}
}

void D3Q27Lattice::set_outlet_bc_macro(LBM_DataHandler& f)
{
	f.rho = f.rho_bc;
	f.uz = -1. + (1./f.rho)*(2.*
			(f.f[7]+f.f[9]+f.f[11]+f.f[13]+f.f[16]+f.f[19]+
					f.f[21]+f.f[23]+f.f[25])+
					(f.f[0]+f.f[1]+f.f[2]+f.f[4]+f.f[5]+
							f.f[14]+f.f[15]+f.f[17]+f.f[18]));
}


