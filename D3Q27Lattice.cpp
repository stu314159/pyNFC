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
ex{0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1,0,0,0,0,1,1,1,1,-1,-1,-1,-1},
ey{0,0,0,1,-1,0,0,1,-1,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1},
ez{0,0,0,0,0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1},
w{8./27.,2./27.,2./27.,2./27.,2./27.,2./27.,2./27.,
    1./54.,1./54.,1./54.,1./54.,1./54.,1./54.,
    1./54.,1./54.,1./54.,1./54.,1./54.,1./54.,
    1./216.,1./216.,1./216.,1./216.,
    1./216.,1./216.,1./216.,1./216.},
bbSpd{0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15,26,25,24,23,22,21,20,19},
Qflat{-1./3.,0,0,0,-1./3.,0,0,0,-1./3., // 0
	  2./3.,0,0,0,-1./3.,0,0,0,-1./3.,  // 1
	  2./3.,0,0,0,-1./3.,0,0,0,-1./3.,// 2
	  -1./3.,0,0,0,2./3.,0,0,0,-1./3.,// 3
	  -1./3.,0,0,0,2./3.,0,0,0,-1./3.,// 4
	  -1./3.,0,0,0,-1./3.,0,0,0,2./3.,// 5
	  -1./3.,0,0,0,-1./3.,0,0,0,2./3.,// 6
	  2./3.,1,0,1,2./3.,0,0,0,-1./3.,// 7
	  2./3.,-1,0,-1,2./3.,0,0,0,-1./3.,//8
	  2./3.,-1,0,-1,2./3.,0,0,0,-1./3.,//9
	  2./3.,1,0,1,2./3.,0,0,0,-1./3.,//10
	  2./3.,0,1,0,-1./3.,0,1,0,2./3.,//11
	  2./3.,0,-1,0,-1./3.,0,-1,0,2./3.,//12
	  2./3.,0,-1,0,-1./3.,0,-1,0,2./3.,//13
	  2./3.,0,1,0,-1./3.,0,1,0,2./3.,//14
	  -1./3.,0,0,0,2./3.,1,0,1,2./3.,//15
	  -1./3.,0,0,0,2./3.,-1,0,-1,2./3.,//16
	  -1./3.,0,0,0,2./3.,-1,0,-1,2./3.,//17
	  -1./3.,0,0,0,2./3.,1,0,1,2./3.,//18
	  2./3.,1,1,1,2./3.,1,1,1,2./3.,//19
	  2./3.,1,-1,1,2./3.,-1,-1,-1,2./3.,//20
	  2./3.,-1,1,-1,2./3.,-1,1,-1,2./3.,//21
	  2./3.,-1,-1,-1,2./3.,1,-1,1,2./3.,//22
	  2./3.,-1,-1,-1,2./3.,1,-1,1,2./3.,//23
	  2./3.,-1,1,-1,2./3.,-1,1,-1,2./3.,//24
	  2./3.,1,-1,1,2./3.,-1,-1,-1,2./3.,//25
	  2./3.,1,1,1,2./3.,1,1,1,2./3.}//26
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
	int sp[9]={5,11,13,15,17,19,21,23,25};
	int bbSp[9]={6,14,12,18,16,26,24,22,20};
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
	f.rho = (1./(1. - f.uz))*(2.*(f.f[6]+f.f[14]+f.f[12]+
			f.f[18]+f.f[16]+f.f[26]+f.f[24]+f.f[22]+f.f[20])+
			(f.f[0]+f.f[1]+f.f[2]+f.f[3]+f.f[4]+
					f.f[7]+f.f[8]+f.f[9]+f.f[10]));
}

void D3Q27Lattice::set_outlet_bc_micro(LBM_DataHandler& f)
{
	int sp[9]={6,14,12,18,16,26,24,22,20};
	int bbSp[9]={5,11,13,15,17,19,21,23,25};
	int numBB = 9;
	for(int s=0;s<numBB;s++)
	{
		f.f[sp[s]]=f.fEq[sp[s]]+f.f[bbSp[s]]-f.fEq[bbSp[s]];
	}
}

void D3Q27Lattice::set_outlet_bc_macro(LBM_DataHandler& f)
{
	f.rho = f.rho_bc;
        //f.ux = 0.; f.uy = 0.;
	f.uz = -1. + (1./f.rho)*(2.*
			(f.f[5]+f.f[11]+f.f[13]+f.f[15]+f.f[17]+f.f[19]+
					f.f[21]+f.f[23]+f.f[25])+
					(f.f[0]+f.f[1]+f.f[2]+f.f[3]+f.f[4]+
							f.f[7]+f.f[8]+f.f[9]+f.f[10]));
}


