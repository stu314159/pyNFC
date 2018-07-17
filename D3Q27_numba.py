#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 08:32:46 2018

@author: stu
"""

import numpy as np
import numba
from numba import cuda
import math

ex27 = np.array([0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1,0,0,0,0,1,1,1,1,-1,-1,-1,-1],
              dtype=np.float32);
ey27 = np.array([0,0,0,1,-1,0,0,1,-1,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1],
              dtype=np.float32);
ez27 = np.array([0,0,0,0,0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1,-1],
              dtype=np.float32);
              
w27 = np.array([8./27.,2./27.,2./27.,2./27.,2./27.,2./27.,2./27.,
                  1./54.,1./54.,1./54.,1./54.,1./54.,1./54.,
                  1./54.,1./54.,1./54.,1./54.,1./54.,1./54.,
                  1./216.,1./216.,1./216.,1./216.,
                  1./216.,1./216.,1./216.,1./216.], dtype=np.float32);

bbSpd27 = np.array([0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15,26,25,24,23,22,21,20,19],
                 dtype=np.int32);
                   
numSpd = 27;

                   
@cuda.jit('void(float32[:],float32,float32,float32,float32,float32[:],float32[:],float32[:],float32[:])',device=True)
def compute_f_eq(f_eq,rho,ux,uy,uz,ex,ey,ez,w):
    """
    
    """
    for spd in range(numSpd):
        cu = 3.*(ex[spd]*ux + ey[spd]*uy + ez[spd]*uz);
        f_eq[spd] = w[spd]*rho*(1. + cu + 0.5*(cu*cu) 
        - 3./2.*(ux*ux + uy*uy + uz*uz));
    
@cuda.jit('float32(float32[:])',device=True)               
def compute_rho(f):
    rho = 0.;
    for i in range(numSpd):
        rho = rho + f[i]
        
    return rho

@cuda.jit('float32(float32[:],float32[:],float32)',device=True)    
def compute_u(f,e,rho):
    """
    compute macroscopic velocity in the x-direction
           
    """
    u = 0.;
    
    for i in range(numSpd):
        u = u + f[i]*e[i];
        
    return u/rho;

@cuda.jit('void(float32[:],float32[:],int32[:])',device=True)
def bounce_back(f_out,f_in,bbSpd):
    """
    bounce-back for solid nodes
    
    """
    for spd in range(numSpd):
        f_out[bbSpd[spd]] = f_in[spd];

@cuda.jit('float32(float32[:],float32)',device=True)     
def set_inlet_bc_macro(f,uz):
    """
    compute rho based on microscopic variables and given bc
    """
    rho = (1./(1.-uz))*(2.*(f[6]+f[14]+f[12]+
            f[18]+f[16]+f[26]+f[24]+f[22]+f[20])+
            (f[0]+f[1]+f[2]+f[3]+f[4]+
             f[7]+f[8]+f[9]+f[10]));
    return rho

@cuda.jit('float32(float32[:],float32)',device=True)
def set_outlet_bc_macro(f,rho):
    """
    compute macroscopic outlet parameters
    """
    uz = -1. + (1./rho)*(2.*
               (f[5]+f[11]+f[13]+f[17]+f[19]+f[21]+f[23]+f[25]) +
               (f[0]+f[1]+f[2]+f[3]+f[4]+f[7]+f[8]+f[9]+f[10]));
    return uz;

@cuda.jit('void(float32[:],float32[:],float32[:],float32[:,:])',device=True)
def compute_Pi1_flat(Pi1_flat,f_in,f_eq,Qflat):
    """
    compute the 3x3 (flattened) array required for the regularization processes
    """
    # initialize Pi1_flat to zero
    for i in range(9):
        Pi1_flat[i] = 0.
    
    for spd in range(numSpd):
        for k in range(9):
            Pi1_flat[k]+=(f_in[spd] - f_eq[spd])*Qflat[spd,k];
        
    
@cuda.jit('void(float32[:],float32[:],float32[:],float32[:,:],float32[:])',device=True)
def regularize_boundary_nodes(f_in, f_eq, Pi1_flat,Qflat,w):
    """
    apply regularization process to inlet and outlet boundary nodes as addapted
    from J. Latt dissertation
    """
    for spd in range(numSpd):
        ft = 0.;
        for k in range(9):
            ft += w[spd]*(9./2.)*Qflat[spd,k]*Pi1_flat[k];
        f_in[spd] = f_eq[spd]+ft;

@cuda.jit('void(float32[:,:],float32[:],float32[:],float32[:],float32[:],float32[:])',device=True)
def compute_strain_tensor(S,ex,ey,ez,f_in,f_eq):
    """
    generates the 3x3 viscous strain tensor for use with the
    1-parameter turbulence model
    """
    for i in range(3):
        for j in range(3):
            S[i,j] = 0.;
    
    for spd in range(numSpd):
        e = (ex[spd],ey[spd],ez[spd]); # I hope this works. (nope)
        for i in range(3):
            for j in range(3):
                S[i,j] += e[i]*e[j]*(f_in[spd] - f_eq[spd]);


@cuda.jit('float32(float32,float32,float32[:,:])',device=True)
def apply_turbulence_model(omega,Cs,S):
    """
    1-parameter Smagorinsky-type turbulence model.  (needs validation)
    """
    nu = ((1./omega) - 0.5)*3.0;
    P = math.sqrt(S[0,0]**2 + S[1,1]**2 + S[2,2]**2 + 
                2.0*(S[0,1]**2 + S[0,2]**2 + S[1,2]**2));
    P *= Cs; P = math.sqrt(P + nu**2) - nu;
    nu_e = P/6.;
    return 1./(3.*(nu+nu_e)+0.5)

@cuda.jit('void(float32[:],float32[:],float32)',device=True)
def relax(f_in, f_eq, omega):
    """
    relax toward equilibrium
    """
    
    for spd in range(numSpd):
        f_in[spd] = f_in[spd] - omega*(f_in[spd] - f_eq[spd]);
        

@cuda.jit('void(float32[:,:],int32[:,:],float32[:],int32)',device=True)
def stream(fOut,adjArray,f_out,tid):
    """
    execute the LBM streaming operation
    """
    for spd in range(numSpd):
        tgtNd = adjArray[tid,spd];
        fOut[tgtNd,spd] = f_out[spd];

@cuda.jit('void(float32[:,:],float32[:,:],int32[:,:],float32[:,:],int32[:],float32,float32,float32,float32,float32[:,:],int32[:],int32)')
def process_node_list(fOut,fIn,adjArray,MacroV,ndType,u_bc,rho_lbm,omega,Cs,Qflat27,theList,N):
    """
    a D3Q27-specific kernel to process LBM nodes
    """
    #ex = cuda.shared.array(27,dtype=numba.float32);
    #ey = cuda.shared.array(27,dtype=numba.float32);
    #ez = cuda.shared.array(27,dtype=numba.float32);
    #w = cuda.shared.array(27,dtype=numba.float32);
    #bbSpd = cuda.shared.array(27,dtype=numba.int32);
    sQflat = cuda.shared.array((27,9),dtype=numba.float32);
    
    ex = cuda.const.array_like(ex27);
    ey = cuda.const.array_like(ey27);
    ez = cuda.const.array_like(ez27);
    w = cuda.const.array_like(w27);
    bbSpd = cuda.const.array_like(bbSpd27);
    #Qflat = cuda.const.array_like(Qflat27); # numSpd x 9 flattened tensor <-- cannot be constant
    
    
# figure out which thread I am...
    tx = cuda.threadIdx.x;
    bx = cuda.blockIdx.x;
    bw = cuda.blockDim.x;
    tid = tx + bx*bw;
    
    if (tid < N):
        # get fIn data
        f_in = cuda.local.array(numSpd,dtype=numba.float32);
        f_out = cuda.local.array(numSpd,dtype=numba.float32);
        S = cuda.local.array((3,3),dtype=numba.float32);
        nd = theList[tid];
        for i in range(numSpd):
            f_in[i] = fIn[nd,i];
            
        ## load constant data into shared arrays
        if (tid < 27):
            for k in range(9):
                sQflat[tid,k]=Qflat27[tid,k];
        #    ex[tx] = ex_c[tx];
        #    ey[tx] = ey_c[tx];
        #    ez[tx] = ez_c[tx];
        #    w[tx] = w_c[tx];
        #    bbSpd[tx] = bbSpd_c[tx];
        #    for s in range(9):
        #        Qflat[tx,s] = Qflat_c[tx,s]
        
        # compute macroscopic data
        rho = compute_rho(f_in);
        ux = compute_u(f_in,ex,rho);
        uy = compute_u(f_in,ey,rho);
        uz = compute_u(f_in,ez,rho);
                     
        
               
        # get the node type
        ndT = ndType[nd];
        
        # set macroscopic BCs:
        if ndT == 1: # solid node
            ux = 0.; uy = 0.; uz = 0.; # consider storing this for return
            bounce_back(f_out,f_in,bbSpd); # just bounce-back and stream
        elif (ndT == 2) or (ndT == 3) or (ndT == 5):
            if (ndT == 2): # inlet node type
                uz = u_bc; ux = 0.; uy = 0.;
                rho = set_inlet_bc_macro(f_in,uz); #re-assigns rho
            if (ndT ==3): # outlet bc node
                rho = rho_lbm; 
                uz = set_outlet_bc_macro(f_in,rho); # re-assign uz
                
        if ndT != 1: # everyone but solid nodes
            # compute equilibribum
            f_eq = cuda.local.array(27,dtype=numba.float32)
            if (ndT == 2) or (ndT == 3): # regularize boundary nodes
                Pi1_flat = cuda.local.array(9,dtype=numba.float32);
                compute_Pi1_flat(Pi1_flat,f_in,f_eq,sQflat);
                regularize_boundary_nodes(f_in, f_eq, Pi1_flat,sQflat,w);
            compute_strain_tensor(S,ex,ey,ez,f_in,f_eq);
            omega_t = apply_turbulence_model(omega,Cs,S);
            f_out = relax(f_in,f_eq,omega_t);
        
            
        
        # save macroscopic dependent variables
        MacroV[nd,0] = rho;
        MacroV[nd,1] = ux;
        MacroV[nd,2] = uy;
        MacroV[nd,3] = uz;
        
        # stream data to neighbors in the fOut array
        stream(fOut,adjArray,f_out,nd);






