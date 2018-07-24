#partition_compare.py
"""
Generate specific problem and generate
partitions using my own partition_suggestion.py
functions as well as using pymetis.

Try to see why one might be better than the other.
"""

#import argparse
import math
import numpy as np
import numba
#from numba import cuda

NO_PYMETIS=0
try:
  from pymetis import part_graph
except ImportError:
  NO_PYMETIS=1
  import sys

#from vtkHelper import saveScalarStructuredGridVTK_binary
#from vtkHelper import saveVelocityAndPressureVTK_binary
from vtkHelper import saveStructuredPointsVTK_ascii

from partition_suggestion import part_advisor
import time


class Partition:
    def __init__(self,xmin,xmax,ymin,ymax,zmin,zmax):
        self.xmin = xmin;
        self.ymin = ymin;
        self.zmin = zmin;
        self.xmax = xmax;
        self.ymax = ymax;
        self.zmax = zmax;

    def inPart(self,x,y,z):
        if ((x < self.xmin) or (x > self.xmax)):
            return False
        if((y < self.ymin) or (y > self.ymax)):
            return False
        if((z < self.zmin) or (z > self.zmax)):
            return False
        return True

def set_geometric_partition(Nx,Ny,Nz,px,py,pz):
    """
    Nx = integer, number of lattice points in the x direction
    Ny = integer, number of lattice points in the y direction
    Nz = integer, number of lattice points in the z direction
    px,py,pz = integers: number of partitions in the x,y and z direction

    returns part_vert for a geometric partition

    """
    partList = []
    bx = Nx/px; by = Ny/py; bz = Nz/pz # block x,y and z minimum sizes
    #part = 0
    zmin = 0
    for z in range(pz):
        zmax = zmin + bz-1
        if ((Nz%pz)>z):
          zmax+=1
        ymin = 0
        for y in range(py):
            ymax = ymin+by-1
            if((Ny%py)>y):
                ymax+=1
            xmin = 0
            for x in range(px):
                xmax = xmin+bx-1
                if((Nx%px)>x):
                    xmax+=1
                partList.append(Partition(xmin,xmax,ymin,ymax,zmin,zmax))
                xmin = xmax+1
            ymin = ymax+1
        zmin = zmax+1

    # now I have my list of partitions, cycle through my nodes get which partition they are in.
    part_vert = []
    for z in range(Nz):
        for y in range(Ny):
            for x in range(Nx):
                for i in range(len(partList)):
                    inPart = partList[i].inPart(x,y,z)
                    if inPart:
                        part_vert.append(i)
                        break
    return part_vert

@numba.jit(nopython=False,parallel=True,fastmath=True)
def set_geometric_partition_improved(Nx,Ny,Nz,px,py,pz):
    """
    Nx = integer, number of lattice points in the x direction
    Ny = integer, number of lattice points in the y direction
    Nz = integer, number of lattice points in the z direction
    px,py,pz = integers: number of partitions in the x,y and z direction

    write the partition information to parts.lbm
    returns part_vert for a geometric partition
    """
    #size of partitions in x y z direction
    xpartsize = np.zeros(px,dtype=np.int32)
    ypartsize = np.zeros(py,dtype=np.int32)
    zpartsize = np.zeros(pz,dtype=np.int32)
    bx = Nx/px; by = Ny/py; bz = Nz/pz # block x,y and z minimum sizes

    for z in range(pz):
        if ((Nz%pz)>z):
            zpartsize[z] = bz + 1
        else:
            zpartsize[z] = bz

    for y in range(py):
        if ((Ny%py)>y):
            ypartsize[y] = by + 1
        else:
            ypartsize[y] = by

    for x in range(px):
        if ((Nx%px)>x):
            xpartsize[x] = bx + 1
        else:
            xpartsize[x] = bx

    part_vert = np.zeros(shape=Nz*Ny*Nx,dtype=np.int32)
    #part_vert = []
    parts = open('parts.lbm','w')

    #write partitions
    index = 0
    for z in range(pz):
        zpartition = z * py * px
        for k in range(zpartsize[z]):
            for y in range(py):
                ypartition = y * px
                for j in range(ypartsize[y]):
                    for x in range(px):
                        xpartition = x
                        for i in range(xpartsize[x]):
                            partition = xpartition + ypartition + zpartition
                            part_vert[index] = partition
                            #part_vert.append(partition)
                            parts.write('%d \n'% partition)
                            index = index + 1

    parts.close()
    return part_vert

@numba.jit(nopython=True,parallel=True,fastmath=True)
def count_cuts(adj,vert_part):
    """
    adj is an iterable object containing an adjacency matrix for a logically graph-like object
    vert_part is a listing of which partition each graph vertex is in.

    returns cut - integer with the number of edges that cross partitions
    """

    edge_cuts = set()
    adjSize = adj.shape[0]
    ngbsSize = adj.shape[1]

    for i in range(adjSize):
        my_part = vert_part[i]
        ngbs = adj[i]
        for j in range(ngbsSize):
            #ngb_part = vert_part[n]
            n = int(ngbs[j])
            if (vert_part[n]!= my_part):
                min_vert = min(i,n)
                max_vert = max(i,n)
                cut_edge = (min_vert,max_vert)
                edge_cuts.add(cut_edge)
    return len(edge_cuts)

@numba.jit(nopython=False,parallel=True,fastmath=True)
def set_adjacency(Nx,Ny,Nz,ex,ey,ez):
    """
    Nx = num of lattice points in X-direction
    Ny = num of lattice points in Y-direction
    Nz = num of lattice points in Z-direction
    ex = lattice speeds in X-direction
    ey = lattice speeds in Y-direction
    ez = lattice speeds in Z-direction

    returns adjDict = dictionary where the keys are the global lattice point numbers
    and the values are lists of neighboring lattice points
    """
    adjDict = np.zeros(shape=(Nx*Ny*Nz,len(ex)),dtype=np.int64)
    #adjDict = {}

    for z in range(Nz):
        for y in range(Ny):
            for x in range(Nx):
                gid = x + y*Nx+z*Nx*Ny
                for spd in range(len(ex)):
                    dx = int(ex[spd]); dy = int(ey[spd]); dz = int(ez[spd]);
                    tx = (x+dx)%Nx; ty= (y+dy)%Ny; tz = (z+dz)%Nz
                    tid = tx+ty*Nx+tz*Nx*Ny
                    adjDict[gid][spd] = tid
                    #adjDict.setdefault(gid,[]).append(tid)

    return adjDict

if __name__=='__main__':
    Ny_divs = 7
    Re = 100.
    # overall domain dimensions
    Lx_p = 6.4 * 2 # "thickness"
    Ly_p = 3. * 2 # "height"
    #Lz_p = 10. # "length"
    Lz_p = 14. * 2


    # describe brick dimensions and location
    h_brick = 1./2.
    #z_brick = 1./2.
    z_brick=4.
    x_brick = 1./2.

    Lo = h_brick;#characteristic length is block height
    R=0; x_c = 0; z_c = 0;# not used.

    Ny = math.ceil((Ny_divs-1)*(Ly_p/Lo))+1
    Nx = math.ceil((Ny_divs-1)*(Lx_p/Lo))+1
    Nz = math.ceil((Ny_divs-1)*(Lz_p/Lo))+1
    nnodes = Nx*Ny*Nz

    # compute geometric data only once
    x = np.linspace(0.,Lx_p,Nx).astype(np.float32);
    y = np.linspace(0.,Ly_p,Ny).astype(np.float32);
    z = np.linspace(0.,Lz_p,Nz).astype(np.float32);
    numEl = Nx*Ny*Nz
    Y,Z,X = np.meshgrid(y,z,x);

    XX = np.reshape(X,numEl)
    YY = np.reshape(Y,numEl)
    ZZ = np.reshape(Z,numEl)

    u = np.zeros_like(XX)
    v = np.zeros_like(XX)
    w = np.zeros_like(XX)

    """
    Now that the geometry is formed, create the adjacency matrix and partition with
    pymetis.

    Also, get a partition with my geometric partitioner and compare the two partitions
    visually.  Perhaps get a surface to volume ratio for both partitions and from
    that make a prediction as to the performance.
    """

    ex = [0.,1.,-1.,0.,0.,0.,0.,1.,-1.,1.,-1.,1.,-1.,1.,-1.]
    ey = [0.,0.,0.,1.,-1.,0.,0.,1.,1.,-1.,-1.,1.,1.,-1.,-1.]
    ez = [0.,0.,0.,0.,0.,1.,-1.,1.,1.,1.,1.,-1.,-1.,-1.,-1.]

    print 'Total lattice points = %d.'%(Nx*Ny*Nz)

    print 'Setting adjacency list'
    adjDict = set_adjacency(int(Nx),int(Ny),int(Nz),ex,ey,ez)

    N_parts = 8

    print 'Nx = %d ' % Nx
    print 'Ny = %d ' % Ny
    print 'Nz = %d ' % Nz

    print 'getting METIS partition'
    if (NO_PYMETIS==1):
      print "pymetis is not available"
      sys.exit()
    cuts, part_vert = part_graph(N_parts,adjDict)

    print 'getting part_advisor partition'
    px,py,pz = part_advisor(Nx,Ny,Nz,N_parts)

    # make sure all of these things are integers...
    Nx = int(Nx); Ny = int(Ny); Nz = int(Nz)
    px = int(px); py = int(py); pz = int(pz)

    start1 = time.time()
    print "set geometric partition 1"
    part_vert_pa1 =  set_geometric_partition(Nx,Ny,Nz,px,py,pz)
    print "set geometric partition 1 done"
    print time.time()-start1

    print "write geometric partition 1"
    parts = open('parts_test.lbm','w')
    for p in part_vert_pa1:
        parts.write('%d \n'% p)
    parts.close()
    print "write done geometric partition 1"
    print time.time()-start1

    start2 = time.time()
    print "set geometric partition 2"
    part_vert_pa2 = set_geometric_partition_improved(Nx,Ny,Nz,px,py,pz)
    print "set geometric partition 2 done"
    print time.time()-start2

    cuts_pa1 = count_cuts(adjDict,part_vert_pa1)
    cuts_pa2 = count_cuts(adjDict,part_vert_pa2)

    print 'cuts pa1 = %d ' % cuts_pa1
    print 'cuts pa2 = %d ' % cuts_pa2

    # cuts_metis = count_cuts(adjDict,part_vert)
    # cuts_pa = count_cuts(adjDict,part_vert_pa)
    # cuts_1D = count_cuts(adjDict,part_vert1D)
    #
    # print 'cuts metis = %d ' % cuts_metis
    # print 'cuts pa = %d ' % cuts_pa
    # print 'cuts_1D = %d ' % cuts_1D
    # print 'writing partition to VTK file'
    #
    # dims = [Nx,Ny,Nz]
    # origin = [0,0,0]
    # dx = x[1]-x[0]; dy = y[1]-y[0]; dz = z[1]-z[0];
    # spacing = [dx,dy,dz]
    #
    # saveStructuredPointsVTK_ascii(part_vert,'partitions','partition_metis.vtk',dims,origin,spacing)
    # saveStructuredPointsVTK_ascii(part_vert_pa,'parititions','partition_pa.vtk',dims,origin,spacing)
    # saveStructuredPointsVTK_ascii(part_vert1D,'partitions','partition_1D.vtk',dims,origin,spacing)
    #saveVelocityAndPressureVTK_binary(part_vert,u,v,w,x,y,z,'partition.vtk',dims)
    #saveScalarStructuredGridVTK_binary(part_vert,'partition',XX,YY,ZZ,'partition.vtk',[Nx,Ny,Nz])
