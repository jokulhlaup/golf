#Contains functions for computing fabric parameter
#Depends on fortran, scipy, f2py, and dolfin/fenics
import numpy as np
import scipy as sp
from dolfin import *
import scipy.linalg as linalg
import opil as opil
import fec2


#Gets the second order orientation tensor for discrete fabrics
#c_axes=c_axes[i,j,k] is third order numpy array. 
#i is spatial index, j is jth observation at point i,
#k=1:3 xyz component of the individual c axis vector.
def get_A(c_axes):
  m=c_axes.shape[0]
  n=c_axes.shape[1]
  
  #find the average outer product
  A[:,1]=tensordot(c_axes[:,:,1],c_axes[:,:,1],(1,1))/n
  A[:,2]=tensordot(c_axes[:,:,2],c_axes[:,:,2],(1,1))/n
  A[:,3]=tensordot(c_axes[:,:,3],c_axes[:,:,3],(1,1))/n
  A[:,4]=tensordot(c_axes[:,:,1],c_axes[:,:,2],(1,1))/n
  A[:,5]=tensordot(c_axes[:,:,1],c_axes[:,:,3],(1,1))/n
  A[:,6]=tensordot(c_axes[:,:,2],c_axes[:,:,1],(1,1))/n


  return

#Alternatively, use kernel smoothing to interpolate at different points.
#In this case, we have need the spatial coordinates of each observation.
#c_axes is a (n,3) array of c-axis angles, and X is (n,3) array of c-axis
#positions, and Y is (n,3) array of positions for the kernel density
# to be evaluated. h is smoothing parameter
#Confusingly, we are doing "second order" kernel smoothing:
#Estimating a pdf for each point, and then interpolating that pdf between points using a gaussian kernel.
#This is implemented in fec.f90


def pars(angles,evs):
  m = angles.shape[0]

  if m != 3:
    print('Function pars(angles,evs) requires nx3 arrays as arguments')
  C=opil.golf.cmat(evs,angles,m)

  return C
 
#Takes the effective strain
def eff_strain(u):
    epsxx = u[0].dx(0)
    epszz = u[1].dx(1)
    epsxz = 0.5*(u[0].dx(1) + u[1].dx(0))
    eps_dot = sqrt(0.5*(epsxx**2 + epszz**2 + 2*epsxz**2))
 
   
    epsxx = u[0].dx(0)
    epsyy = u[1].dx(1)
    epszz = u[2].dx(2)
    epsxy = 0.5*(u[0].dx(1) + u[1].dx(0))
    epsxz = 0.5*(u[0].dx(2) + u[2].dx(0))
    epsyz = 0.5*(u[1].dx(2) + u[2].dx(1))
    eps_dot = sqrt(0.5*(epsxx**2 + epsyy**2 + epszz**2 + 2*epsxy**2 + 2*epsxz**2 + 2*epsyz**2))
 
    return np.array(eps_dot)

def ReadVa(NetaI):
   
   etaI=opil.golf.ReadVa(NetaI)
   return etaI
