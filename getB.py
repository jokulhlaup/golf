import numpy as np
import scipy as sp
import fec2.fec2 as fec2
from inspect import isfunction

#Base class for fabric distributions. This takes an argument
#for the ODF density = density(X, theta,y) continuous function over
#R3 and orientation angle (theta, phi) \in ([0,2pi],[0,pi])
#X must be a [n,3] numpy array to handle an arbitrary number of points
#one call. 
#Use class FabricDistArrDen to give array density as 
#2D numpy array density[spatial index, spherical index] or Dolfin/Fenics GenericVector density[finite element index].
#Default is isotropic density.
class FabricDist
   def __init__(self,
                pts=None,
                density = (lambda X,theta,phi: (4*3.14159)**-1
                *np.ones([len(X[:,i]),len(theta)]),
                path='lebedev_data.txt',
                lebedev_data=np.zeros([1,3])
		x0=np.array([1,1,1])):
      #Set defaults
      #Error tolerence for computation of Carlson form Rd.
      self.errtol=1e-4
      self.path=path
      self.x0=x0

      if A=None: 
         ABinit(self, density) 

   def ABinit(self,density):
      if (leb_data in globals()) == False:
         #Read in data
         fl=open(path,'w')
         flstr=fl.readlines()
         fl.close()
         nl=len(flstr)
         leb_data=empty([nl,3])
         for i in range(nl)
            leb_data[i,:] = np.array(float(flstr[i].split()))


      #If density is supplied as function of (X \in R3, theta, phi),
      #user must specify a (nx,3) numpy array of R3 points to evaluate
      #the function on.
      #If a numpy array(nx,nl) of densities on (spatial point, lebedev point)
      #is supplied, there is no interaction between spatial points. 
      if isfunction(density) == True:
         if pts=None: 
	    print 'Error: FabricDist object initialized
	    with density=pyfunction and no spatial 
	    points "X" supplied'
	    return
  
         density = density(pts,leb_data[:,1],leb_data[:,2])
      #Case for density supplied as numpy array[nx,nl]
      else
         if len(density[:,1]) != nl:
	    print 'Error: Dimension 1 of supplied array "density"
	    not equal to dimension 0 of "leb_data"'
	    return

         GetA(self)
         GetB(self)

   def GetA(self,density):
      self.A=fec2.GetA(density,leb_data,m,n)
   
   def GetB(self): 
      #Diagonalize the A(:,3,3) matrix
      #Q <- (m,3,3) rotation matrix
      #W <- (:,3) eigenvalues
      self.Q,self.W=fec2.dsyar(self.A,m)
      cur= lambda X: objective(self,X) 
      
      for i in range(m):
         #Find b11 b22 b33 in diagonalized reference
         #frame.
         Bw=scipy.optimize.root(cur,FabricDist.x0,
	                        method='hybr',jac=None)
         B[i,:]=fec2.RotSym(Bw,self.Q)
      return B

   #Compute Aii - Rd, where Rd is the Carlson
   #symmetric form.
   def objective(self,X):
      Y[0]=self.W[i,0] - FabricDist.RD(self,X)
      #Permute Xi
      Y[1]=self.W[i,1] - FabricDist.RD(self,np.array([X[1],X[0],X[2]))
      Y[2]=self.W[i,2] - FabricDist.RD(self,np.array([X[2],X[0],X[1]))
      return Y

   def RD(self,X):
      return fec2.RD(X,self.errtol)
      




