#How it works::
#Takes inputs of 3 fabric eigenvalues, as well as material symmetry Euler angles
#(ie, from fabric eigenvectors)
#etaI is the other input, which gives viscosities as a function of eigenvalues
#based on an arbitrary micro-macro model (.va file, normally)

#Be sure to have it in the same directory.
#This imports the other stuff also
import fab_utils as fab_utils
import numpy as np
import scipy as sp
from dolfin import *
import  opil as opil
#Define function spaces
#ProblemMesh=UnitCube(10,10,10)
###########################################################
qs_length = 40.0        # side length of ice column, used in BC definition
qs_height = 120.0       # height of ice column to be used in BC definition           

# Load the mesh from a gmsh generated file
ProblemMesh = UnitCube(20,20,20)#Mesh("2cyl.xml.gz")



############################################################

#######################
V = VectorFunctionSpace(ProblemMesh, "DG", 2)# + VectorFunctionSpace(ProblemMesh, "B", 4)
Q = FunctionSpace(ProblemMesh, "CG",1)


Y = V * Q

tol = 1.E-6

# Boundaries
def right(x, on_boundary): return x[0] > (1.0 - DOLFIN_EPS)
def left(x, on_boundary): return x[0] < DOLFIN_EPS
def top_bottom(x, on_boundary):
    return x[1] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS

# No-slip boundary condition for velocity
noslip = Constant((0.0, 0.0, 0.0))
bc0 = DirichletBC(Y.sub(0), noslip, top_bottom)

sc=Constant(1/(3600*24*365))
# Inflow boundary condition for velocity
inflow = Expression(("sin(x[1]*pi)", "0.0", "0.0"))
bc1 = DirichletBC(Y.sub(0), inflow, right)

# Boundary condition for pressure at outflow
zero = Constant(0)
bc2 = DirichletBC(Y.sub(1), zero, left)

# Collect boundary conditions
bc = [bc0, bc1, bc2]






#Declare C




C00=Function(Q)
C01=Function(Q)
C02=Function(Q)
C03=Function(Q)
C04=Function(Q)
C05=Function(Q)

C10=Function(Q)
C11=Function(Q)
C12=Function(Q)
C13=Function(Q)
C14=Function(Q)
C15=Function(Q)

C20=Function(Q)
C21=Function(Q)
C22=Function(Q)
C23=Function(Q)
C24=Function(Q)
C25=Function(Q)

C30=Function(Q)
C31=Function(Q)
C32=Function(Q)
C33=Function(Q)
C34=Function(Q)
C35=Function(Q)

C40=Function(Q)
C41=Function(Q)
C42=Function(Q)
C43=Function(Q)
C44=Function(Q)
C45=Function(Q)

C50=Function(Q)
C51=Function(Q)
C52=Function(Q)
C53=Function(Q)
C54=Function(Q)
C55=Function(Q)











#Get vertex (DG = more) count
m=C00.vector()
m=m.size()
#
#
##############################################
##Use icetools as a control##################
##############################################
##################################################################################
## Define the body force f, i.e. gravity. Here I use a rotated coordinate
## system similar to the Cuffey & Paterson p. 295 Fig. 8.5a, just in 3D 
#alpha = 10.0            # inclination of glacier
#alphar = alpha*pi/180
#g = -9.81               # gravitational constant
#rho = 917.0             # ice density
#f_x0 = sin(alphar)*g*rho
#f_x2 = cos(alphar)*g*rho
### A note here on how Fenics defines coordinates
### x[0] == x
### x[1] == y
### x[2] == z
#f = Constant((f_x0, 0, f_x2))
#
#Aglen = 2.4e-24         # Glen flow parameter for temperate ice (Cuffey & Paterson,2010 p. 73)
#nglen = 3.0             # Glen's n
#nu = Constant(8e13)     # initial viscosity of ice, before viscosity iteration
#
## Load the mesh from a gmsh generated file
##ProblemMesh = Mesh("column3D.xml.gz")
## or select Fenics generated mesh. To use this mesh, uncomment the line below.
## Note that nx, ny, nz define the number of mesh points in each dimension.
## nx = 4
## ny = nx
## nz = 3*nx
##ProblemMesh = Box(0, 0, 0, qs_length, qs_length, qs_height, nx, ny, nz)
#
## Define the sub domains for the boundary conditions
#def NoslipBoundary(x, on_boundary):
#    return x[2] < DOLFIN_EPS and on_boundary
#
## Define the periodic boundary on the vertical faces in X direction
#class PeriodicBoundary_x(SubDomain):
#
#    def inside(self, x, on_boundary):
#        return x[0] == 0 and on_boundary
#
#    def map(self, x, y):
#        y[0] = x[0] - cs_length
#        y[1] = x[1]
#        y[2] = x[2]
## Define the periodic boundary on the vertical faces in Y direction
#class PeriodicBoundary_y(SubDomain):
#
#    def inside(self, x, on_boundary):
#        return x[1] == 0 and on_boundary
#
#    def map(self, x, y):
#        y[1] = x[1] - cs_length
#        y[0] = x[0]
#        y[2] = x[2]   

# Apply a no-slip boundary condition for velocity
#noslip = Constant((0,0,0))
#bc0 = DirichletBC(Y.sub(0), noslip, NoslipBoundary)
## Apply the periodic boundary condition in X
#pbc_x = PeriodicBoundary_x()
#bc1 = PeriodicBC(Y.sub(0), pbc_x)
## Apply the periodic boundary condition in Y
#pbc_y = PeriodicBoundary_y()
#bc2 = PeriodicBC(Y.sub(0), pbc_y)

# Collect boundary conditions
#bc = [bc0, bc1, bc2]



###############################################################
###############################################################
#End###############################

(u, p) = TrialFunctions(Y)
(v, q) = TestFunctions(Y)
#####################
nu=Constant(8e4)
nu_lim= 1e25 #Pa/s #1e15
f=Constant(("0.0","0.0","0.0"))

W=np.empty([m,6])
W[:,0]=3.5*10**-25 #Pa^3 s
W[:,1]=3.0
W[:,2]=6.4*10**4 #J/mol
W[:,3]=11.5*10**4 #J/mol
W[:,4]=268 #K
W[:,5]=263 #263+2*10**-8 *Pressure

RelError=1

RelTol=0.02/(365*24*3600)

MaxViscIter=10
nViscIter=0
nu=Constant(8e14)
U_sol=np.array([])


Angle=np.zeros((m,3))

ai=0.33333*np.ones((m,3))
#ai=np.zeros((m,3))
#ai[:,1]=1*np.ones(m)


i=Index()
j=Index()

#Define forms


Up=Function(Y)
U=Function(Y)
U_l=Function(V)
#Define the "C" matrix (Gillet-Chaulet 2006)
C=np.array(fab_utils.pars(Angle,ai),ndmin=3)
tf=Function(Q)
tf.vector()[:]=C[:,2,2]
C00.vector()[:]=C[:,0,0]
C10.vector()[:]=C[:,1,0]
C20.vector()[:]=C[:,2,0]
C30.vector()[:]=C[:,3,0]
C40.vector()[:]=C[:,4,0]
C50.vector()[:]=C[:,5,0]

C01.vector()[:]=C[:,0,1]
C11.vector()[:]=C[:,1,1]
C21.vector()[:]=C[:,2,1]
C31.vector()[:]=C[:,3,1]
C41.vector()[:]=C[:,4,1]
C51.vector()[:]=C[:,5,1]

C12.vector()[:]=C[:,1,2]
C22.vector()[:]=C[:,2,2]



C32.vector()[:]=C[:,3,2]
C42.vector()[:]=C[:,4,2]
C52.vector()[:]=C[:,5,2]
C02.vector()[:]=C[:,0,2]

C13.vector()[:]=C[:,1,3]
C23.vector()[:]=C[:,2,3]
C33.vector()[:]=C[:,3,3]
C43.vector()[:]=C[:,4,3]
C53.vector()[:]=C[:,5,3]
C03.vector()[:]=C[:,0,3]

C14.vector()[:]=C[:,1,4]
C24.vector()[:]=C[:,2,4]
C34.vector()[:]=C[:,3,4]
C44.vector()[:]=C[:,4,4]
C54.vector()[:]=C[:,5,4]
C04.vector()[:]=C[:,0,4]

C15.vector()[:]=C[:,1,5]
C25.vector()[:]=C[:,2,5]
C35.vector()[:]=C[:,3,5]
C45.vector()[:]=C[:,4,5]
C55.vector()[:]=C[:,5,5]
C05.vector()[:]=C[:,0,5]

print np.array(C22.vector()[:]).shape
print np.array(C22.vector()[:])


def nrvisc(u,W):
    
    epsxx = (u[0].dx(0))
    epsyy = (u[1].dx(1))
    epszz = (u[2].dx(2))
    epsxy = (0.5*(u[0].dx(1) + u[1].dx(0)))
    epsxz = (0.5*(u[0].dx(2) + u[2].dx(0)))
    epsyz = (0.5*(u[1].dx(2) + u[2].dx(1)))
    eps0 = (0.5*(epsxx**2 + epsyy**2 + epszz**2 + 2*epsxy**2 + 2*epsxz**2 + 2*epsyz**2))
 
    #nu=opil.golf.nrvisc(W,eps_e) 
    #Converting function vals np array gets error with f2py. Do in Python:
    R=8.131
    Q=((W[:,5] >= W[:,4])*W[:,2] + (W[:,5] > W[:,4])*W[:,3]) 
    W=W 
    nu=0.5*W[:,0]*np.exp(-Q/R*(1/W[:,4]-1/W[:,5]))*pow(eps0,((1-W[:,1])/W[:,1]))

    print nu.shape

    return nu

class GolfFlow(NonlinearVariationalProblem):
   def __init__(self,a,L,bcs):
      self.L=L
      self.a=A
      self.reset_sparsity=True
   def F(self,b,x):
      assemble(self.L,tensor=b)
   def J(self,A,x):
      assemble(self.a,tensor=A,reset_sparsity=self.reset_sparsity)
      self.reset_sparsity=False
      #Sparsity should be the same. Don't redo it after the first time.


nu=1
eps=10e-6
#Iterate for viscosity
while (RelError>RelTol) & (nViscIter<=MaxViscIter):

   #Copy the old u to u_l

#a=   nu*(u[i]*v[i])*dx #test
#This takes a long time to assemble.
   a=    (v[0].dx(0)*(C00*u[0].dx(0)+C01*u[1].dx(1)+C02*u[2].dx(2)+C03*u[1].dx(2)+C04*u[2].dx(0)+C05*u[0].dx(1))   \
       +v[2].dx(2)*(C20*u[0].dx(0)+C21*u[1].dx(1)+C22*u[2].dx(2)+C23*u[1].dx(2)+C24*u[2].dx(0)+C25*u[0].dx(1))   \
       +(v[1].dx(0)+v[0].dx(1))*(C50*u[0].dx(0)+C51*u[1].dx(1)+C52*u[2].dx(2)+C53*u[1].dx(2)+C54*u[2].dx(0)+C55*u[0].dx(1)) \
       +(v[2].dx(0)+v[0].dx(2))*(C40*u[0].dx(0)+C41*u[1].dx(1)+C42*u[2].dx(2)+C43*u[1].dx(2)+C44*u[2].dx(0)+C45*u[0].dx(1)) \
       +v[1].dx(1)*(C10*u[0].dx(0)+C11*u[1].dx(1)+C12*u[2].dx(2)+C13*u[1].dx(2)+C14*u[2].dx(0)+C15*u[0].dx(1))   \
       +(v[2].dx(1)+v[1].dx(2))*(C30*u[0].dx(0)+C31*u[1].dx(1)+C32*u[2].dx(2)+C33*u[1].dx(2)+C34*u[2].dx(0)+C35*u[0].dx(1)))*dx \
       + v[i].dx(i)*p*dx + q*u[i].dx(i)*dx + 1e-9*p*q*dx 
   C00f = File("C22.pvd")
   C00f = C22
   #This is the isotropic Stokes flow case
#   a=C22*(u[j].dx(i)*v[j].dx(i))*dx +  v[i].dx(i)*p*dx + q*u[i].dx(i)*dx # 10e-16*p*q*dx
   L=f[i]*v[i]*dx
   
   #Assemble the preconditioner matrix
   #Might want to figure out how to do this for the anisotropic case.
   #Can always use standard preconditioner
#   bs=nu*u[i].dx(j)*v[i].dx(j)*dx + p*q*dx
#   bs=  nu*(v[0].dx(0)*(C00*u[0].dx(0)+C01*u[1].dx(1)+C02*u[2].dx(2)+C03*u[1].dx(2)+C04*u[2].dx(0)+C05*u[0].dx(1))   \
#   +v[2].dx(2)*(C20*u[0].dx(0)+C21*u[1].dx(1)+C22*u[2].dx(2)+C23*u[1].dx(2)+C24*u[2].dx(0)+C25*u[0].dx(1))   \
#   +2*v[0].dx(1)*(C50*u[0].dx(0)+C51*u[1].dx(1)+C52*u[2].dx(2)+C53*u[1].dx(2)+C54*u[2].dx(0)+C55*u[0].dx(1)) \
#   +2*v[0].dx(2)*(C40*u[0].dx(0)+C41*u[1].dx(1)+C42*u[2].dx(2)+C43*u[1].dx(2)+C44*u[2].dx(0)+C45*u[0].dx(1)) \
#   +v[1].dx(1)*(C10*u[0].dx(0)+C11*u[1].dx(1)+C12*u[2].dx(2)+C13*u[1].dx(2)+C14*u[2].dx(0)+C15*u[0].dx(1))   \
#   +2*v[1].dx(2)*(C30*u[0].dx(0)+C31*u[1].dx(1)+C32*u[2].dx(2)+C33*u[1].dx(2)+C34*u[2].dx(0)+C35*u[0].dx(1)))*dx

   #Assemble main system
   A,B=assemble_system(a,L,bc)

   #Assemble preconditioner
#   P,b=assemble_system(bs,L,bc)
   

   nViscIter=nViscIter+1
   print "Viscosity iteration:", nViscIter
#   FlowProblem = LinearVariationalProblem(a,L,Up, bcs=bc)
   solve(A,Up.vector(),B,"tfqmr","amg")
#   isolver.set_operators(A,P)
#   solver.solve(Up.vector(),B) 


#   solve(a == L, Up, bcs=bc)

   #Count
#   FlowProblem = LinearVariationalProblem(a,L,Up, bcs=bc)
#   solver=LinearVariationalSolver(FlowProblem)
#   solver.parameters["linear_solver"] =  "cholesky"
#   solver.parameters["preconditioner"] = "ilu"
#   solver.solve()
   
   (U,P)=Up.split(deepcopy=True)
   #Compute the new relative error
   if nViscIter > 1:
      RelError = errornorm(U,U_l,mesh=ProblemMesh) 
      
   
   if RelError > RelTol:
      #Now iterate for viscosity
      #First get the effective strain rate
      #V_v=FunctionSpace(ProblemMesh,"CG",3)
      nut=TrialFunction(Q)
      r=TestFunction(Q)
      nu = Function(Q)
      nu=nrvisc(U,W)

#      av=inner(nut,r)*dx
#      bv=inner(nrvisc(U,W),r)*dx
#      AV=assemble(av)
#      BV=assemble(bv)
#      solve(AV,nu.vector(),BV)
      
      #solve(av==bv,nu)
      #assemble_system(AV,BV,av,bv)
      #solve(AV,nu.vector(),BV,"ilu")
      
#nu_out=nu.vector < nu_lim
      #nu.vector[nu_out]=nu_lim
#      print 'max',  np.max(nu.vector)
#      print 'min',  np.min(nu.vector)
      #nu0=fab_utils.nrvisc(U,W)
       


      #eta_e=assemble(F)

      #ViscosityProblem = LinearVariationalProblem(a_v,L_v,nu)
      #solver=LinearVariationalSolver(ViscosityProblem)
      #solver.parameters["linear_solver"] = "gmres"
      #solver.parameters["preconditioner"] = "ilu"
      #solver.solve

      #Compute the nr viscosity
      # Calculate the non-linear viscosity as a variational form  
#      a_s = inner(w_s,v_s)*dx
#      L_s = inner(nu_ice(u), v_s)*dx
#      nu = Function(V_s)
      # Solve for viscosity
#      solve(a_s==L_s,nu)      
   #nViscIter=50 
   #Append the new velocity
   U_sol=np.append(U_sol,U)
   #Set the U to U_l
   U_l=U
#   nuf=File("nu.pvd")
#   nuf << nu
#Save U to file
Uend = File("vel_fin.pvd")
Uend << U
Pend = File("p_end.pvd")
Pend << P
