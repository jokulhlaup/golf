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
V = VectorFunctionSpace(ProblemMesh, "CG", 2) #+ VectorFunctionSpace(ProblemMesh, "B", 4)
Q = FunctionSpace(ProblemMesh, "CG",1)
parameters.form_compiler.quadrature_degree=2
Y = V * Q

tol = 1.E-6

# Boundaries
def right(x, on_boundary): return x[0] > (1.0 - DOLFIN_EPS) and on_boundary
def left(x, on_boundary): return x[0] < DOLFIN_EPS and on_boundary
def top_bottom(x, on_boundary):
    return x[1] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS and on_boundary
# No-slip boundary condition for velocity
noslip = Constant((0.0, 0.0, 0.0))
bc0 = DirichletBC(Y.sub(0), noslip, top_bottom)



c=Constant(1/(3600*24*365))
# Inflow boundary condition for velocity
inflow = Expression(("-sin(x[1]*pi)", "0.0", "0.0"))
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
U=Function(Y)
u,p=split(U)
#For Picard
(u_p,p_p)=TrialFunctions(Y)
(v, q) = TestFunctions(Y)

#Define test/trial functions for viscosity
#
R=FunctionSpace(ProblemMesh, "DG",1)
v_v=TestFunction(R)
w_v=Trialfunction(R)
#####################
nu=Constant(8e12)
nu_lim= 1e25 #Pa/s #1e15
f=Constant(("0.0","0.0","0.0"))


W0=Constant(3.5*10**-25) #Pa^3 s
W1=Constant(3.0)
W2=Constant(6.4*10**4) #J/mol
W3=Constant(11.5*10**4) #J/mol
W4=Constant(268) #K
W5=Constant(263) #263+2*10**-8 *Pressure

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

class GolfFlow(NonlinearProblem):
   def __init__(self,a,L,bcs):
      NonlinearProblem.__init__(self,a,L)
      self.L=L
      self.a=a
      self.reset_sparsity=True
   def F(self,b,x):
      assemble(self.L,tensor=b)
   def J(self,A,x):
      assemble(self.a,tensor=A,reset_sparsity=self.reset_sparsity)
      self.reset_sparsity=False
      #Sparsity should be the same. Don't redo it after the first time.

    

#############################
##Define nonlinear viscosity#
#############################
nlvisc(u,W0,W1,W2,W3,W4,W5):
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
   Qe=conditional(gt(W4,W5),W5*W2,W4*W2) + conditional(gt(W5,W4),W5*W3,W3*W4)
   return 0.5*W0*exp(-Qe/R*(1/W4-1/W5))*(eps0**((1-W1)/W1))
nrvisc = lambda u: nlvisc(u,W0,W1,W2,W3,W4,W5)
###############
##Define forms
##############
eps=10e-6
f=Constant((0,0,0))
#This takes a long time to assemble.
F_p=    nu_p*(v[0].dx(0)*(C00*u_p[0].dx(0)+C01*u_p[1].dx(1)+C02*u_p[2].dx(2)+C03*u_p[1].dx(2)+C04*u_p[2].dx(0)+C05*u_p[0].dx(1))   \
>>>>>>> 04712097b920beb71cfc18eaaf11c9449a222670
    +v[2].dx(2)*(C20*u_p[0].dx(0)+C21*u_p[1].dx(1)+C22*u_p[2].dx(2)+C23*u_p[1].dx(2)+C24*u_p[2].dx(0)+C25*u_p[0].dx(1))   \
    +(v[1].dx(0)+v[0].dx(1))*(C50*u_p[0].dx(0)+C51*u_p[1].dx(1)+C52*u_p[2].dx(2)+C53*u_p[1].dx(2)+C54*u_p[2].dx(0)+C55*u_p[0].dx(1)) \
    +(v[2].dx(0)+v[0].dx(2))*(C40*u_p[0].dx(0)+C41*u_p[1].dx(1)+C42*u_p[2].dx(2)+C43*u_p[1].dx(2)+C44*u_p[2].dx(0)+C45*u_p[0].dx(1)) \
    +v[1].dx(1)*(C10*u_p[0].dx(0)+C11*u_p[1].dx(1)+C12*u_p[2].dx(2)+C13*u_p[1].dx(2)+C14*u_p[2].dx(0)+C15*u_p[0].dx(1))   \
    +(v[2].dx(1)+v[1].dx(2))*(C30*u_p[0].dx(0)+C31*u_p[1].dx(1)+C32*u_p[2].dx(2)+C33*u_p[1].dx(2)+C34*u_p[2].dx(0)+C35*u_p[0].dx(1)))*dx \
    + v[i].dx(i)*p*dx + q*u_p[i].dx(i)*dx - f[i]*v[i]*dx

l_p=f[i]*v[i]

#This is the isotropic Stokes flow case
#F=nu*(u[j].dx(i)*v[j].dx(i))*dx +  v[i].dx(i)*p*dx + q*u[i].dx(i)*dx -f[i]*v[i]*dx+ 10e-16*p*q*dx

#F=u[i].dx(j)*v[i].dx(j)*dx + v[i].dx(i)*p*dx + q*u[i].dx(i)*dx -f[i]*v[i]*dx 
#Compute Jacobian


class PicardSolver
   #init method
   #F -> bilinear form
   #l -> linear form
   #u -> Function() to solve for.
   #nu0 -> initial viscosity
   def __init__(self,F,lf=None,u,nu0=Constant(1.0),max_iter=10,max_u_err=0.005):
      self.F=F
      self.l=l
      self.u=u
      self.nu=nu0
      self.max_iter=max_iter
      self.max_u_err=max_u_err
      if lf=None: self.lf=Constant((0.0,0.0,0.0))*self.v
      self.solver_parameters

   def solve(self)

      ############################
      ####Picard Iteration########
      ############################
      #Initialize nu.
      nu=Constant(1)
      while  i < maxiter
         U_k=Uv
         #Define bilinear and linear form. Use viscosity nu. Is constant if first iteration.
         #If not first iteration, is found by solving this variational form for viscoisty.
         replace(L,nu:nu)
         
      
         solve(F==lf,U,bc)
         Uv,P=self.U.split()
         
         #Compute the error. Break if good.
         if i>0: u_err = errornorm(Uk,Uv,"L2")
         if u_err < max_u_err: break
      
         #Compute new viscosity
         F_v=inner(w_v,v_v)*dx
         l_v=inner(nrvisc(Uv),v_v)*dx
         if i=0: nu=Function(R)
         solve(F_v==l_v, nu)
         ################################
         ###Back
         ########################################

   



  
   

J=derivative(F,U)
golfproblem=NonlinearVariationalProblem(F,U,bcs=bc,J=J)
solver=AdaptiveNonlinearVariationalSolver(golfproblem)
#solver.parameters["linear_solver"]="lu"
#solver.parameters["preconditioner"]="ilu"
solver.solve()
#solver.solve(Up.vector(),B) 


#solve(a == 0, U, bc,solver_parameters={'newton_solver':
#                                          {'relative_tolerance':1e-6}})

#Compute the new relative error
if nViscIter > 1:
   RelError = errornorm(U,U_l,mesh=ProblemMesh) 
   

#Set the U to U_l
#   nuf=File("nu.pvd")
#   nuf << nu
#Save U to file
u,p=U.split(deepcopy=True)
Uend = File("vel_fin.pvd")
Uend << u
Pend = File("p_end.pvd")
Pend << p
