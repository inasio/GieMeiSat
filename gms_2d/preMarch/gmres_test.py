import numpy as NY
import scipy as SY
import scipy.sparse as SS
import scipy.sparse.linalg as SL
import scipy.linalg as SLA
import matplotlib.pyplot as PP
import mpl_toolkits.mplot3d.axes3d as P3
import time as TT

def rand_Jacobian_matrix(n):
  
  I = SS.eye(n,n)
  data = NY.random.rand(3,n)
  B = SS.spdiags(data,[0,-1,1],n,n)
  A = SS.kron(B,I) + SS.kron(I,B)
  return A
 
if __name__ == "__main__":
  
  N = 500
  A = rand_Jacobian_matrix(N)
  #A = A + 50*SS.eye(N**2,N**2)
  #NY.random.seed(90)
  B = NY.random.rand(N**2); B = B/NY.linalg.norm(B)
  T3 = TT.clock()
  X = SL.dsolve.spsolve(A,B)
  T4 = TT.clock(); Tspsolve = T4 - T3
  print Tspsolve, NY.linalg.norm(X)
  #epsi = 1e-5
  #X0 = X + NY.random.rand(N**2)*epsi
  #T1 = TT.clock()
  #[X1,flag] = SL.isolve.gmres(A,B,X0,1e-7,restrt=N**2,maxiter=N**2)
  #T2 = TT.clock(); Tgmres = T2 - T1
  #print Tgmres
  #Erro1 = NY.linalg.norm(X1-X)
  #Erro2 = NY.linalg.norm(X2-X)