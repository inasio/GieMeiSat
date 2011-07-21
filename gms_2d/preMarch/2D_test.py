import numpy as NY
import scipy.sparse as SS
import scipy.sparse.linalg as SL
import scipy.linalg as SLA
import matplotlib.pyplot as PP
import mpl_toolkits.mplot3d.axes3d as P3

def Jacobian_matrix(n):
  
  #dxdy = 1./(n-1)**2
  I = SS.eye(n,n)
  data = NY.ones((3,n))
  data[0,:] = -2.
  data[1,:] = 1.
  data[2,:] = 1.
  B = SS.spdiags(data,[0,-1,1],n,n)
  A = SS.kron(B,I) + SS.kron(I,B)
  ###four corners
  A[0,1] = 2.; A[0,n] = 2.
  A[n-1,n-2] = 2.; A[n-1,2*n-1] = 2.
  A[n*(n-1),n*(n-2)] = 2.; A[n*(n-1),n*(n-1)+1] = 2.
  A[n**2-1,n**2-2] = 2.; A[n**2-1,n*(n-1)-1] = 2.
  ###boundaries
  #for i in NY.arange(n-2):
	#A[i+1,i+1] = 3.
	#A[n*(n-1)+1+i,n*(n-1)+1+i] = 3.
	#A[n*(i+1),n*(i+1)] = 3.
	#A[n*(i+2)-1,n*(i+2)-1] = 3.
  return A
  
def initial_condition(n):
  b = NY.zeros([n,n])
  for i in NY.arange(n):
	for j in NY.arange(n):
	  #b[i,j] = NY.cos(NY.pi*(i+1)/(n+1))*NY.cos(NY.pi*(j+1)/(n+1))
	  b[i,j] = NY.random.rand()
  return b.flatten()
  
def plotty(xx,yy,dx,n):
  z1 = xx.reshape(n,n)
  z2 = yy.reshape(n,n)
  x = NY.arange(dx,1.-dx**2,dx)
  y = NY.arange(dx,1.-dx**2,dx)
  xc, yc = NY.meshgrid(x,y)
  fig = PP.figure(1)
  ax = P3.Axes3D(fig)
  ax.contourf3D(xc,yc,z1,60)
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  fig = PP.figure(2)
  ax = P3.Axes3D(fig)
  ax.contourf3D(xc,yc,z2,60)
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')
  PP.show()

if __name__ == "__main__":

  RelError = NY.zeros(6)
  Numpoints = NY.zeros(6)
  #for i in NY.arange(6):
	
	#N = 4*2**i 
  N = 4
  Lx = 1.; Ly = 1.
  Dx = Lx/(N-1); Dxdy = Dx**2; Dt = 1

  J = Jacobian_matrix(N)
  print J.todense()
  