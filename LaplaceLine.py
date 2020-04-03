# LaplaceLine.py: Solve Laplace’s eqn, 3D matplot, close shell to quit

from numpy import *
import matplotlib.pylab as p;
from mpl_toolkits.mplot3d import Axes3D

#p.xkcd()                  # for fun, uncomment for xkcd atyle plots
print("Initializing")
Nmax = 100;               # grid points Nx, Ny 
Niter = 70;               # sweeps/iterations to perform
V = zeros((Nmax, Nmax), float)  

print( "Working hard, wait for the figure while I count to 60")
for k in range(0, Nmax-1): V[k,0] = 100.0  # line at 100V

for iter in range(Niter):  # iterations over algorithm
    if iter%10 == 0: print(iter)
    for i in range(1, Nmax-2):
        for j in range(1,Nmax-2): V[i,j] = 0.25*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1])
x = range(0, Nmax-1, 2); y = range(0, 50, 2)    # plot every other point
X, Y = p.meshgrid(x,y)

def functz(V):    # Function returns V(x, y)
    z=V[X,Y]
    return z

Z = functz(V)
fig = p.figure()  # Create figure
ax = Axes3D(fig)  # plot axes
ax.plot_wireframe(X, Y, Z, color = 'r')  # red wireframe
ax.set_xlabel('X')                       # label axes
ax.set_ylabel('Y')
ax.set_zlabel('Potential')
print("Control-C to exit")
p.show()                                 # display fig
