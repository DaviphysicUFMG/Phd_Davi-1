import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

x,y,z = np.loadtxt('output_001.dat',unpack=True)
N = int(len(z)**0.5)
z = z.reshape(N,N)

plt.imshow(z,extent=(np.amin(x),np.amax(x),np.amin(y),np.amax(y)),cmap='RdGy',interpolation='spline36')

#plt.imshow(z,extent=(-10.,10.,-10.,10.),cmap=cm.inferno,interpolation='nearest')

#plt.contour(X,Y,Z,colors='black')

plt.clim(0,0.505)
plt.title('T = 15.')
plt.xlabel('$q_x(r.l.u)$')
plt.ylabel('$q_y(r.l.u)$')
plt.colorbar()
plt.show()
