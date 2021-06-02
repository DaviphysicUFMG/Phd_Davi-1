import numpy as np

def rec(a1x,a1y,a2x,a2y):
    A = 2.0*np.pi/(a1x*a2y - a1y*a2x)
    b1x = A*a2y
    b1y = -A*a2x
    b2x = -A*a1y
    b2y = A*a1x
    return (b1x,b1y,b2x,b2y)

a1x = 2.*np.cos(60.*np.pi/180.)
a1y = 2.*np.sin(60.*np.pi/180.)
a2x = 2.*1.
a2y = 0.

b1x,b1y,b2x,b2y = rec(a1x,a1y,a2x,a2y)

print(b1x,b1y)
print(b2x,b2y)
