#! /usr/bin/env python2
import sys
sys.path.insert(0, '..')
from afmm import *

def lin(left, right, width):
    return lambda x: x*(right-left)*1.0/width + left

lamb = 975.0
air = material(1.0**2)
substrate = material((1.5)**2, color='yellow')
#non-magnetic absorber
#a1 = material(lin(1.0**2, (1.0+1.0j)**2, lamb/4),color='grey')
#a2 = material(lin((2.9+1.0j)**2, 2.9**2, lamb/4),color='grey')

#PML
#mu = 5*(1.0+1.0j)
#a1 = material(1.0/mu, mu, mu, color='grey')
#a2 = material(2.9**2/mu, 2.9**2*mu, mu, color='grey')

def fr(gamma, qq):
    q = 2*qq
    return lambda x: (1.0-gamma*np.sin(np.pi*x/q)**2)*np.cos(np.pi*x/q)**2

def fl(gamma, qq):
    q = 2*qq
    return lambda x: (1.0-gamma*np.sin(np.pi*(qq-x)/q)**2)*np.cos(np.pi*x/q)**2

gamma = 1.0/(1.0-1.0j)
kk = 1000.0
a1l = pml(air,fl(gamma,kk), 0.5)
a1r = pml(air,fr(gamma,kk), 0.5)
a2l = pml(substrate,fl(gamma,kk), 0.5)
a2r = pml(substrate,fr(gamma,kk), 0.5)

# convention: a+bj for absorbing
#w = 1.0
wwout = 3000.0
a = layer([[a1l, kk], [air, wwout], [air, 300.0], [air, wwout], [a1r, kk]])
b = layer([[a2l, kk], [substrate, wwout], [air, 300.0], [substrate, wwout], [a2r, kk]])
c = layer([[a2l, kk], [substrate, wwout], [substrate, 300.0], [substrate, wwout], [a2r, kk]])

geometry = stack([[a,2000.0],[b,500.0],[c,5000.0]][::-1])
#geometry.paint().show()
o = 50
waveguide = afmm(geometry, order=o, lambda0=lamb)
waveguide.compute()
#lamb_fundTM0=min(waveguide.modes(2),key=lambda x: x.real)
#n_fundTM0=list(waveguide.modes(2)).index(lamb_fundTM0)
#waveguide.inputmode([0]*n_fundTM0+[1.0]+[0]*(2*o-n_fundTM0))
waveguide.inputfunc(lambda x: np.exp(-((x-4150.0)/1000.0)**2))
# plot real input
f = lambda x: np.exp(-((x-4150.0)/1000.0)**2)
frac = 800
x = [i*1.0*waveguide.w/frac for i in range(frac)]
y = [None]*frac
for i in range(frac):
    y[i] = f(x[i])
plt.plot(x,y)
plt.show()

#waveguide.inputfunc(lambda x: 1.0)
waveguide.plotinput()
waveguide.plotreflect()
waveguide.plotW(1)
waveguide.plotHy_xz()

