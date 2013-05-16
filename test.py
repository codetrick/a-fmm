#! /usr/bin/env python2
from afmm import *

def lin(left, right, width):
    return lambda x: x*(right-left)*1.0/width + left

lamb = 975.0
air = material(1.0**2)
core = material(3.5**2, color='black')
substrate = material(2.9**2, color='yellow')
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
kk = 500.0
a1 = pml(air,fr(gamma,kk), 0.5)
a2 = pml(substrate,fl(gamma,kk), 0.5)

# convention: a+bj for absorbing
#w = 1.0
wwout = 700.0
a = layer([[a2, kk], [substrate, wwout], [core, 300.0], [air, wwout], [a1, kk]])
b = layer([[a2, kk], [substrate, wwout], [air, 300.0], [air, wwout], [a1, kk]])
c = layer([[a2, kk], [substrate, wwout], [core, 300.0], [air, wwout], [a1, kk]])
d = layer([[a2, kk], [substrate, wwout], [air, 300.0], [air, wwout], [a1, kk]])
e = layer([[a2, kk], [substrate, wwout], [core, 300.0], [air, wwout], [a1, kk]])

geometry = stack([[a,200.0],[b,150.0],[c,150.0],[d,150.0],[e,200.0]][::-1])
#geometry.paint().show()
o = 50
waveguide = afmm(geometry, order=o, lambda0=lamb)
waveguide.compute()
lamb_fundTM0=min(waveguide.modes(4),key=lambda x: x.real)
n_fundTM0=list(waveguide.modes(4)).index(lamb_fundTM0)
waveguide.inputmode([0]*n_fundTM0+[1.0]+[0]*(2*o-n_fundTM0))
waveguide.plotinput()
waveguide.plotreflect()
waveguide.plotW(1)
waveguide.plotHy_xz()
#waveguide.plotHy_z(651)
#waveguide.plotHy_z(653)
#waveguide.plotHy_z(501)
print abs(waveguide.u[4][n_fundTM0])**2

#wg = waveguide
#print np.dot(wg.W[1], np.add(np.dot(wg.X[1], wg.u[1]), np.dot(lg.inv(wg.X[1]), wg.d[1])))
#print np.dot(wg.W[2], np.add(wg.u[2], wg.d[2]))
#print np.dot(wg.W[0], np.add(np.dot(wg.X[0], wg.u[0]), np.dot(lg.inv(wg.X[0]), wg.d[0])))
#print np.dot(wg.W[1], np.add(wg.u[1], wg.d[1]))
#print np.dot(wg.V[1], np.add(np.dot(wg.X[1], wg.u[1]), np.negative(np.dot(lg.inv(wg.X[1]), wg.d[1]))))
#print np.dot(wg.V[2], np.add(wg.u[2], np.negative(wg.d[2])))
