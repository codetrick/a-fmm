#! /usr/bin/env python2
import scipy.integrate as inte
from afmm import *
from cmath import *

def lin(left, right, width):
    return lambda x: x*(right-left)*1.0/width + left

lamb = 800.0
n_d = 1.0
n_m = 0.18+5.13j
eps_m = n_m**2
eps_d = n_d**2

air = material(eps_d)
substrate = material(eps_m, color='yellow')
#non-magnetic absorber
#a1 = material(lin(1.0**2, (1.0+1.0j)**2, lamb/4),color='grey')
#a2 = material(lin((2.9+1.0j)**2, 2.9**2, lamb/4),color='grey')

#PML
#mu = 5*(1.0+1.0j)
#a1 = material(1.0/mu, mu, mu, color='grey')
#a2 = material(2.9**2/mu, 2.9**2*mu, mu, color='grey')

def fr(gamma, q):
    return lambda x: (1.0-gamma*np.sin(np.pi*x/q)**2)*np.cos(np.pi*x/q)**2

def fl(gamma, q):
    return lambda x: (1.0-gamma*np.sin(np.pi*(q-x)/q)**2)*np.cos(np.pi*x/q)**2

gamma = 1.0/(1.0-1.0j)
kk = 2*lamb
a1l = pml(air,fl(gamma,kk), 0.5)
a1r = pml(air,fr(gamma,kk), 0.5)
a2l = pml(substrate,fl(gamma,kk), 0.5)
a2r = pml(substrate,fr(gamma,kk), 0.5)

# convention: a+bj for absorbing
w = 10.0
wasym = 0.0
a = layer([[a1l, kk], [air, (2*w+0.25)*lamb], [a1r, kk]])
b = layer([[a2l, kk], [substrate, w*lamb], [air, (0.25+wasym)*lamb], [substrate, (w-wasym)*lamb], [a2r, kk]])
c = layer([[a2l, kk], [substrate, w*lamb], [air, 0.25*lamb], [substrate, w*lamb], [a2r, kk]])

geometry = stack([[a,8*lamb],[b,0.125*lamb],[c,2*lamb]])
#geometry.paint().show()
o = 150
waveguide = afmm(geometry, order=o, lambda0=lamb)
waveguide.compute()
lamb_fundTM0=min(waveguide.modes(2),key=lambda x: x.real)
n_fundTM0=list(waveguide.modes(2)).index(lamb_fundTM0)
waveguide.inputmode([0]*n_fundTM0+[1.0]+[0]*(2*o-n_fundTM0))

waveguide.plotinput()
waveguide.plotreflect()
waveguide.plotW(1)
waveguide.plotHy_xz()

# calculating overlap integral
n_sp = np.sqrt(eps_m*eps_d/(eps_m+eps_d))
k = 2*np.pi/lamb
xi_m = np.sqrt(eps_m-n_sp**2)
xi_d = np.sqrt(eps_d-n_sp**2)
if xi_m.imag<0:
    xi_m = -xi_m
if xi_d.imag<0:
    xi_d = -xi_d
z0 = 8*lamb
x0 = kk+(w+0.125)*lamb
def psi_sp_plus(x,z):
    ans = np.exp(1j*k*n_sp*(x-x0))
    if z>z0:
        ans = ans*np.exp(1j*k*xi_m*(z-z0))
        return [ans, -ans*xi_m/eps_m, ans*n_sp/eps_m]
    else:
        ans = ans*np.exp(-1j*k*xi_d*(z-z0))
        return [ans, ans*xi_d/eps_d, ans*n_sp/eps_d]

def psi_sp_minus(x,z):
    ans = np.exp(-1j*k*n_sp*(x-x0))
    if z>z0:
        ans = ans*np.exp(1j*k*xi_m*(z-z0))
        return [ans, -ans*xi_m/eps_m, -ans*n_sp/eps_m]
    else:
        ans = ans*np.exp(-1j*k*xi_d*(z-z0))
        return [ans, ans*xi_d/eps_d, -ans*n_sp/eps_d]

def overlap_integral(zz, a, b, c, d):
    e = a*b - c*d
    return inte.simps(e,zz)

XX = np.linspace(1.2*kk, (2*w+0.25)*lamb+kk*0.8,10)
ZZ = np.linspace(0*lamb, 10.125*lamb, 1000, endpoint=False)
alpha_plus = [None]*len(XX)
alpha_minus = [None]*len(XX)

for i in range(len(XX)):
    x = XX[i]
    E_z = np.array([waveguide.Ez(x,z) for z in ZZ])
    H_y = np.array([waveguide.Hy(x,z) for z in ZZ])
    H_y_spp_minus = np.array([psi_sp_minus(x,z)[0] for z in ZZ])
    E_z_spp_minus = np.array([psi_sp_minus(x,z)[2] for z in ZZ])
    H_y_spp_plus = np.array([psi_sp_plus(x,z)[0] for z in ZZ])
    E_z_spp_plus = np.array([psi_sp_plus(x,z)[2] for z in ZZ])
    N_sp = 1j/k*(n_d**4-n_m**4)/(n_d**3*n_m**3)
    alpha_plus[i] = overlap_integral(ZZ,H_y,E_z_spp_minus,E_z,H_y_spp_minus)/N_sp
    alpha_minus[i] = overlap_integral(ZZ,E_z,H_y_spp_plus,H_y,E_z_spp_plus)/N_sp

plt.plot(XX,alpha_plus,XX,alpha_minus)
plt.show()
plt.plot(XX,[180*phase(x)/np.pi for x in alpha_plus],XX,[180*phase(x)/np.pi for x in alpha_minus])
plt.show()
print [180*phase(x)/np.pi for x in alpha_plus]

