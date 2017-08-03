#! /usr/bin/env python2
import numpy as np
import scipy as sp
import scipy.linalg
import numpy.linalg as lg
import matplotlib.pyplot as plt

lamb = 800.0
n_d = 1.0
n_m = 0.18+5.13j
eps_m = n_m**2
eps_d = n_d**2
kk = 2*lamb
w = 8.0

# calculating overlap integral
n_sp = np.sqrt(eps_m*eps_d/(eps_m+eps_d))
print n_sp
k = 2*np.pi/lamb
xi_m = np.sqrt(eps_m-n_sp**2)
xi_d = np.sqrt(eps_d-n_sp**2)
if xi_m.imag<0:
    xi_m = -xi_m
if xi_d.imag<0:
    xi_d = -xi_d
z0 = 4*lamb
x0 = kk+(w+0.125)*lamb
def psi_sp_plus(x,z):
    ans = np.exp(1j*k*n_sp*(x-x0))
    if z>z0:
        ans = ans*np.exp(1j*k*xi_m*(z-z0))
        return [ans, ans*xi_m/eps_m, -ans*n_sp/eps_m]
    else:
        ans = ans*np.exp(-1j*k*xi_d*(z-z0))
        return [ans, -ans*xi_d/eps_d, -ans*n_sp/eps_d]

def psi_sp_minus(x,z):
    ans = np.exp(-1j*k*n_sp*(x-x0))
    if z>z0:
        ans = ans*np.exp(1j*k*xi_m*(z-z0))
        return [ans, ans*xi_m/eps_m, ans*n_sp/eps_m]
    else:
        ans = ans*np.exp(-1j*k*xi_d*(z-z0))
        return [ans, -ans*xi_d/eps_d, ans*n_sp/eps_d]

def overlap_integral(dz, a, b, c, d):
    l = len(a);
    ans = 0.0;
    for i in range(l):
        ans += (a[i]*b[i]-c[i]*d[i])*dz
    return ans

XX = np.linspace(kk, (2*w+0.25)*lamb+kk,10)
ZZ = np.linspace(0, 6.125*lamb, 500, endpoint=False)
x = XX[6]
H_y_spp_minus = [psi_sp_minus(x,z)[0] for z in ZZ]
plt.plot(ZZ, H_y_spp_minus)
plt.show()
a = np.array([psi_sp_minus(x,z)[0] for z in ZZ])
b = np.array([psi_sp_minus(x,z)[2] for z in ZZ])
c = np.array([psi_sp_plus(x,z)[0] for z in ZZ])
d = np.array([psi_sp_plus(x,z)[2] for z in ZZ])
print a.shape
e = a*b-c*d
print e.shape
