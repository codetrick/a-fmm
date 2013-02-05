from afmm import *

air = material(1.0)
water = material(1.2-10J, color='yellow')
# convention: a-bj for absorbing
a = layer([[air, 10.0]])
b = layer([[air, 1.0], [water, 9.0]])
#b = layer([[water, 10.0]])
c = layer([[air, 10.0]])

geometry = stack([[a,1],[b,1.7],[c,1]][::-1])
geometry.paint().show()
waveguide = afmm(geometry, order=1, lambda0=1.0)
waveguide.compute()
waveguide.inputmode([0]*1+[1.0]+[0]*1)
#print sum(np.power(np.abs(waveguide.d[0]),2))


#print waveguide.d[0]
#for i in range(10):
#    print waveguide.Hy(5, i)
#print waveguide.u[2]
#print waveguide.d[2]
#print waveguide.W[2]
#print waveguide.Hy(0.5,10.99)
#print waveguide.T_dd[0]
#print 'W1=',waveguide.W[1]
#print 'W0=',waveguide.W[0]
print 'd[0]=', waveguide.d[0]
print 'd[1]=', waveguide.d[1]
print 'd[2]=', waveguide.d[2]
print 'u[0]=', waveguide.u[0]
print 'u[1]=', waveguide.u[1]
print 'u[2]=', waveguide.u[2]

##print 'T_dd[0]=', waveguide.T_dd[0]
wg = waveguide
#print 'X_1=', wg.X[1]
#print 'xu+x-1d=', np.add(np.dot(wg.X[1], wg.u[1]), np.dot(lg.inv(wg.X[1]), wg.d[1]))
print np.dot(wg.W[1], np.add(np.dot(wg.X[1], wg.u[1]), np.dot(lg.inv(wg.X[1]), wg.d[1])))
print np.dot(wg.W[2], np.add(wg.u[2], wg.d[2]))
print np.dot(wg.W[0], np.add(np.dot(wg.X[0], wg.u[0]), np.dot(lg.inv(wg.X[0]), wg.d[0])))
print np.dot(wg.W[1], np.add(wg.u[1], wg.d[1]))
print np.dot(wg.V[1], np.add(np.dot(wg.X[1], wg.u[1]), np.negative(np.dot(lg.inv(wg.X[1]), wg.d[1]))))
print np.dot(wg.V[2], np.add(wg.u[2], np.negative(wg.d[2])))

#print waveguide.R_ud[1]
#print waveguide.T_dd[1]
#print waveguide.modes(0)
