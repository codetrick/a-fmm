#! /usr/bin/env python2
# Aperiodic Fourier Modal Method written by Fan Yang

import numpy as np
import scipy as sp
import scipy.linalg 
import numpy.linalg as lg
import matplotlib.pyplot as plt
import Image
import ImageColor

class material:
    def __init__(self, eps_xx, eps_zz=None, mu=1.0, color='white', name=None):
        if eps_zz==None:
            eps_zz = eps_xx
        if callable(eps_xx):
            self.eps_xx = eps_xx
        else:
            self.eps_xx = lambda x: eps_xx

        if callable(eps_zz):
            self.eps_zz = eps_zz
        else:
            self.eps_zz = lambda x: eps_zz

        if callable(mu):
            self.mu = mu
        else:
            self.mu = lambda x: mu

        if callable(color):
            self.color = color
        else:
            self.color = lambda x: ImageColor.getrgb(color)

class layer:
    def __init__(self, pieces): # pieces = [ [material, width], ... ]
        self.pieces = pieces

        self.xmesh = []
        wc = 0.0
        self.xmesh.append(wc)
        for m in pieces:
            wc += m[1]
            self.xmesh.append(wc)
        self.w = self.xmesh[-1]
        self.xmesh = np.array(self.xmesh)

    def eps_xx(self, x):
        n = self.xmesh.searchsorted(x,'right')-1
        return self.pieces[n][0].eps_xx(x-self.xmesh[n])

    def eps_zz(self, x):
        n = self.xmesh.searchsorted(x,'right')-1
        return self.pieces[n][0].eps_zz(x-self.xmesh[n])
    
    def mu(self, x):
        n = self.xmesh.searchsorted(x,'right')-1
        return self.pieces[n][0].mu(x-self.xmesh[n])

    def color(self, x):
        n = self.xmesh.searchsorted(x,'right')-1
        return self.pieces[n][0].color(x-self.xmesh[n])

    def __add__(a, b):
        return layer(a.pieces+b.pieces)

class stack:
    def __init__(self, layers): # layers = [ [layer, thick], ... ]
        self.layers = map(lambda x:x[0], layers)
        self.thick = map(lambda x:x[1], layers)
        self.w = layers[0][0].w
        self.h = sum(self.thick)

    def f_mu(self):
        return [layer.mu for layer in self.layers]

    def f_eps_xx(self):
        return [layer.eps_xx for layer in self.layers]

    def f_eps_zz(self):
        return [layer.eps_zz for layer in self.layers]

    def getcolor(self, n_layer, x):
        return self.layers[n_layer].color(x)

    def paint(self,d=None):
        T = self.w
        H = self.h
        # default resolution
        if d == None:
            d = self.w/100.0
        slices = self.thick
        ymesh = []
        y = 0
        for h in slices:
            y += h
            ymesh += [y]
        ymesh = np.array(ymesh)
        ny = int(ymesh[-1] / float(d)) - 1
        nx = int(T/float(d)) -1 
        data = np.zeros([nx,ny,3], dtype = np.uint8)
        y = d
        for i in range(ny):
            n_layer = ymesh.searchsorted(y)
            x = d
            for k in range(nx):
                data[k,i,:] = self.getcolor(n_layer, x)
                x=x+d
            y=y+d

        return Image.frombuffer('RGB',(ny,nx),data,'raw','RGB',0,1).transpose(2)

    def __add__(a, b):
        return stack([[a.layers[n], a.thick[n]] for n in range(len(a.layers))] + [[b.layers[n], b.thick[n]] for n in range(len(b.layers))])

class afmm:
    def __init__(self, stack, order=12, lambda0=600.0, eps0=1.0, mu0=1.0):
        self.stack = stack
# width of slices
        self.w = stack.w
# electric & magnetic susceptibility
        self.f_eps_xx = stack.f_eps_xx()
        self.f_eps_zz = stack.f_eps_zz()
        self.f_mu = stack.f_mu()
        self.thick = stack.thick
        self.lambda0 = lambda0
        self.order = order
        self.K = 2*np.pi/self.w
        self.eps0 = eps0
        self.mu0 = mu0
        self.k0 = 2*np.pi/self.lambda0
        zmesh = [0.0]
        z = 0
        for h in self.thick:
            z += h
            zmesh += [z]
        self.zmesh = np.array(zmesh)
        
    def toeplitz(self, func, width, order):
        n = order*2+1
# FFT convention, see numpy documentation
        a = np.fft.fft([func(m*self.w/n) for m in range(n)])/n
        a = np.concatenate((a, [a[0]]))
        return sp.linalg.toeplitz(a[0:order+1], a[n+1:order:-1])

    def compute(self):
        k0 = self.k0
        w = self.w
        order = self.order
        E = []
        A = []
        B = []
        Kx = []
        S = []
        W = []
        V = []
        Lamb = []
        R_ud = []
        T_dd = []
        self.X = []

        for i in range(len(self.thick)):
            E.append(self.toeplitz(self.f_eps_zz[i], w, order*2))
            A.append(self.toeplitz(lambda x: 1.0/self.f_eps_xx[i](x), w, order*2))
            B.append(self.toeplitz(self.f_mu[i], w, order*2))
            Kx.append(np.diag(range(-order, order+1))*2.0*np.math.pi/w/k0)
            S.append(np.dot(lg.inv(A[i]), np.subtract(np.dot(np.dot(Kx[i],lg.inv(E[i])),Kx[i]),B[i])))
            lam, ww = lg.eig(S[i])
            Lamb.append(np.sqrt(lam))
            for j in range(order*2+1):
                if Lamb[i][j].imag>0 and abs(Lamb[i][j].real)<abs(Lamb[i][j].imag):
                    #print 'Lamb neg ', i, ' is ', Lamb[i][j]
                    Lamb[i][j] = -Lamb[i][j]
            W.append(ww)
            V.append(np.dot(np.dot(A[i],W[i]),np.diag(Lamb[i])))
        #    print 'Lamb=',Lamb[i]
           # print np.divide(np.dot(S[i], ww), np.dot(ww, np.diag(lam)))
           # print 'S_',i, '=', S[i]
           # print 'ok'

        #print 'E=',E
        #print 'A=',A
        #print 'B=',B
        #print 'Kx',Kx

        for i in range(len(self.thick)-1):
            X = np.diag(np.exp(np.negative(Lamb[i]*self.thick[i])*k0))
            self.X.append(X)
        #    print 'X_',i,'=', X
            #print 'Lambda_',i,'=', Lamb[i]
            #print 'W_',i,'=', W[i]
            wwvv_inv = lg.inv(np.add(np.dot(lg.inv(W[i]),W[i+1]),np.dot(lg.inv(V[i]),V[i+1])))
            t_uu = 2*np.dot(wwvv_inv, X)
            r_ud = np.dot( wwvv_inv, \
                           np.subtract(np.dot(lg.inv(V[i]),V[i+1]), np.dot(lg.inv(W[i]),W[i+1])))
            vvww_inv = lg.inv(np.add(np.dot(lg.inv(W[i+1]),W[i]),np.dot(lg.inv(V[i+1]),V[i])) )
        #    print 'vvww_inv=', vvww_inv
            r_du = np.dot(np.dot(X, np.dot( vvww_inv, \
                           np.subtract(np.dot(lg.inv(V[i+1]),V[i]),np.dot(lg.inv(W[i+1]),W[i])) )), X)
            t_dd = 2*np.dot(X, vvww_inv)
        #    print 't_dd',i,'=',t_dd
        #    print 'r_du',i,'=',r_du
        #    print 'r_ud',i,'=',r_ud
        #    print 't_uu',i,'=',t_uu
            if i==0:
                R_ud.append(r_ud)
                T_dd.append(t_dd)
                continue
            #RT_inv = lg.inv(np.add(np.identity(2*order+1),np.negative(np.dot(r_du,R_ud[i-1]))))
            RT_inv = np.subtract(np.identity(2*order+1),np.dot(r_du,R_ud[i-1]))
        #    print 'RT=', RT_inv
            RT_inv = lg.inv(RT_inv)
        #    print 'RT_inv=', RT_inv
            R_ud.append(np.add(r_ud ,np.dot(np.dot(np.dot(t_uu, R_ud[i-1]), RT_inv), t_dd)))
            T_dd.append(np.dot(np.dot(T_dd[i-1],RT_inv), t_dd))

        self.W = W
        self.V = V
        self.Lamb = Lamb
        self.R_ud = R_ud
        self.T_dd = T_dd

    def inputmode(self, d_in):
        n = len(self.thick)
        d = [None] * n
        u = [None] * n
        d[n-1] = d_in
        u[0] = np.array([0.0]*(self.order*2+1))
        d[0] = np.dot(self.T_dd[n-2], d[n-1])
        u[n-1] = np.dot(self.R_ud[n-2], d[n-1])
        for i in range(n-2):
            d[i+1] = np.dot(lg.inv(self.T_dd[i]), d[0])
            u[i+1] = np.dot(self.R_ud[i], d[i+1])
        self.d = d
        self.u = u

    def input(self, Hy):
        pass

    def modes(self, n):
        return self.Lamb[n]

    def Hy(self, x, z):
        n = self.zmesh.searchsorted(z,'right')-1
        zz = self.zmesh[n]
        #print 'n=',n
        #print 'zz=', zz
        #print 'z=', z
        #print 'd=',self.d[n]
        exp_u = np.multiply(self.u[n], np.exp(self.Lamb[n]*(zz-z)*self.k0))
        exp_d = np.multiply(self.d[n], np.exp(self.Lamb[n]*(z-zz)*self.k0))
        #print 'exp_d=', exp_d
        U = np.dot(self.W[n], np.add(exp_u, exp_d))
        #print 'U=', U
        hy = np.dot(U, np.exp(np.arange(-self.order, self.order+1)*x*1J*self.K))
        return hy

    def plotHy(self, z, frac=200):
        x = [i*1.0*self.w/frac for i in range(frac)]
        y = [None]*frac
        for i in range(frac):
            y[i] = self.Hy(x[i], z)
        plt.plot(x, y)
        plt.show()

    def plotinput(self, frac=200):
        x = [i*1.0*self.w/frac for i in range(frac)]
        y = [None]*frac
        n = len(self.thick)-1
        for i in range(frac):
            y[i] = np.dot(np.dot(self.W[n], self.d[n]), np.exp(np.arange(-self.order, self.order+1)*x[i]*1J*self.K))
        plt.plot(x, y)
        plt.show()

    def plotreflect(self, frac=200):
        x = [i*1.0*self.w/frac for i in range(frac)]
        y = [None]*frac
        n = len(self.thick)-1
        for i in range(frac):
            y[i] = np.dot(np.dot(self.W[n], self.u[n]), np.exp(np.arange(-self.order, self.order+1)*x[i]*1J*self.K))
        plt.plot(x, y)
        plt.show()

    def plotW(self, order, frac=200):
        x = [i*1.0*self.w/frac for i in range(frac)]
        y = [None]*frac
        n = len(self.thick)-1
        for i in range(frac):
            y[i] = np.dot(np.array(self.W[n])[:,order], np.exp(np.arange(-self.order, self.order+1)*x[i]*1J*self.K))
        plt.plot(x, y)
        plt.show()

    def Ex(self, x, z):
        n = self.zmesh.searchsorted(z, 'right')-1
        zz = self.zmesh[n]
        exp_u = np.multiply(self.u[n], np.exp(self.Lamb[n]*(zz-z)*self.k0))
        exp_d = np.multiply(self.d[n], np.exp(self.Lamb[n]*(z-zz)*self.k0))
        S = np.dot(self.V[n], np.subtract(exp_u, exp_d))
        ex = np.dot(S, np.exp(np.arange(-self.order, self.order+1)*x*1J*self.K))
        ex = ex * 1.0j * np.sqrt(self.mu0/self.eps0)
        return ex

