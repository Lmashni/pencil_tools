#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 17:15:52 2021

@author: lyth

This routine should print out the radial difusivity. For now it assumes a 3d simu
-lation with periodic boundries. To be called in run directory after the simulation
has ran e.g:

python rad_diff.py s st n

parameters:
    s:     int pvar file to start from
    st:    int pvar file to stop at 
    n:     int number of particles to track

 
the routine tracks the variance of the particle displacment over time and 
fits it linearly assuming fickian diffusion. the diffusivity is then given by:
    
    H * c_s * delta_x = 0.5 * d(sig^2)/dt

where sig^2 is the variance of the particle displacment at a given time and 
d(sig^2)/dt is for fickian diffusion the coefficient of the linear fit. A plot 
is preduced of the fit and the particle distrebution as it spreads out with time.
The routine needs suffciently many PVARs.

Nutzung auf eigene Gefahr, eltern haften für ihre Kinder.
"""

import numpy as np
import matplotlib.pyplot as plt
import pencil as pc
import os 
import sys



s =int(sys.argv[1])
st =int(sys.argv[2])
n  =int(sys.argv[3])

var = pc.read_var(varfile='VAR1' )
dt = var.t

def patry(v=-1):
    if v == -1:
        pvar = pc.read_pvar()
    else:
        pvar = pc.read_pvar(varfile='PVAR'+str(v))
    i,x,y,z,vx,vy,vz =pvar.ipars,pvar.xp,pvar.yp,pvar.zp,pvar.vpx,pvar.vpy,pvar.vpz
    p=np.array([i,x,y,z,vx,vy,vz])
    return p

def get_party(s,n):
    #P = np.zeros((n,patry().shape))
    p = patry()
    P = np.zeros((n-s, p.shape[0],p.shape[1]))
    for i in range(s,n):
        P[i-s] = patry(i)
        print(i)
    #np.save('./P',P)
    return P

def deconv(x,lx ):
    dx=x-np.roll(x,1)
    dx[0]=0
    for i in range(len(x)):
        if dx[i] > lx/2:
            x[i:] -= lx
        if dx[i] < -lx/2:
            x[i:] += lx
    return x

def X_(s,st,n,p,lx):
    inds=np.random.randint(0,1000000,n)
#    p = np.load('./diag/party.npy')
    I,X = p[s:st,0,:],p[s:st,1,:]
    X_ = np.zeros((st-s,len(inds)))
    for i in range(len(inds)):
        x = deconv(X[I==inds[i]],lx)
        #print(len(X[I==inds[i]]))
        #print(len(x ))  
        X_[:,i] = x-x[0]
        print(i)
    return X_

param=pc.read_param()
lx = param.lxyz[0]

p = get_party(s,st )

xx = X_(0,st-s,n,p,lx )

sig2a = np.var(xx,1 )

from sklearn.linear_model import LinearRegression
    #sig2a = np.var(e,1)
reg = LinearRegression(fit_intercept=False).fit(dt*np.arange(len(sig2a)).reshape((-1,1)),sig2a)

fig ,ax = plt.subplots(1,2)
ax[0].plot(dt*np.arange(len(sig2a)),sig2a ,'.',label= r'$\delta_{x}=$'+str(0.5*reg.coef_[0])+r'$H^{-}c_{s}^{-}$')
ax[0].plot(dt*np.arange(len(sig2a)) ,reg.intercept_+reg.coef_*dt*np.arange(len(sig2a)))

i = st -s -1

while i > 0:
        ax[1].hist(xx[i,:],label= str(int(i*dt))+r'$\Omega^{-}$' )
        i = i - int((st-s)/4 )

ax[0].set_xlabel(r'$t[\Omega^{-}]$')
ax[0].set_ylabel(r'$\sigma^{2}[H^{2}]$')
ax[0].set_title('radial displacement variance in time' )

ax[1].set_title('radial displacement distribution')
ax[1].set_xlabel('x(t)-x(0)[H] ')
ax[1].set_ylabel(r'$N_{p}$' ) 
plt.legend()
print('diffusivity:'+'c_s*H*delta_x =',reg.coef_[0]*0.5)
#plt.plot(dt*np.arange(len(sig2t)),sig2t ) 
plt.show()