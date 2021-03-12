#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 16:25:46 2021

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

s =int(sys.argv[1])     #VAR file to start from
st =int(sys.argv[2])    #VAR file to stop at
n  =int(sys.argv[3])    #number of particles to pick out

var = pc.read_var(varfile='VAR1' )  #read in fisrt var 
dt = var.t              # gets time step to use latter

def patry(v=-1):        # this function stores all info of v-th pvar on array p
    if v == -1:
        pvar = pc.read_pvar()
    else:
        pvar = pc.read_pvar(varfile='PVAR'+str(v))
    i,x,y,z,vx,vy,vz =pvar.ipars,pvar.xp,pvar.yp,pvar.zp,pvar.vpx,pvar.vpy,pvar.vpz
    p=np.array([i,x,y,z,vx,vy,vz])
    return p

def get_party(s,st):   # sotres partilce data for selcted time slices          
    p = patry()
    P = np.zeros((st-s, p.shape[0],p.shape[1]))
    for i in range(s,st):
        P[i-s] = patry(i)
        print(i)
    return P

def deconv(x,lx ):      # for given array of particle positions x and domain length lx this function take cares of periodic bcs
    dx=x-np.roll(x,1)
    dx[0]=0
    for i in range(len(x)):
        if dx[i] > lx/2:
            x[i:] -= lx
        if dx[i] < -lx/2:
            x[i:] += lx
    return x

'''
this function picks out n particles at random. And uses there indexes to track them 
from pvar to pvar. It returns an array where the i-th coloum is the displacment of 
the i-th particle over time. from this the variance of the dicplacemnt(over the particles) 
can be found at any given time. also the mean drift can be found.
'''
def X_(s,st,n,p,lx):    
    inds=np.random.randint(0,1000000,n)
    I,X = p[s:st,0,:],p[s:st,1,:]
    X_ = np.zeros((st-s,len(inds)))
    for i in range(len(inds)):
        x = deconv(X[I==inds[i]],lx)
        X_[:,i] = x-x[0]
        print(i)
    return X_

param=pc.read_param()   
lx = param.lxyz[0]      # here we get the domain size

p = get_party(s,st )  # get particle info for desired time steps

xx = X_(0,st-s,n,p,lx )  # xx is the array of displacments

sig2a = np.var(xx,1 )  # this is the variacne of the displacment over the particles.
                        # by fitting this we get the diffusivity
                        
meanxt = np.mean(xx,1 ) # by fitting this we get the mean drift


'''
this creates directory diffusio_measerments and saves both time series from
 above for latter use or refitting
'''
if not os.path.isdir('./diffusion_measurement'):
    os.mkdir('./diffusion_measurement')
np.save('./diffusion_measurement/variance_of_displacment_timeseries',sig2a) 
np.save('./diffusion_measurement/mean_displacment',meanxt)
np.save('./diffusion_measurement/ts',np.arange(len(sig2a))*dt)


''' here the fit is done to get the diffusivity from sig2a. if sig2a is not 
changing linearly this values is useless. In which case the saved data can be 
used to refit to a better time span or what ever. 
'''
from sklearn.linear_model import LinearRegression
    #sig2a = np.var(e,1)
reg = LinearRegression(fit_intercept=False).fit(dt*np.arange(len(sig2a)).reshape((-1,1)),sig2a)


''' these two plots can be used to check if the fitt makes sense. A
 histogram of the displacments is shown for several time steps. these
 should be gaussian for fickian diffusion
 
'''
fig ,ax = plt.subplots(1,2)
ax[0].plot(dt*np.arange(len(sig2a)),sig2a ,'.',label = r'$VAR_p(x(t)-x(0))[H^{2}]$' )
ax[0].plot(dt*np.arange(len(sig2a)) ,reg.intercept_+reg.coef_*dt*np.arange(len(sig2a)),label='lin.fit' )

ax[0].legend( title = r'$\delta_{x}=$'+str(0.5*reg.coef_[0])+r'$H^{-}c_{s}^{-}$'   )


i = st -s -1

while i > 0:
        ax[1].hist(xx[i,:],label= str(int(i*dt))+r'$\Omega^{-}$' )
        i = i - int((st-s)/4 )

ax[0].set_xlabel(r'$t[\Omega^{-}]$')
ax[0].set_ylabel(r'$VAR_p(x(t)-x(0))[H^{2}]$')
ax[0].set_title('radial displacement variance in time' )

ax[1].set_title('radial displacement distribution')
ax[1].set_xlabel('(x(t)-x(0))[H] ')
ax[1].set_ylabel(r'$N_{p}$' ) 
ax[1].legend()
print('diffusivity:'+'c_s*H*delta_x =',reg.coef_[0]*0.5)
print('This value is valid only if change of displacement variance is linear. check plots to make sure that is the case. Otherwise refit using data saved in directory "Diffusion measurement", where the data for mean partilce drift has also been saved ')

plt.show()
