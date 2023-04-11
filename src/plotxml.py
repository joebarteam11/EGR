#this script loads cantera freeflame xml files and plot them
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import sys
def flamme_thickness(f):
    #get the temperature profile
    T=f.T
    #get the grid
    X=f.grid
    #get the temperature at the flame base
    T0=T[0]
    #get the temperature at the flame tip
    T1=T[-1]
    #compute the maximal spatial temperature gradient along the flame
    dTdx=np.max(np.gradient(T,X))
    #compute the thickness
    thickness=(T1-T0)/dTdx
    print('Thicknes (m)',thickness)
    return thickness

file = 'egr0.1_phi1.2_T300.0_P1.0.xml'
path = os.getcwd()
print('Current folder: ',path)
gas=ct.Solution('Aramco13.cti')
f=ct.FreeFlame(gas,width=0.03)
f.restore(file,loglevel=0)
flamme_thickness(f)
X=f.grid
Y=f.T
#function that compute the flamme thickness


plt.plot(X,Y,'.-')

file = 'egr0.3_phi1.2_T300.0_P1.0.xml'
f.restore(file,loglevel=0)
flamme_thickness(f)
X=f.grid
Y=f.T
plt.plot(X,Y,'.-')

plt.savefig(path+'/img/testxml.png')