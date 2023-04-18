#this script loads cantera freeflame xml files and plot them
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import sys

path = os.getcwd()
print('Current folder: ',path)
species = ['CH4','OH','H2','O2','CO','CO2','H2O']
xaxis = np.arange(len(species))
color = ['b','r','g','k','m','c','y']
files=[
      'egr0.0_phi1.2_T300.0_P5.0.xml',
      'egr0.1_phi1.2_T300.0_P5.0.xml',
      'egr0.3_phi1.2_T300.0_P5.0.xml',
      'egr0.5_phi1.2_T300.0_P5.0.xml',
     ]
filepath=[path+'/src/data/'+file for file in files]

#legend = ['EGR=0.0','EGR=0.1','EGR=0.3','EGR=0.5']
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

gas=ct.Solution('schemes/Aramco13.cti')
f=ct.FreeFlame(gas,width=0.03)
#legend = ['%CO2:{}%;{}:{}bar {}'.format(round(x[1]*100,0), f.columns.names[0] ,round(x[0]/100000,0),'Aramco1.3') for x in f.columns]

#create an offset function to slide each bar
step=1/(len(files)+1)
offset=-step*len(files)/2+step/2
legend=[]

phis = [float(file.split('_')[1].split('phi')[-1]) for file in files]
egr =[float(file.split('_')[0].split('egr')[-1]) for file in files]
pressures=[float(file.split('_')[3].split('P')[-1][:-4]) for file in files]
print(pressures)

for idx,file in enumerate(filepath):
    file = path+'/src/data/egr0.1_phi{}_T300.0_P1.0.xml'.format(phi)
    f.restore(file,loglevel=0)
    #index = [f.gas.species_index(specie) for specie in species]
    legend.append('%CO2:{}%;P:{}bar {}'.format(round(egr[idx]*100,0), round(pressures[idx],0),'Aramco1.3'))
    
    #offset=offset+step
plt.plot(X,f.outlet.phase.cp_mass,width=step,color=color[idx])
print(legend)
#extract a list of phi values in files names


title = 'Molar fractions (in outlet) of species at Phi={}'.format(phis[0])+' ('+r"$\bf{"+'T_{in}CO2:'+str(300)+'K'+ "}$"+')'
#flamme_thickness(f)
#X=f.grid
#Y=f.T
#function that compute the flamme thickness
#plt.plot(X,Y,'.-')

##file = path+ '/src/data/egr0.3_phi1.0_T300.0_P1.0.xml'
#f.restore(file,loglevel=0)
#index = [f.gas.species_index(specie) for specie in species]
#plt.bar(xaxis+0.2,f.X[index,-1],0.4)
#flamme_thickness(f)
#X=f.grid
#Y=f.T
#plt.plot(X,Y,'.-')

#plt.savefig(path+'/img/testxml.png')
plt.title(title)
plt.legend(legend)
plt.xticks(xaxis,species)
plt.ylabel('Molar fraction [-]')
plt.grid()
plt.savefig(path+'/img/barplot_phi={}_P={}.png'.format(phis[0],pressures[0]))