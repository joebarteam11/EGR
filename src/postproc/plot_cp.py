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

gas=ct.Solution('schemes/Aramco13.cti')
f=ct.FreeFlame(gas,width=0.03)

phis=[round(i,2) for i in np.arange(0.8,1.21,0.05)]
print(phis)
pressures=[1.0,5.0]
egr =[0.0,0.1,0.3,0.5]


for p in pressures :
    

    _,ax = plt.subplots(1,1,figsize=(10,10))
    #legend = ['%CO2:{}%;{}:{}bar {}'.format(round(x[1]*100,0), f.columns.names[0] ,round(x[0]/100000,0),'Aramco1.3') for x in f.columns]

    #create an offset function to slide each bar
    #step=1/(len(files)+1)
    #offset=-step*len(files)/2+step/2
    legend=[]
    for idx,rate in enumerate(egr):
        X=phis
        Y=[]
        files=[
            'egr{}_phi{}_T300.0_P{}.xml'.format(rate,phis[i],p) for i in range(len(phis))
            ]

        filepath=[path+'/src/data/'+file for file in files]
        #X.append(rate)
        legend.append('%CO2:{}%;P:{}bar {}'.format(round(rate*100,0), round(p,0),'Aramco1.3'))
        print(legend)
        for file in filepath:
            f=ct.FreeFlame(gas,width=0.03)
            #file = path+'/src/data/egr0.1_phi1.0_T300.0_P1.0.xml'
            f.restore(file,loglevel=0)
            #index = [f.gas.species_index(specie) for specie in species]
            print('Cp',f.cp_mass[-1]*f.density_mass[-1]/1000)
            Y.append(f.cp_mass[-1]*f.density_mass[-1]/1000)
            #offset=offset+step
        ax.plot(X,Y,'.-',color=color[idx])
        
         #extract a list of phi values in files names

    
    title = 'Mean volumetric heat capacity (outlet) '+' ('+r"$\bf{"+'T_{in}CO2:'+str(300)+'K'+ "}$"+')'
    fs=20

    ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.rcParams.update({'font.size': fs})
    plt.title(title,fontsize=fs)
    plt.legend(legend,fontsize=fs)
    #plt.xticks(xaxis,species)
    plt.xlabel('Phi',fontsize=fs)
    plt.ylabel('rho*Cp [kJ/(m3.K)]',fontsize=fs)
    plt.grid()
    plt.tight_layout()
    plt.savefig(path+'/img/plot_cp_P={}.png'.format(p),dpi=300,)