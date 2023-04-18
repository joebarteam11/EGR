import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import sys

path = os.getcwd()
print('Current folder: ',path)
species1 = ['CO2','H2O','O2']##
species2 = ['CH4','OH','H2','CO']
allspecies = [species2,species1]

fig,ax = plt.subplots(1,2, figsize=(10,10),gridspec_kw={'width_ratios': [1.2, 1]})
#ax1  = fig.add_subplot(121)
#ax2  = fig.add_subplot(122)



color = ['b','r','g','k','m','c','y']

gas=ct.Solution('schemes/Aramco13.cti')
f=ct.FreeFlame(gas,width=0.03)

phis=[i for i in np.arange(0.8,1.21,0.1)]
pressures=[1.0,5.0]

for i,phi in enumerate(phis) :
    for j,p in enumerate(pressures):
        #_,ax = plt.subplots(2,1,figsize=(10,10))
        for n,species in enumerate(allspecies):
            xaxis = np.arange(len(species))

            files=[
                'egr0.0_phi{}_T300.0_P{}.xml'.format(phi,p),
                'egr0.1_phi{}_T300.0_P{}.xml'.format(phi,p),
                'egr0.3_phi{}_T300.0_P{}.xml'.format(phi,p),
                'egr0.5_phi{}_T300.0_P{}.xml'.format(phi,p),
                ]

            filepath=[path+'/src/data/'+file for file in files]
            
            #legend = ['%CO2:{}%;{}:{}bar {}'.format(round(x[1]*100,0), f.columns.names[0] ,round(x[0]/100000,0),'Aramco1.3') for x in f.columns]

            #create an offset function to slide each bar
            step=1/(len(files)+1)
            offset=-step*len(files)/2+step/2
            legend=[]

            egr =[float(file.split('_')[0].split('egr')[-1]) for file in files]


            for idx,file in enumerate(filepath):
                #file = path+'/src/data/egr0.1_phi1.0_T300.0_P1.0.xml'
                f.restore(file,loglevel=0)
                index = [f.gas.species_index(specie) for specie in species]
                #legend.append('%CO2:{}%;P:{}bar {}'.format(round(egr[idx]*100,0), round(p,0),'Aramco1.3'))
                legend.append('%CO2:{}%;{}'.format(round(egr[idx]*100,0),' Aramco1.3'))

                if(n==0):
                    ax[n].bar(xaxis+offset,f.X[index,-1],width=step,color=color[idx])
                elif(n==1):
                    ax[n].bar(xaxis+offset,f.X[index,-1],width=step,color=color[idx])
                #plt.bar(xaxis+offset,f.Y[index,-1],width=step,color=color[idx])
                offset=offset+step

            print(legend)
            #extract a list of phi values in files names


            title = 'Mole fractions (outlet) of species at Phi={} '.format(phi)+' ('+r"$\bf{"+'T_{in}CO2:'+str(300)+'K'+ "}$"+')'+" - P:{}bar".format(p)
            fs=20

            if(n==0):
                #ax1.s
                ax[n].set_xticks(xaxis,species,fontsize=fs)
            elif(n==1):
                ax[n].set_xticks(xaxis,species,fontsize=fs)

        plt.rcParams.update({'font.size': fs})
        [a.tick_params(axis='both', which='major', labelsize=fs) for a in ax]
        #ax2.tick_params(axis='both', which='major', labelsize=fs)
        fig.suptitle(title,fontsize=fs)
        ax[0].legend(legend,fontsize=15,loc='upper right')
        [a.set_ylabel('Mole fraction [-]',fontsize=fs) for a in ax]
        #ax2.set_ylabel('',fontsize=fs)
        #ax1.set_yticks()
        #plt.ylim(0.001,0.25)
        #plt.yscale('log')
        #ax.set_yticklabels(['0.001','0.01','0.1'])
        #ax.set_y
        # = ['0.005','0.010','0.05','0.1','0.15','0.2','0.25']
        fig.tight_layout()
        [a.grid() for a in ax]
        #ax2.grid()
        fig.savefig(path+'/img/mole_barplot_CO2H2O_phi={}_P={}.png'.format(phis[i],pressures[j]), bbox_inches='tight')