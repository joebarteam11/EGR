import sys,os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter,AutoMinorLocator
import matplotlib 
import numpy as np

sys.path.append(os.getcwd())
path = os.getcwd()
fs=20

print('Current folder: ',path)
print(f"Running Matplotlib version: {matplotlib.__version__}")

files=[
    #'plan_partiel_0D_dilutionKP_20230516-114346.csv',
    #'plan_partiel_0D_dilutionKP_BFER_20230516-132035.csv',
    #'plan_total_equilibrium_20230516-163815.csv',
    #'plan_total_equilibrium_20230516-164522.csv',
    'plan_total_equilibrium_KP_fco2.csv',

]
files=[path+'/results/'+file for file in files]

inputs=pd.concat([pd.read_csv(file).round({'P':-1,'EGR':2,'T':-1}) for file in files])

var_to_plot=[
            'O2',
            'CO',
            'CO2',
             #'CO',
             
            ]
ylabels = (  '$X_{O2}$',
             '$X_{CO}$',
             #'$X_{CO}$',
             '$X_{CO2}$',
          )
titles = ['Molar fraction',
          'Molar fraction',
          'Molar fraction',
         ]

mech=['Aramco1.3','BFER','GRIMech3.0']
colors=['b','r','g','k','m','y','c']
style=['o-','x--']
X_var='EGR'
mydata=inputs.pivot_table(index=X_var,columns=['P','T',],values=var_to_plot)
#print(mydata.columns.get_level_values('EGR').unique())

for idx,rate in enumerate(mydata.columns.get_level_values('T').unique()):
    #Figure and subplots
    fig,ax1=plt.subplots(1,1,figsize=(15,10))
    fig.subplots_adjust(left=0.3)
    axes=[ax1]
    for i in enumerate(var_to_plot[1:]):
        axes.append(ax1.twinx())

    #data and min max for scaling graphs
    tempdata=mydata.loc[:,mydata.columns.get_level_values('T')==rate]
    mins=[tempdata[[var]].min().min() for var in var_to_plot]
    maxs=[tempdata[[var]].max().max() for var in var_to_plot]

    #Plotting
    tempdataP=[]
    graphs=[]
    tempdataP+=[tempdata.loc[:,tempdata.columns.get_level_values('P')==P] for P in tempdata.columns.get_level_values('P').unique()]
    graphs+=[data.plot(y=var,ax=axes[i],style=style[j],color=colors[i],legend=False) for i,var in enumerate(var_to_plot) for j,data in enumerate(tempdataP)]
    

    #Y axis handling
    [ax.spines["left"].set_position(("axes", -0.17*(i+1)-(i-1)*0.00)) for i,ax in enumerate(axes[1:])]
    [ax.set_ylim(mins[i]-abs(mins[i]-maxs[i])*0.05,maxs[i]+abs(mins[i]-maxs[i])*0.05) for i,ax in enumerate(axes)]
    [ax.set_ylabel(ylabels[i],fontsize=fs) for i,ax in enumerate(axes)]
    [ax.yaxis.label.set_color(colors[j]) for j,ax in enumerate(axes)]

    tkw = dict(size=7, width=1.5)
    [ax.tick_params(axis='y', colors=colors[j], labelsize=fs-5, **tkw) for j,ax in enumerate(axes)]
    [ax.yaxis.set_minor_locator(AutoMinorLocator())for ax in axes]
    [ax.tick_params(which='minor', length=4) for ax in axes]
    [ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0),useMathText=True) for ax in axes]
    #[ax.yaxis.set_major_formatter(FormatStrFormatter('%.5e')) for ax in axes]
    
    ax1.tick_params(axis='x', labelsize=fs-5, **tkw)
    ax1.ticklabel_format(axis="x", style="plain")
    ax1.set_xlabel('$X_{CO2in}$',fontsize=fs)
    #ax1.set_xlim(mydata[[X_var]].min().min(),mydata[[X_var]].max().max())

    [ax.spines["left"].set_visible(True) for ax in axes[1:]]
    [ax.yaxis.set_label_position('left') for ax in axes[1:]]
    [ax.yaxis.set_ticks_position('left') for ax in axes[1:]]
    [ax.get_yaxis().get_offset_text().set_x(-0.17*(i+1)-(i-1)*0.00) for i,ax in enumerate(axes[1:])]

    
    #Legend handling
    labels=['P:{}bar'.format(round(x/100000,0)) for x in tempdata.columns.get_level_values('P').unique()]
    ax1.legend(labels,loc='center right',bbox_to_anchor=(1.07, 0.5),fontsize=fs-2)
    legs = ax1.get_legend()
    [leg.set_color('black') for leg in legs.legend_handles]

    ax1.grid()
    plt.rcParams.update({'font.size': fs-5})
    plt.title('Equilibrium computation - Scheme: '+mech[0]+' ('+r"$\bf{"+'T_{in}CO2+Air:'+str(1700)+'K'+ "}$"+')',fontsize=fs)
    #plt.show(block=True)
    plt.savefig('AIR_CO2_{}_KP_ARAMCO13.png'.format(1700),dpi=300)