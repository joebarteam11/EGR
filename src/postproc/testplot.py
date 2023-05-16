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
    'plan_total_equilibrium_20230516-164522.csv',

]
files=[path+'/results/'+file for file in files]

inputs=pd.concat([pd.read_csv(file).round({'P':-1,'EGR':2,'T':-1}) for file in files])
#papier=pd.read_csv(path+'/results/'+'data-handmade.csv',delimiter=';').round({'EGR':1,'phi':2})
#papier=pd.read_csv(path+'/plan_total_dilution_BFERUNITY_20230412-170317.csv').round({'P':1,'EGR':1,'phi':2})
#print(input)

var_to_plot=['CO',
             'CO2',
             'O2',
            ]
ylabels = (  '$X_{CO}$',
             '$X_{CO2}$',
             '$X_{O2}$',
          )
titles = ['Molar fraction',
          'Molar fraction',
          'Molar fraction',
         ]

mech=['Aramco1.3','BFER']
colors=['b','r','g','k','m','y','c']
style=['o-','x--']

mydata=inputs.pivot_table(index='T',columns=['P','EGR',],values=var_to_plot)
print(mydata.columns.get_level_values('EGR').unique())

for idx,rate in enumerate(mydata.columns.get_level_values('EGR').unique()):
    fig,ax1=plt.subplots(1,1,figsize=(15,10))
    fig.subplots_adjust(left=0.3)

    # if(idx<2):
    #     idx1=0
    #     idx2=idx
    # else:
    #     idx1=1
    #     idx2=idx-2

    tempdata=mydata.loc[:,mydata.columns.get_level_values('EGR')==rate]
    print(tempdata)
    
    minCO=tempdata[[var_to_plot[0]]].min().min()
    minCO2=tempdata[[var_to_plot[1]]].min().min()
    minO2=tempdata[[var_to_plot[2]]].min().min()

    maxCO=tempdata[[var_to_plot[0]]].max().max()
    maxCO2=tempdata[[var_to_plot[1]]].max().max()
    maxO2=tempdata[[var_to_plot[2]]].max().max()
    # print(maxCO)
    # print(maxCO2)
    # print(maxO2)

    twin1 = ax1.twinx()
    twin2 = ax1.twinx()
    twin1.spines["left"].set_position(("axes", -0.15))
    twin2.spines["left"].set_position(("axes", -0.30))

    #if(maxCO>0):
    offset = (maxCO - minCO) * 0.05
    ax1.set_ylim(minCO-offset,maxCO+offset)
    #if(maxCO2>0):
    offset = (maxCO2 - minCO2) * 0.05
    twin1.set_ylim(minCO2-offset,maxCO2+offset)
    #if(maxO2>0):
    offset = (maxO2 - minO2) * 0.05
    twin2.set_ylim(minO2-offset,maxO2+offset)

    
    #twin1.set_ylabel(ylabels[1])
    #twin2.set_ylabel(ylabels[2])

    axes=[ax1,
          twin1,twin2,
          ]
    coloridx=len(var_to_plot)

    tempdataP=[]
    graphs=[]
    tempdataP+=[tempdata.loc[:,tempdata.columns.get_level_values('P')==P] for P in tempdata.columns.get_level_values('P').unique()]
    graphs+=[data.plot(y=var,ax=axes[i],style=style[j],color=colors[i],legend=False) for i,var in enumerate(var_to_plot) for j,data in enumerate(tempdataP)]
    [ax.set_ylabel(ylabels[i]) for i,ax in enumerate(axes)]
    [ax.yaxis.label.set_color(colors[j]) for j,ax in enumerate(axes)]


    tkw = dict(size=4, width=1.5)
    [ax.tick_params(axis='y', colors=colors[j], **tkw) for j,ax in enumerate(axes)]
    [ax.yaxis.set_minor_locator(AutoMinorLocator())for ax in axes]
    [ax.tick_params(which='minor', length=4) for ax in axes]
    #[ax.set_yticks(ticks) for ticks in yticks for ax in axes]
    [ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0),useMathText=True) for ax in axes]
    #[ax.yaxis.set_major_formatter(FormatStrFormatter('%.5e')) for ax in axes]
    ax1.tick_params(axis='x', **tkw)
    ax1.set_xlabel('Temperature [K]')

    twin1.spines["left"].set_visible(True)
    twin1.yaxis.set_label_position('left')
    twin1.yaxis.set_ticks_position('left')  
    twin1.get_yaxis().get_offset_text().set_x(-0.2) 
    #twin1.get_yaxis().get_offset_text().set_y(1.2) 
    twin2.spines["left"].set_visible(True)
    twin2.yaxis.set_label_position('left')
    twin2.yaxis.set_ticks_position('left')
    twin2.get_yaxis().get_offset_text().set_x(-0.4)
    #twin2.get_yaxis().get_offset_text().set_y(2) 
    plt.grid()

    labels=['P:{}bar'.format(round(x/100000,0)) for x in tempdata.columns.get_level_values('P').unique()]
    ax1.legend(labels,loc='best',bbox_to_anchor=(1.1, 1.1))
    legs = ax1.get_legend()
    [leg.set_color('black') for leg in legs.legendHandles]

    plt.title('$X_{CO2in}$'+': {}% - '.format(rate*100)+ mech[1]+' scheme')
    #plt.show(block=True)
    plt.savefig('AIR_CO2_{}_KP_Aramcoeq.png'.format(rate),dpi=300)