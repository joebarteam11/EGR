import sys,os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
import matplotlib 
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
path = os.getcwd()

#from lib_egr_260 import show_graphs
print('Current folder: ',path)
print(f"Running Matplotlib version: {matplotlib.__version__}")
Ncurves = 3
files=[
    '/1DNO_EGR_Aramco13.csv',
    '/1DCO2_Aramco13.csv',
    '/1DREAL_EGR_Aramco13.csv',
    '/1DNO_EGR_AP.csv',
    '/1DCO2_AP.csv',
    '/1DREAL_EGR_AP.csv',
    #'/1DREAL_EGR_EGR_0.1.csv',
    #'/1DCO2_EGR_0.1.csv',
    #'/CH4_15_256_9_AP.csv',
]
files=[path+
       '/results/'+
       file for file in files]

inputs=[pd.read_csv(file).round({'P':1,'EGR':1,'FB':1,'phi':2}) for file in files]


var_to_plot=['dF',
             'u',
             'T',
            ]
ylabels = ('Flame Thickness [um]',
            'SL0[m.s-1]',
            'T[K]',
          )
titles = ['Flame thickness',
          'Laminar flame speed',
          'Adiabatic flame temperature',
         ]

mech=[#'Polimi',
      #'Real EGR (dried) ARC',
      'No EGR',
      'CO2',
      'Real EGR',
      #'No EGR ARC',
      #'No EGR detailed',
      #'ARC(AP)',
      ]
symbols = ['o-','x-','s-','o--','x--','s--']
markers = ['o','x','s']

def update_prop(handle, orig):
    handle.update_from(orig)
    handle.set_marker("")
    handle.set_linestyle("-")
    handle.set_linewidth(15.0)

def update_prop2(handle, orig):
    handle.update_from(orig)
    handle.set_marker("")
    handle.set_linewidth(2.0)

plt.rcParams.update({'font.size': 15,
                        'lines.markersize': 10,
                        })

# for each mech (i.e. each input in inputs) and each variable in var to plot, plot one graph with phi as x axis and the variable as y axis
# each couple of EGR and FB (fuel blend) is a subplot
for i,var in enumerate(var_to_plot):
    fig,ax=plt.subplots(1,1,figsize=(10,10))
    labels=[]
    for j,input in enumerate(inputs):
        if(j==Ncurves):
            plt.gca().set_prop_cycle(None)
        data=input.pivot_table(index=['phi'],columns=['EGR','FB'],values=var)
        print(data)
        
        label = [r'$\ X_{EGR}^{fuel}$'+':'+str(val[0])+' '+'('+mech[j%3]+')'
                 #str(data.columns.names[1])+ 
                 #X,'H2' as index and 'fuel' as exponent
                 #+r'$\ X_{H2}^{fuel}$'+':'+str(val[1]) 
                 for val in data.columns.values]
        
        if(j<Ncurves):
            labels+=[label[i]for i in range(len(label))]

        if(var=='dF'):
            data=data*1e6
        d = data.plot(#x='phi',y=var,
                  ax=ax,style=symbols[j],
                  label=label,#color=colors[i*j]
                  linewidth=3.0,
                 )
        
    ax.legend(labels,loc='upper left',
                handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)},
                )
        
    plt.xlabel('Effective equivalence ratio')
    plt.ylabel(ylabels[i])
    plt.title(titles[i]
              +' - ('+r"$\bf{"+'T_{in}EGR:'+str(round(input['Tin'][0],1))+'K'+ "}$"+')'
              )

    #add a second legend to show each marker for each mech
    ax2 = ax.twinx()
    # for k,_ in enumerate(inputs):
    #     #try:
    #     if(i<Ncurves):    
    ax2.plot([],[],symbols[0],
            label='Detailed', 
            c='black',
            )
    ax2.plot([],[],symbols[3],
        label='ARC (AP)', 
        c='black',
        )
        #except:
        #    pass
    ax2.get_yaxis().set_visible(False)
    ax2.legend(loc='best',handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop2)},)
    ax.grid()
    
    #plt.legend(label,loc='upper left')
    # plt.rcParams.update({'font.size': 15,
    #                      'lines.markersize': 10,
    #                      })
    # if(var=='dF'):
    #     plt.rcParams.update({'font.size': 25,
    #                         'lines.markersize': 15,
    #                         })
    plt.savefig(path+'/img/'+
                'COMP_EGRtypes_GRI30_ARC'+
                '_'+var+'.png')

    
