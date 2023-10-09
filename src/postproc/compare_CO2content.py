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
    # '/1DNO_EGR_Aramco13.csv',
    # '/1DCO2_Aramco13.csv',
    # '/1DREAL_EGR_Aramco13.csv',
    #'/1DCO2_Poli1bar_M12.csv',
    #'/1DCO2_Poli5bar_M12.csv',
    #'/1DCO2_AP1bar_M12.csv',
    '/1DREAL_EGR_AP1bar0_M12.csv',
    '/1DREAL_EGR_AP1bar1_M12.csv',
    '/1DREAL_EGR_AP1bar2_M12.csv',
    #'/1DCO2_AP5bar_M12.csv',
    # '/1DCO2_AP.csv',
    # '/1DREAL_EGR_AP.csv',
    #'/1DNO_EGR_H2_AP.csv',
    #'/1DREAL_EGR_EGR_0.1.csv',
    #'/1DCO2_AP5bar_M12.csv',
    #'/CH4_15_256_9_AP.csv',
]
files=[path+
       '/results'+
       file for file in files]

inputs=[pd.read_csv(file).round({'P':1,'EGR':1,'FB':1,'phi':2}) for file in files]


var_to_plot=[#'dF',
             #'u',
             #'T',
             'CO2'
            ]
ylabels = (#'Flame Thickness [um]',
           # r"${"+'S_{L0}[m.s^{-1}]'+ "}$",
           # 'T[K]',
            r"$\ X_{CO_2}$",
          )
titles = [#'Flame thickness',
          #'Laminar flame speed',
          #'Adiabatic flame temperature',
          'CO2 content in flue gas',
         ]

mech=[#'Polimi',
      #'Real EGR (dried) ARC',
      #'No EGR',
      #'CO2',
      #'Real EGR',
      #'No EGR ARC',
      #'No EGR detailed',
      #'ARC',
      '1bar',
      '5bar',
      ]
symbols = ['o-','o-','o-','x--']
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
            print("restart color cycle")
        data=input.pivot_table(index=['phi'],columns=['EGR','FB'],values=var)
        #data=data.loc[:,data.columns.get_level_values('P')<0.1]
        #data = data.replace({np.nan: np.inf})
        #data = data.dropna(subset=['EGR'])
        print(data)
        
        label = [
                 r'$\ \%EGR_{mol}$'+':'+str(val[0]*100)+' '
                #  +'('+
                #  mech[j%2]
                #  +')'
                 #str(data.columns.names[1])+ 
                 #X,'H2' as index and 'fuel' as exponent
                 #+r'$\ X_{H2}^{fuel}$'+':'+str(val[1]) 
                 for val in data.columns.values]
        
        if(j<Ncurves):
            labels+=[label[i]for i in range(len(label))]

        if(var=='dF'):
            data=data*1e6
        d = data.plot(#y=var,
                  ax=ax,style=symbols[j],
                  label=label,#color=colors[i*j]
                  linewidth=3.0,
                  #use_index=True,
                 )
        # plot data with line even if there is no data for a given phi
        
    ax.legend(labels,loc='upper left',
                handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)},
                )
        
    plt.xlabel('Effective equivalence ratio')
    plt.ylabel(ylabels[i])
    plt.title(titles[i]
              +' - ('+r"$\bf{"+'T_{in}EGR:'
              #+str(round(input['Tin'][0],1))
              +str(573)
              +'K'+ "}$"+') - 1bar'
              )

    #add a second legend to show each marker for each mech
    ax2 = ax.twinx()
    # for k,_ in enumerate(inputs):
    #     #try:
    #     if(i<Ncurves):    
    # ax2.plot([],[],symbols[0],
    #         label='CO2', 
    #         c='black',
    #         )
    ax2.plot([],[],symbols[1],
        label='EGR', 
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
                'M12_REALEGR_1bar'+
                '_'+var+'.png')

    
