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
Ncurves = 2
files=[

    '/1DCO2_AP.csv',
    '/1DFCO2_AP.csv',
    '/1DCO2_8AP.csv',
    '/1DFCO2_8AP.csv',

]
files=[path+
       '/results'+
       file for file in files]

inputs=[pd.read_csv(file).round({'P':1,'EGR':1,'FB':1,'phi':2}) for file in files]


var_to_plot=[#'dF',
             'u',
             #'T',
             #'CO2'
            ]
ylabels = (#'Flame Thickness [um]',
           #r"${"+'S_{L0}[m.s^{-1}]'+ "}$",
           r"${\frac{S_{L0_F}-S_{L0}}{S_{L0}}"+'[-]'+ "}$",
           # 'T[K]',
            #r"$\ X_{CO_2}$",
          )
titles = [#'Flame thickness',
          'Laminar flame speed relative error (CO2-FCO2)',
          #'Adiabatic flame temperature',
          #'CO2 content in flue gas',
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
symbols = ['o-','x--']
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

plt.rcParams.update({'font.size': 20,
                        'lines.markersize': 10,
                        'font.family': 'sans-serif',
                        'text.usetex': True,
                        })

# for each mech (i.e. each input in inputs) and each variable in var to plot, plot one graph with phi as x axis and the variable as y axis
# each couple of EGR and FB (fuel blend) is a subplot
for i,var in enumerate(var_to_plot):
    fig,ax=plt.subplots(1,1,figsize=(10,10))
    # plot curves representing the relative errors between two inputs

    data=[]
    for j,input in enumerate(inputs):
        # if(j==1):
        #     plt.gca().set_prop_cycle(None)
        #     print("restart color cycle")
        data.append(input.pivot_table(index=['phi'],columns=['EGR',],values=var))
    
    #create a new dataframe with the relative errors
    #relative_err = pd.DataFrame()
    # for j in range(len(data)-1):
        #add phi in the new dataframe
    for m in [0,2]:
        plt.gca().set_prop_cycle(None)
        relative_err = pd.DataFrame()
        relative_err['phi']=data[m].index.values
        relative_err['EGR']=[data[m].columns.values[0] for _ in range(len(data[m].index.values))]
        relative_err['error']=(abs(data[m+1]-data[m])/data[m]).values.flatten()

        relative_err=relative_err.pivot_table(index=['phi'],columns=['EGR'],values='error')
        #plot the relative errors
        relative_err.plot(ax=ax,style=symbols[m//2],linewidth=3.0,)

    # garbages
    labels=[]
    # for j,input in enumerate(inputs):
    #     if(j==Ncurves):
    #         plt.gca().set_prop_cycle(None)
    #         print("restart color cycle")
    #     data=input.pivot_table(index=['phi'],columns=['EGR','FB'],values=var)
    #     #data=data.loc[:,data.columns.get_level_values('P')<0.1]
    #     #data = data.replace({np.nan: np.inf})
    #     #data = data.dropna(subset=['EGR'])
    #     print(data)
    val = [0.2] #data.columns.values
    label = [
                r'$\ \%EGR_{mol}$'+':'+str(val[0]*100)+' ',
                #  +'('+
                #  mech[j%2]
                #  +')'
                #str(data.columns.names[1])+ 
                #X,'H2' as index and 'fuel' as exponent
                #+r'$\ X_{H2}^{fuel}$'+':'+str(val[1]) 
                #for val in [0.2] #data.columns.values
                ]
        
    #if(j<Ncurves):
    labels+=[label[i]for i in range(len(label))]

    #     if(var=='dF'):
    #         data=data*1e6
    #     d = data.plot(#y=var,
    #               ax=ax,style=symbols[j],
    #               label=label,#color=colors[i*j]
    #               linewidth=3.0,
    #               #use_index=True,
    #              )
    #     # plot data with line even if there is no data for a given phi
        
    ax.legend(labels,loc='upper left',
                handler_map={plt.Line2D:HandlerLine2D(update_func=update_prop)},
                )
        
    plt.xlabel(r'Effective equivalence ratio')
    plt.ylabel(ylabels[i])
    plt.title(titles[i]
              +' - ('+r"$\bf{"+'T_{in}EGR:'
              #+str(round(input['Tin'][0],1))
              +str(393)
              +'K'+ "}$"+')',
                fontsize=20
              )

    #add a second legend to show each marker for each mech
    ax2 = ax.twinx()
    #for k,_ in enumerate(inputs):
        #try:
    #if(i<Ncurves):    
    ax2.plot([],[],symbols[0],
            label='1bar', 
            c='black',
            )
    ax2.plot([],[],symbols[1],
        label='8bar', 
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
                'FCO2_relative_err'+
                '_'+var+'.png')

    
