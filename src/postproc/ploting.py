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

files=[
    '/CH4_16_250_10_QC.csv',
    '/CH4_18_444_12_DL.csv',
    #'/CH4_15_256_9_AP.csv',
]
files=[path+
       #'/EGR0D'+
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
      'ARC(QC16)',
      'ARC(DL)',
      #'ARC(AP)',
      ]
symbols = ['o-','x-','s-',]



# for each mech (i.e. each input in inputs) and each variable in var to plot, plot one graph with phi as x axis and the variable as y axis
# each couple of EGR and FB (fuel blend) is a subplot
for i,var in enumerate(var_to_plot):
    fig,ax=plt.subplots(1,1,figsize=(10,10))
    
    for j,input in enumerate(inputs):
        plt.gca().set_prop_cycle(None)
        data=input.pivot_table(index=['phi'],columns=['EGR','FB'],values=var)
        print(data)
        
        label = [str(data.columns.names[0])+':'+str(val[0])+' '+
                 #str(data.columns.names[1])+ 
                 #X,'H2' as index and 'fuel' as exponent
                 r'$\ X_{H2}^{fuel}$'+
                 ':'+str(val[1]) for val in data.columns.values]
        d = data.plot(#x='phi',y=var,
                  ax=ax,style=symbols[j],
                  label=label,#color=colors[i*j]
                 )
        ax.legend(label,loc='upper left',handler_map={d: HandlerLine2D(numpoints=0)})
        
    plt.xlabel('Equivalence ratio')
    plt.ylabel(ylabels[i])
    plt.title(titles[i]+' - ('+r"$\bf{"+'T_{in}CO2:'+str(input['Tin'][0])+'K'+ "}$"+')')

    #add a second legend to show each marker for each mech
    ax2 = ax.twinx()
    for k,_ in enumerate(inputs):
        ax2.plot([],[],symbols[k],
                 label=mech[k], 
                 c='black',
                )
    ax2.get_yaxis().set_visible(False)
    ax2.legend(loc='upper right')
    ax.grid()
    
    #plt.legend(label,loc='upper left')
    plt.savefig(path+'/img/'+
                'ARC2'+
                '_'+var+'.png')

    
