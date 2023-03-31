import sys,os
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append(os.getcwd())
path = os.getcwd()

from lib_egr_260 import *
print('Current folder: ',path)


input = pd.read_csv('plan_complet_20230331-112013.csv')

var_to_plot=['T','u']
labels = ['Tad[K]','u[m.s-1]',]
titles = ['Flame temperature', 'Laminar flame velocity',]

for idx,var in enumerate(var_to_plot):
    data=input.pivot_table(index='phi',columns='EGR',values=var)

    title='(1D) '+titles[idx]+' vs equivalence ratio (Tin_EGR:'+str(300)+'K)'
    human_labels = labels
    xlabel='Phi'
    ylabel=labels[idx]

    show_graphs(data,title,human_labels,xlabel,ylabel,subplot=1,save=True,path=path+'/img/')
