import sys,os
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append(os.getcwd())
path = os.getcwd()

from lib_egr_260 import *
print('Current folder: ',path)

files=[
    'plan_partiel_dilution_0.0_20230406-120038.csv',
    'plan_partiel_dilution_0.1_20230406-131614.csv',
    'plan_partiel_dilution_0.3_20230406-201540.csv',
    'plan_partiel_dilution_0.5_20230407-121510.csv',
]
files=[path+'/results/'+file for file in files]

inputs=pd.concat([pd.read_csv(file).round({'P':1,'EGR':1,'phi':2}) for file in files])
papier=pd.read_csv(path+'/results/'+'data-papier.csv',delimiter=';').round({'P':1,'EGR':1,'phi':2})

#print(input)

var_to_plot=['T','u']
ylabels = ('Tad[K]','u[m.s-1]',)
titles = ['Flame temperature', 'Laminar flame velocity',]


for idx,var in enumerate(var_to_plot):

    mydata=inputs.pivot_table(index='phi',columns=['P','EGR',],values=var)
    paper=papier.pivot_table(index='phi',columns=['P','EGR',],values=var)

    #get only P > 200000 in df2 and df
    mydata=mydata.loc[:,mydata.columns.get_level_values('P')>200000]
    paper=paper.loc[:,paper.columns.get_level_values('P')>200000]

    #get only EGR < 0.7 in df3 and df1
    mydata=mydata.loc[:,mydata.columns.get_level_values('EGR')<0.7]
    paper=paper.loc[:,paper.columns.get_level_values('EGR')<0.7]

    #create a string label for each P and EGR in df3 with format '%EGR: ; P:'
    #Create a string label based on df3 columns and names
    labels = ['%s: %s; P: %s' % (mydata.columns.names[1], x[1], x[0]) for x in mydata.columns]

    
    
    #print(df3.columns.names)

    title='(1D) '+titles[idx]+' vs equivalence ratio (Tin_EGR:'+str(300)+'K)'
    human_labels = labels
    xlabel='Phi'
    ylabel=ylabels[idx]

    ax=show_graphs(mydata,title,human_labels,xlabel,ylabel,subplot=1,save=True,path=path+'/img/')
    show_graphs(paper,title,human_labels,xlabel,ylabel,subplot=1,ax=ax,plot=True,path=path+'/img/',style='x')
    
