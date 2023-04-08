import sys,os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib 

sys.path.append(os.getcwd())
path = os.getcwd()

from lib_egr_260 import show_graphs
print('Current folder: ',path)
print(f"Running Matplotlib version: {matplotlib.__version__}")

files=[
    'plan_partiel_0D_dilution_20230408-194012.csv'
    # 'plan_partiel_dilution_0.0_20230406-120038.csv',
    # 'plan_partiel_dilution_0.1_20230406-131614.csv',
    # 'plan_partiel_dilution_0.3_20230406-201540.csv',
    # 'plan_partiel_dilution_0.5_20230407-121510.csv',
]
files=[path+'/results/'+file for file in files]

inputs=pd.concat([pd.read_csv(file).round({'P':1,'EGR':1,'phi':2}) for file in files])
papier=pd.read_csv(path+'/results/'+'data-papier.csv',delimiter=';').round({'P':1,'EGR':1,'phi':2})

#print(input)

var_to_plot=['T',
            # 'u'
            ]
ylabels = ('Tad[K]',
           # 'u[m.s-1]',
          )
titles = ['Flame temperature',
          #'Laminar flame velocity',
         ]


for idx,var in enumerate(var_to_plot):

    
    for i in range(2):
        mydata=inputs.pivot_table(index='phi',columns=['P','EGR',],values=var)
        paper=papier.pivot_table(index='phi',columns=['P','EGR',],values=var)

        #get only P > 200000 in df2 and df
        if(i==0):
            mydata=mydata.loc[:,mydata.columns.get_level_values('P')>200000]
            paper=paper.loc[:,paper.columns.get_level_values('P')>200000]

            #get only EGR < 0.7 in df3 and df1
            mydata=mydata.loc[:,mydata.columns.get_level_values('EGR')<0.7]
            paper=paper.loc[:,paper.columns.get_level_values('EGR')<0.7]
        else:
            mydata=mydata.loc[:,mydata.columns.get_level_values('P')<200000]
            paper=paper.loc[:,paper.columns.get_level_values('P')<200000]

            #get only EGR < 0.7 in df3 and df1
            mydata=mydata.loc[:,mydata.columns.get_level_values('EGR')<0.7]
            paper=paper.loc[:,paper.columns.get_level_values('EGR')<0.7]

        #create a string label for each P and EGR in df3 with format '%EGR: ; P:'
        #Create a string label based on df3 columns and names
        labels = ['%CO2:{}%;{}:{}bar'.format(round(x[1]*100,0), mydata.columns.names[1] ,round(x[0]/100000,0)) for x in mydata.columns]

        
        
        #print(df3.columns.names)

        title='(0D) '+titles[idx]+' vs equivalence ratio (Tin_CO2:'+str(300)+'K)'
        human_labels = labels
        xlabel='Phi'
        ylabel=ylabels[idx]

        _,ax = plt.subplots(1,1,figsize=(10,10))
        mydata.plot(ax=ax,style='.-',title=title,legend=False)
        paper.plot(ax=ax,style='x',title=title,legend=False)
        plt.grid()
        ax.legend(human_labels,loc='best')
        plt.tight_layout()
        #plt.show()
        plt.savefig(path+'/img/'+var+str(i)+'_0D.png', dpi=300, bbox_inches='tight')
        #plt.close()

    #show_graphs(mydata,title,human_labels,xlabel,ylabel,subplot=1,plot=False,save=False,path=path+'/img/')
    #show_graphs(paper,title,human_labels,xlabel,ylabel,subplot=1,ax=ax,save=True,path=path+'/img/',style='x')
    
