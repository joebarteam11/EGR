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
    #'plan_total_dilution_BFERMIX_20230412-154440.csv',
    'CRASHTEST_20230608-164516.csv',
    # 'plan_partiel_dilution_0.0_20230406-120038.csv',
    # 'plan_partiel_dilution_0.1_20230406-131614.csv',
    # 'plan_partiel_dilution_0.3_20230406-201540.csv',
    # 'plan_partiel_dilution_0.5_20230407-121510.csv',
]
files=[path+'/results/'+file for file in files]

inputs=pd.concat([pd.read_csv(file).round({'P':1,'EGR':1,'phi':2,'Tin':1}) for file in files])
#papier=pd.read_csv(path+'/results/'+'data-handmade.csv',delimiter=';').round({'EGR':1,'phi':2})
#papier=pd.read_csv(path+'/plan_total_dilution_BFERUNITY_20230412-170317.csv').round({'P':1,'EGR':1,'phi':2})
#print(input)

var_to_plot=['dF',
             'u',
             'T',
            ]
ylabels = ('Flame Thickness [um]',
            'SL0[cm.s-1]',
            'T[K]',
          )
titles = ['Flame thickness',
          'Laminar flame speed',
          'Adiabatic flame temperature',
         ]

mech=['SchemaTest','[1]']
colors=['b','r','g','k','m','c','y']

for idx,var in enumerate(var_to_plot):

    for i in range(2):
        mydata=inputs.pivot_table(index='phi',columns=['Tin','P',],values=var)
        try:
            paper=papier.pivot_table(index='phi',columns=['EGR',],values=var)
        except:
            pass

        #get only P > 200000 in df2 and df
        if(i==0):
            mydata=mydata.loc[:,mydata.columns.get_level_values('Tin')<1000]
            try:
                paper=paper.loc[:,paper.columns.get_level_values('Tin')<1000]
            except:
                pass

            #get only EGR < 0.7 in df3 and df1
            mydata=mydata.loc[:,mydata.columns.get_level_values('P')<200000]
            try:
                paper=paper.loc[:,paper.columns.get_level_values('P')<200000]
            except:
                pass
        else:

            mydata=mydata.loc[:,mydata.columns.get_level_values('Tin')<1000]
            try:
                paper=paper.loc[:,paper.columns.get_level_values('Tin')<1000]
            except:
                pass
            #get only EGR < 0.7 in df3 and df1
            mydata=mydata.loc[:,mydata.columns.get_level_values('P')>200000]
            try:
                paper=paper.loc[:,paper.columns.get_level_values('P')>200000]
            except:
                pass

        #create a string label for each P and EGR in df3 with format '%EGR: ; P:'
        #Create a string label based on df3 columns and names
        labels = ['%CO2:{}%;{}:{}bar {}'.format(round(x[1]*100,0), mydata.columns.names[0] ,round(x[0]/100000,0),mech[0]) for x in mydata.columns]

        
        
        #print(df3.columns.names)

        title='(1D) '+titles[idx]+' vs equivalence ratio ('+r"$\bf{"+'T_{in}CO2:'+str(300)+'K'+ "}$"+')'
        human_labels = labels#+['hand-made calculations (constant Cp)']
        xlabel='Phi'
        ylabel=ylabels[idx]

        fs=20

        _,ax = plt.subplots(1,1,figsize=(10,10))
        mydata.plot(ax=ax,style='o-',legend=False,color=colors[:len(mydata.columns)])
        try:
            paper.plot(ax=ax,style='x',legend=False,color=colors[:len(mydata.columns)],markersize=10)
        except:
            pass

        #increase font size of ticks
        ax.tick_params(axis='both', which='major', labelsize=fs)
        #increase font size of labels x and y 
        ax.set_xlabel(xlabel,fontsize=fs)
        ax.set_ylabel(ylabel,fontsize=fs)
        ax.set_title(title,fontsize=20)

        plt.rcParams.update({'font.size': fs})
        plt.grid()
        ax.legend(human_labels,loc='best',numpoints=1,fontsize=fs)
        #ax.legend('[1]')
        plt.tight_layout()
        #plt.show()
        plt.savefig(path+'/'+var+str(i)+'_GRI30.png', dpi=300, bbox_inches='tight')
        #plt.close()

    #show_graphs(mydata,title,human_labels,xlabel,ylabel,subplot=1,plot=False,save=False,path=path+'/img/')
    #show_graphs(paper,title,human_labels,xlabel,ylabel,subplot=1,ax=ax,save=True,path=path+'/img/',style='x')
    
