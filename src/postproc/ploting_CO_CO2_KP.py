import sys,os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib 

sys.path.append(os.getcwd())
path = os.getcwd()


print('Current folder: ',path)
print(f"Running Matplotlib version: {matplotlib.__version__}")

files=[
    #'plan_total_dilution_BFERMIX_20230412-154440.csv',
    'plan_partiel_0D_dilutionKP_20230514-130426.csv',
    # 'plan_partiel_dilution_0.0_20230406-120038.csv',
    # 'plan_partiel_dilution_0.1_20230406-131614.csv',
    # 'plan_partiel_dilution_0.3_20230406-201540.csv',
    # 'plan_partiel_dilution_0.5_20230407-121510.csv',
]
files=[path+'/results/'+file for file in files]

inputs=pd.concat([pd.read_csv(file).round({'P':-1,'EGR':1,'T':-1}) for file in files])
#papier=pd.read_csv(path+'/results/'+'data-handmade.csv',delimiter=';').round({'EGR':1,'phi':2})
#papier=pd.read_csv(path+'/plan_total_dilution_BFERUNITY_20230412-170317.csv').round({'P':1,'EGR':1,'phi':2})
#print(input)

var_to_plot=['O2',
             'CO2',
             'CO',
            ]
ylabels = (  '$X_{O2}$',
             '$X_{CO2}$',
             '$X_{CO}$',
          )
titles = ['Molar fraction',
          'Molar fraction',
          'Molar fraction',
         ]

mech=['Aramco1.3','[1]']
colors=['b','r','g','k','m','y','c']
style=['o-','x-']

for idx,var in enumerate(var_to_plot):
    _,ax = plt.subplots(1,1,figsize=(10,10))
    labels=[]
    for i in range(2):
        mydata=inputs.pivot_table(index='T',columns=['P','EGR',],values=var)
        try:
            paper=papier.pivot_table(index='phi',columns=['EGR',],values=var)
        except:
            pass
            
        #get only P > 200000 in df2 and df
        if(i==0):
            mydata=mydata.loc[:,mydata.columns.get_level_values('P')>15_00000]
            try:
                paper=paper.loc[:,paper.columns.get_level_values('P')>15_00000]
            except:
                pass

            #get only EGR < 0.7 in df3 and df1
            mydata=mydata.loc[:,mydata.columns.get_level_values('EGR')<0.7]
            try:
                paper=paper.loc[:,paper.columns.get_level_values('EGR')<0.7]
            except:
                pass
        else:
            mydata=mydata.loc[:,mydata.columns.get_level_values('P')<15_00000]
            try:
                paper=paper.loc[:,paper.columns.get_level_values('P')<15_00000]
            except:
                pass
            #get only EGR < 0.7 in df3 and df1
            mydata=mydata.loc[:,mydata.columns.get_level_values('EGR')<0.7]
            try:
                paper=paper.loc[:,paper.columns.get_level_values('EGR')<0.7]
            except:
                pass
        #print(mydata.columns)
        #create a string label for each P and EGR in df3 with format '%EGR: ; P:'
        #Create a string label based on df3 columns and names
        labels=[('$X_{CO2in}$'+':{}%,P:{}bar'.format(round(x[1]*100,0),round(x[0]/100000,0))) for x in mydata.columns]

        
        
        #print(df3.columns.names)

        title='(0D) '+titles[idx]+' vs temperature with '+mech[0]+' scheme'
        human_labels = labels
        xlabel='Phi'
        ylabel=ylabels[idx]

        fs=20

        
        mydata.plot(ax=ax,style=style[i],legend=False,color=colors[:len(labels)])
        try:
            paper.plot(ax=ax,style='x',legend=False,color=colors[:len(labels)],markersize=25)
        except:
            pass

    #increase font size of ticks
    ax.tick_params(axis='both', which='major', labelsize=fs)
    #increase font size of labels x and y 
    ax.set_xlabel(xlabel,fontsize=fs)
    ax.set_ylabel(ylabel,fontsize=fs)
    ax.set_title(title,fontsize=20)

    #ax.set_xscale('log')
    #ax.set_yscale('log')

    plt.rcParams.update({'font.size': fs})
    plt.grid()
    ax.legend(human_labels,loc='best',numpoints=1,fontsize=fs)
    #ax.legend('[1]')
    plt.tight_layout()
    #plt.show()
    plt.savefig(path+'/img/'+var+'_0D_KP.png', dpi=300, bbox_inches='tight')
    #plt.close()

    #show_graphs(mydata,title,human_labels,xlabel,ylabel,subplot=1,plot=False,save=False,path=path+'/img/')
    #show_graphs(paper,title,human_labels,xlabel,ylabel,subplot=1,ax=ax,save=True,path=path+'/img/',style='x')
    
