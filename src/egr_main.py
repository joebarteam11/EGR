from egr import *
import time
from egr_multiproc import *
import pandas as pd

fig, ax = plt.subplots(1,1)

if __name__ == '__main__':
    # get the start time
    st = time.time()

    config = case('CH4:1.',                   #fuel compo
                'O2:1. N2:3.76',  #N2:3.76            #ox compo
                'CO2:0.5',                    #egr compo
                3960.0,                       #thermal output power NOT IMPLEMENTED YET
                0.0001,                       #egr rate
                'mol'                         #egr rate unit
                )
    dfs=[]
    pressures = [100000,500000]
    for p in pressures:
        #set reservoirs thermo-state
        config.res.fuel = create_reservoir(config.compo.fuel,'Aramco13.cti', 300.0, p)
        config.res.ox = create_reservoir(config.compo.ox,'air.xml', 300.0, p)
        config.res.egr = create_reservoir(config.compo.egr,'Aramco13.cti', 300.0, p)

        #range of computation
        egr_percentages = np.arange(0.0,0.11,0.05)
        df = pd.DataFrame()
        for egr_rate in egr_percentages:
            config.egr_rate = egr_rate #override config.egr_rate set during object instanciation
            phi_range = np.arange(0.7,1.35,0.05)
            reactor,results,pdresult = compute_solutions_0D(config,phi_range,power_regulation=False)
            #plt.plot(phi_range,results[:,4],label='EGR reacteurs:'+str(round(config.egr_rate*100,1)),marker='o')
            #subplot_data(phi_range,results,'Phi',['T[K]','HRR[W/m3]','Y_O2','Y_CO2'],'EGR rate (%):'+str(round(config.egr_rate*100,1))+'%')
            df=pd.concat([df,pdresult])
        
        df = df.pivot_table(index='phi',columns='EGR',values='T')
        print_reactor(df)

        # get the end time
        et = time.time()

        # get the execution time
        elapsed_time = et - st
        print('Execution time:', elapsed_time, 'seconds')
        dfs.append(df)

    dfsc=pd.concat(dfs, axis = 1, keys = pressures)
    print(dfsc)

    human_labels = [str(round(p/100000,1))+' bar, '+str(round(e*100,1))+"%EGRmol" for p in pressures for e in egr_percentages]
    print(human_labels)
    

    dfsc.plot(ax=ax, style='-o',title='(0D)Flame temperature vs equivalence ratio (Tin_EGR:800K)',xlabel='Equivalence ratio',ylabel='T [K]')

    equilibrate_data=main()
    equilibrate_data=equilibrate_data.pivot_table(columns='EGR',index='phi',values='T')
    equilibrate_data.plot(ax=ax, style='--x',title='Temperature vs equivalence ratio',xlabel='Equivalence ratio',ylabel='T',legend=False, label='EGR Reactor')

    plt.grid()
    plt.legend(loc='best')
    ax.legend(human_labels)
    plt.show()

    #see_graphs('Mode:'+config.mode
    #          +' | fuel:'+str(round(config.res.fuel.thermo.P,1)) +' Pa, '+ str(round(config.res.fuel.thermo.T,1))+' K'
    #          +' | ox:'+str(round(config.res.ox.thermo.P,1)) +' Pa, '+ str(round(config.res.ox.thermo.T,1))+' K'
    #          +' | egr:'+str(round(config.res.egr.thermo.P,1)) +' Pa, '+ str(round(config.res.egr.thermo.T,1))+' K'
    #          )