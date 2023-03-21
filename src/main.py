from lib_egr import *
import time
import pandas as pd

if __name__ == '__main__':
    # get the start time
    st = time.time()

    config = case('CH4:1.',                     #fuel compo
                  'O2:1. N2:3.76',              #ox compo
                  'CO2:1.',                     #egr compo
                  [0.8,0.9,1.0,1.1,1.2],         #phi range
                  [0.0,0.1],                     #egr range
                  'mole'                         #egr rate unit
                 )
    dfs=[]
    pressures = [100000,500000] #Pa
    for p in pressures:
        #set reservoirs thermo-state
        config.res.fuel = create_reservoir(config.compo.fuel,'gri30.cti', 300.0, p)
        config.res.ox = create_reservoir(config.compo.ox,'air.xml', 300.0, p)
        config.res.egr = create_reservoir(config.compo.egr,'gri30.cti', 300.0, p)

        reactor,pdresult = compute_solutions_0D(config,real_egr=False)
        
        dfs.append(pdresult)
    
    print(dfs)

    # get the end time
    et = time.time()
    # get the execution time
    elapsed_time = et - st

        
    #select the columns to plot (X,Xbis,Y)    
    dfsc=[df.pivot_table(index='phi',columns='EGR',values='T') for df in dfs]
    dfsc=pd.concat(dfsc, axis = 1, keys = pressures)
    
    print(dfsc)
    print('Execution time:', elapsed_time, 'seconds')

    title='(0D) Flame temperature vs equivalence ratio (Tin_EGR:'+str(config.res.egr.thermo.T)+'K)'
    human_labels = [str(round(p/100000,1))+' bar, '+str(round(e*100,1))+"%EGR"+config.egr_unit for p in pressures for e in config.egr_range]
    xlabel='Equivalence ratio'
    ylabel='Tad [K]'

    show_graphs(dfsc,title,human_labels,xlabel,ylabel)