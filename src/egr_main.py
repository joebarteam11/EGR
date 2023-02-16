from egr import *
import time
from egr_multiproc import *

fig, ax = plt.subplots(1,1)

if __name__ == '__main__':
    # get the start time
    st = time.time()

    config = case('CH4:1.',                   #fuel compo
                'O2:1. N2:3.76',  #N2:3.76            #ox compo
                'CO2:0.5',                    #egr compo
                3960.0,                       #thermal output power NOT IMPLEMENTED YET
                0.0001,                       #egr rate
                'vol'                         #egr rate unit
                )

    #set reservoirs thermo-state
    config.res.fuel = create_reservoir(config.compo.fuel,'gri30.xml', 300.0, 100000)
    config.res.ox = create_reservoir(config.compo.ox,'air.xml', 300.0, 100000)
    config.res.egr = create_reservoir(config.compo.egr,'gri30.xml', 350.0, 100000)

    #range of computation
    egr_percentages = np.arange(0.0,0.11,0.05)
    df = pd.DataFrame()
    for egr_rate in egr_percentages:
        config.egr_rate = egr_rate #override config.egr_rate set during object instanciation
        phi_range = np.arange(0.6,1.95,0.05)
        reactor,results,pdresult = compute_solutions(config,phi_range,power_regulation=False)
        #plt.plot(phi_range,results[:,4],label='EGR reacteurs:'+str(round(config.egr_rate*100,1)),marker='o')
        #subplot_data(phi_range,results,'Phi',['T[K]','HRR[W/m3]','Y_O2','Y_CO2'],'EGR rate (%):'+str(round(config.egr_rate*100,1))+'%')
        df=pd.concat([df,pdresult])
    
    df = df.pivot_table(index='phi',columns='EGR',values='T')
    print_reactor(df)

    """config.pow = 10000
    for egr_rate in egr_percentages:
        config.egr_rate = egr_rate #override config.egr_rate set during object instanciation
        phi_range = np.arange(0.6,1.95,0.05)
        reactor,results = compute_solutions(config,phi_range,power_regulation=True)
        #print(reactor.volume)
        plt.plot(phi_range,results[:,4],label='EGR reacteurs:'+str(round(config.egr_rate*100,1)),marker='o')
        #subplot_data(phi_range,results,'Phi',['T[K]','HRR[W/m3]','Y_O2','Y_CO2'],'EGR rate (%):'+str(round(config.egr_rate*100,1))+'%',symbol='--o')
        #print_reactor(reactor)
    """
    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')
    df.plot(ax=ax, style='-o',title='Temperature vs equivalence ratio',xlabel='Equivalence ratio',ylabel='T',legend=False, label='EGR Reactor')

    equilibrate_data=main()
    equilibrate_data=equilibrate_data.pivot_table(columns='EGR',index='phi',values='T')
    equilibrate_data.plot(ax=ax, style='--x',title='Temperature vs equivalence ratio',xlabel='Equivalence ratio',ylabel='T',legend=False, label='EGR Reactor')
    plt.grid()
    plt.legend(loc='best')
    plt.show()
    #see_graphs('Mode:'+config.mode
    #          +' | fuel:'+str(round(config.res.fuel.thermo.P,1)) +' Pa, '+ str(round(config.res.fuel.thermo.T,1))+' K'
    #          +' | ox:'+str(round(config.res.ox.thermo.P,1)) +' Pa, '+ str(round(config.res.ox.thermo.T,1))+' K'
    #          +' | egr:'+str(round(config.res.egr.thermo.P,1)) +' Pa, '+ str(round(config.res.egr.thermo.T,1))+' K'
    #          )