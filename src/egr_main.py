#import cantera as ct
from egr import *
import time
from egr_multiproc import *
#import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,1)

def main():
    # start the process pool
    with ProcessPoolExecutor(max_workers=4)as executor:
        # submit many tasks

        egr_range = np.arange(0.0,0.21,0.1)
        futures = [executor.submit(task,phi_range,egr) for egr in egr_range]
        print('Waiting for tasks to complete...')
        # update each time a task finishes
        for _ in as_completed(futures):
            # report the number of remaining tasks
            print(f'About {len(executor._pending_work_items)} tasks remain')
    print('Done')

    df = pd.DataFrame([{'EGR':item[0], 'phi':item[1], 'T':item[2]} for future in futures for item in future.result()])
    df=df.pivot_table(columns='EGR',index='phi',values='T')
    print(df)
    return df

if __name__ == '__main__':
    # get the start time
    st = time.time()

    config = case('CH4:1.',                   #fuel compo
                'O2:1. N2:3.76',  #N2:3.76            #ox compo
                'CO2:1.0',                    #egr compo
                3960.0,                       #thermal output power NOT IMPLEMENTED YET
                0.0001,                       #egr rate
                'vol'                         #egr rate unit
                )

    #set reservoirs thermo-state
    config.res.fuel_gas, config.res.fuel = create_reservoir(config.compo.fuel,'gri30.xml', 300.0, 100000)
    config.res.ox_gas, config.res.ox = create_reservoir(config.compo.ox,'air.xml', 300.0, 100000)
    config.res.egr_gas, config.res.egr = create_reservoir(config.compo.egr,'gri30.xml', 300.0, 100000)

    #range of computation
    egr_percentages = np.arange(0.0,0.21,0.1)
    for egr_rate in egr_percentages:
        config.egr_rate = egr_rate #override config.egr_rate set during object instanciation
        phi_range = np.arange(0.6,2.1,0.05)
        reactor,results = compute_solutions(config,phi_range,power_regulation=False)
        plt.plot(phi_range,results[:,0],label='EGR reacteurs:'+str(round(config.egr_rate*100,1)),marker='x')
        #subplot_data(phi_range,results,'Phi',['T[K]','HRR[W/m3]','Y_O2','Y_CO2'],'EGR rate (%):'+str(round(config.egr_rate*100,1))+'%')
        #print_reactor(reactor)
    '''
    config.pow = 10000
    for egr_rate in egr_percentages:
        config.egr_rate = egr_rate #override config.egr_rate set during object instanciation
        phi_range = np.arange(0.55,2.1,0.05)
        phi_bilger,reactor,results = compute_solutions(config,phi_range,print_report=False,real_egr=True)
        #print(reactor.volume)
        subplot_data(phi_range,results,'Phi',['T[K]','HRR[W/m3]','Y_O2','Y_CO2'],'EGR rate (%):'+str(round(config.egr_rate*100,1))+'%',symbol='--o')
        #print_reactor(reactor)
    '''
    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')

    data=main()
    plot(data,fig,ax)
    #plt.legend()

    #see_graphs('Mode:'+config.mode
    #          +' | fuel:'+str(round(config.res.fuel.thermo.P,1)) +' Pa, '+ str(round(config.res.fuel.thermo.T,1))+' K'
    #          +' | ox:'+str(round(config.res.ox.thermo.P,1)) +' Pa, '+ str(round(config.res.ox.thermo.T,1))+' K'
    #          +' | egr:'+str(round(config.res.egr.thermo.P,1)) +' Pa, '+ str(round(config.res.egr.thermo.T,1))+' K'
    #          )