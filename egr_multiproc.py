import numpy as np
import matplotlib.pyplot as plt
from matplotlib import *
import cantera as ct
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
import pandas as pd
import time

# Use reaction mechanism GRI-Mech 3.0. For 0-D simulations,
# no transport model is necessary.

fuel = {'CH4': 1.0}         # Methane composition
oxidizer = {'O2': 1.0, 'N2':3.76} #, 'N2':3.76
egr = {'CO2':0.5}               # EGR composition

phi_range = np.arange(0.7,1.35,0.05)

def task(phi_range, egr_rate,p, scheme='gri30.cti'):
    gas = ct.Solution(scheme)
    data=[]
    df = pd.DataFrame(['egr','phi','T'])
    for phi in phi_range:
        gas.TP = 300.0, p
        gas.set_equivalence_ratio(phi, fuel, oxidizer, basis="mole",
                                diluent=egr,fraction={"diluent":egr_rate})
        #print(phi,gas['CH4'].X,gas['O2'].X,gas['N2'].X,gas['CO2'].X)
        #print(gas['CO2'].X*gas.volume_mole)
        corr_phi = float(2.0*(gas['CH4'].X+gas['H2'].X)/gas['O2'].X)
        #print('corrected_phi : ',corr_phi)
        gas.equilibrate('HP')
        data.append([egr_rate,corr_phi,gas.T, gas['CH4'].X, gas['H2'].X, gas['O2'].X, gas['CO2'].X, gas['H2O'].X]) #gas['O2'].Y[0]
    return data

def main(pressures):
    
    dfs=[]
    
    for p in pressures:
        # start the process pool
        with ProcessPoolExecutor(max_workers=4)as executor:
            # submit many tasks
            egr_range = [0.0,0.1] #np.arange(0.0,0.11,0.1)
            futures = [executor.submit(task,phi_range,egr,p) for egr in egr_range]
            print('Waiting for tasks to complete...')
            # update each time a task finishes
            for _ in as_completed(futures):
                # report the number of remaining tasks
                print(f'About {len(executor._pending_work_items)} tasks remain')
        print('Done')
        df = pd.DataFrame([{'EGR':item[0], 'phi':item[1], 'T':item[2], 'CH4':item[3], 'H2':item[4], 'O2':item[5], 'CO2':item[6], 'H2O':item[7]} for future in futures for item in future.result()]).astype(float)
        #df=df.pivot_table(columns='EGR',index='phi',values='T')
        dfs.append(df)
    #dfsc=pd.concat(dfs, axis = 1, keys = pressures)
    #print(dfsc)
    return dfs

def plot(df,fig,ax,pressures):

    egr_percentages = [0.0,0.1]
    human_labels = [str(round(p/100000,1))+' bar, '+str(round(e*100,1))+"%EGRmol" for p in pressures for e in egr_percentages]
    df.plot(ax=ax, style='x',title='Temperature vs equivalence ratio',xlabel='Equivalence ratio',ylabel='T',legend=False, label='EGR equilibrate')
    fig.tight_layout()
    plt.grid()
    plt.legend()
    ax.legend(human_labels)
    plt.show()


if __name__ == '__main__':
    #get start time
    st=time.time()

    fig, ax = plt.subplots(1, 1)

    pressures = [100000,500000]
    data=main(pressures)
    data=[df.pivot_table(index='phi',columns='EGR',values='T') for df in data]
    datac=pd.concat(data, axis = 1)
    print(datac)
    #data.plot(ax=ax, style='--x',title='Temperature vs equivalence ratio',xlabel='Equivalence ratio',ylabel='T',legend=False, label='EGR Reactor')


    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')
    plot(datac,fig,ax,pressures)

