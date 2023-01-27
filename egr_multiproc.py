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
egr = {'CO2':1.0}               # EGR composition

phi_range = np.arange(0.6,2.1,0.05)

def task(phi_range, egr_rate, scheme='gri30.xml'):
    gas = ct.Solution(scheme)
    data=[]
    df = pd.DataFrame(['egr','phi','T'])
    for phi in phi_range:
        gas.TP = 300.0, 100000.0
        gas.set_equivalence_ratio(phi, fuel, 
                                    oxidizer,
                                    basis="mole",diluent=egr,fraction={"diluent":egr_rate})
        print(phi,gas['CH4'].X,gas['O2'].X,gas['N2'].X,gas['CO2'].X)
        gas.equilibrate('HP')
        
        data.append([egr_rate,phi,gas.T ]) #gas['O2'].Y[0]
    return data

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

def plot(df,fig,ax):
    df.plot(ax=ax, style='-o',title='Temperature vs equivalence ratio',xlabel='Equivalence ratio',ylabel='Temperature (K)',legend=False)
    fig.tight_layout()
    plt.grid()
    plt.show()


if __name__ == '__main__':
    #get start time
    st=time.time()

    fig, ax = plt.subplots(1, 1)
    data=main()

    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')
    plot(data,fig,ax)

