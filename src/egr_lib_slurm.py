import sys,os
import time
import numpy as np
import cantera as ct
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing as mp
from packaging import version
#from alive_progress import alive_bar
from tqdm.dask import TqdmCallback
from dask_jobqueue import SLURMCluster
from dask_mpi import initialize
#initialize()
from dask.distributed import Client,progress
#from dask.diagnostics import ProgressBar
from lib_egr_260 import *

sys.path.append(os.getcwd())
path = os.getcwd()

print('Time :',time.strftime("%Y%m%d-%H%M%S"))
print('Current folder: ',path)
print(f"Running Cantera version: {ct.__version__}")
#print(f"Running Python version: {sys.version}")
print(f"Running Matplotlib version: {matplotlib.__version__}")

# get the start time
st = time.time()

config = case('CH4:1.',                     #fuel compo
                [300],                    #tin fuel
                [1e5,5e5],                        #pin fuel
                'O2:1. N2:3.76',              #ox compo
                [300],                    #tin ox
                [1e5,5e5],                        #pin ox
                'CO2:1.',                     #egr compo
                [300],                    #tin egr
                [1e5,5e5],                        #pin egr
                [0.8],        #phi range
                [0.0],            #egr range
                'mole',                       #egr rate unit
                'gri30.cti'                   #scheme
                )

Tins = [[t[i] for t in config.tin] for i in range(len(config.tin.fuel))]
Pins = [[p[i] for p in config.pin] for i in range(len(config.pin.fuel))]

items = [(config,phi,egr,Tin,Pin) for phi in config.phi_range for egr in config.egr_range for Tin in Tins for Pin in Pins]
#print(items)

cluster = SLURMCluster( header_skip=['--mem','-n','--cpus-per-task'],
                    queue='debug',
                    processes=len(items),
                    #n_workers=36,
                    cores=1,
                    memory='80GB',
                    walltime='00:10:00',
                    job_extra=['--nodes=1',
                               '--ntasks-per-node=36',
                               ],
                    )
cluster.scale(jobs=1) #the number of nodes to request
print(cluster.job_script())

ncpu = mp.cpu_count()
print('nCPU :',ncpu)
print('nCase :',len(items))
print([i for i in zip(items)])

solve=True

if(solve):
    
    dask_client = Client(cluster)
    #dask_client.register_worker_callbacks(update)
    #with TqdmCallback():
    #lazy_results=[]
    #for item in items:
    res = dask_client.map(compute_solutions_1D, *zip(*items))#, callback=update)
    progress(res)
    #    lazy_results.append(res)
    results = dask_client.gather(res)

    

    data = pd.concat(results,axis=0)

    print('Closing client...')
    #dask_client.close()

    print('Saving data...')
    data.to_csv(path+'/plan_complet_'+time.strftime("%Y%m%d-%H%M%S")+'.csv',index=False)

    print(data)
