from lib_egr_260 import * #(compute_equilibrium,compute_solutions_0D,compute_solutions_1D)
import time
import pandas as pd
from mpi_func import *


if __name__ == '__main__':

    # MPI Init.
    path = os.getcwd()
    comm,ncpu,myrank,rank_0=initialize_MPI()

    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
    
    if rank_0:
        
        mpiprint('Current folder: '+path)
        mpiprint(f"Running Cantera version: {ct.__version__}")
        # get the start time
        st = time.time()

        temptlist = [i for i in np.arange(290,305,100.0)]
        presslist= [i for i in np.arange(1E5,1.4E5,0.2E5)]
        phirange = [i for i in np.arange(0.7,1.3,0.2)]
        fuelblendrange = [i for i in np.arange(0.0,0.3,0.1)]
        egrrange = [i for i in np.arange(0.0,0.1,0.2)]
        config = case(['CH4:1.','H2:1'],                     #fuel compo
                    temptlist,                    #tin fuel
                    presslist,                        #pin fuel
                    'O2:1. N2:3.76',              #ox compo
                    temptlist,                    #tin ox
                    presslist,                        #pin ox
                    'CO2:1.',                     #egr compo
                    temptlist,                    #tin egr
                    presslist,                        #pin egr
                    phirange,   #[i for i in np.arange(0.60,2.51,0.1)],        #phi range
                    fuelblendrange,
                    egrrange,   #[i for i in np.linspace(0.0,0.6,30)],#[0.0,0.1,0.15,0.2],            #egr range
                    'mole',                       #fuelblend rate unit mole / mass
                    'mole',                       #egr rate unit mole / mass
                    'schemes/CH4_15_256_9_AP.cti',               #path to scheme
                    'AVBP', #transport model
                    'ARC',  #is an ARC chemistry ? 'ARC' = yes, other = no
                    )
        

        Tins = [[t[i] for t in config.tin] for i in range(len(config.tin.fuel))]
        Pins = [[p[i] for p in config.pin] for i in range(len(config.pin.fuel))]
        
        items = [[config,phi,Tin,Pin,EGR,FB] for phi in config.phi_range for EGR in config.egr_range for FB in config.fuelblend_range for Tin in Tins for Pin in Pins ]
        #print(items)

        # Get start time for intermediate result saving
        

        mpiprint('nCPU :'+str(ncpu))
        mpiprint('nItems :'+str(len(items)))

        #Progress bar declaration
        # pbar=tqdm(total=len(items),file=sys.stdout) 
        # def update(*a):
        #     pbar.update()
    else:
        items=[]


    dim='1D'
    time_formated = time.strftime("%Y%m%d-%H%M%S")
    optimise_mpi_flame_order = True
    species = ['O2','CO','CO2']
    real_egr = False
    restart_rate = None
    save_file_name = path + "/" + dim + "_" + time_formated + ".csv"


    if ncpu==1:
        PRINT_MONO_CPU_WARNING()
        MONO_CPU_CALCULATION(items,species,save_file_name,dim,real_egr,restart_rate)
    else:
        MPI_CALCULATION(rank_0,items,comm,ncpu,optimise_mpi_flame_order,save_file_name,species,dim,restart_rate,real_egr)

    # get the execution time
    if rank_0:
        et = time.time()
        elapsed_time = et - st

        print('Execution time:', elapsed_time, 'seconds')
        sys.stdout.flush()
        MPI.COMM_WORLD.Abort(0)


