import time
from mpi_func import *



if __name__ == '__main__':

    try: 
        from mpi4py import MPI
        mpiprint("mpi4py properly installed, // available ",priority="info",file=sys.stdout)
        MPI_LOADED=True
    except:
        MPI_LOADED=False
        mpiprint("mpi4py not installed, can only run on one CPU",priority="warning",file=sys.stdout)
        pass

    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

    # MPI Init.
    path = os.path.dirname(os.path.dirname(__file__))+'/'
    if MPI_LOADED:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD 
        ncpu,myrank,rank_0=initialize_MPI(comm)
    else:
        ncpu=1
        myrank=0
        rank_0=True

    ct.suppress_thermo_warnings() 

    if rank_0:
        
        mpiprint('Current folder: '+path,file=sys.stdout)
        mpiprint(f"Running Cantera version: {ct.__version__}",file=sys.stdout)
        mpiprint(f"Running on {ncpu} cores",file=sys.stdout)
        # get the start time
        st = time.time()

        temptlist = [300,305,310,320] #[i for i in np.arange(290,305,100.0)]
        presslist= [1E5] #[i for i in np.arange(1E5,1.4E5,0.2E5)]
        phirange = [0.85] #[0.7,0.8,0.9,1.0,1.05,1.1,1.2]#[i for i in np.arange(0.705,1.305,0.100)] 
        fuelblendrange = [i for i in np.arange(0.0,0.301,0.100)] # [0]#
        egrrange = [i for i in np.arange(0.0,0.301,0.1)]
        config = case(['CH4:1.0','H2:1.0'],         #fuel compo
                    temptlist,                    #tin fuel
                    presslist,                    #pin fuel
                    'O2:1.0 N2:3.76',              #ox compo
                    temptlist,                    #tin ox
                    presslist,                    #pin ox
                    'CO2:1.0',                     #egr compo
                    temptlist,                    #tin egr
                    presslist,                    #pin egr
                    phirange,                     #[i for i in np.arange(0.60,2.51,0.1)],        #phi range
                    fuelblendrange,
                    egrrange,                     #[i for i in np.linspace(0.0,0.6,30)],#[0.0,0.1,0.15,0.2],            #egr range
                    'mole',                       #fuelblend rate unit mole / mass
                    'mole',                       #egr rate unit mole / mass
                    path+'schemes/CH4_15_256_9_AP.cti', #mechanism file
                    'AVBP', #transport model
                    'ARC',  #is an ARC chemistry ? 'ARC' = yes, other = no
                    )
        

        Tins = [[t[i] for t in config.tin] for i in range(len(config.tin.fuel))]
        Pins = [[p[i] for p in config.pin] for i in range(len(config.pin.fuel))]
        
        items = [[config,phi,Tin,Pin,EGR,FB] for phi in config.phi_range for EGR in config.egr_range for FB in config.fuelblend_range for Tin in Tins for Pin in Pins ]
        
        init_cases(config)
    else:
        items=[]

    dim='1D'
    time_formated = time.strftime("%Y%m%d-%H%M%S")
    optimise_mpi_flame_order = False
    species_bg_output = ['O2','CO','CO2']
    real_egr = False
    dry=True
    T_reinj=300

    restart_rate = None
    if (real_egr):
        save_file_name = path + "/results/" + dim + "REAL_EGR" + ".csv"
    else:
        save_file_name = path + "/results/" + dim + "NO_EGR" + ".csv"

    if ncpu==1:
        if MPI_LOADED:
            PRINT_MONO_CPU_WARNING()
        else: 
            PRINT_MONO_CPU_WARNING_AND_MPI_NOT_LOADED()
            for i in range(50):
                mpiprint('\n WARNING, mpi4py is not installed in you env. Please do not use mpirun -n ...',file=sys.stdout)
                

        MONO_CPU_CALCULATION(items,save_file_name,dim,restart_rate,real_egr,dry,T_reinj,species_bg_output)
    else:
        MPI_CALCULATION(rank_0,items,comm,ncpu,optimise_mpi_flame_order,save_file_name,dim,restart_rate,real_egr,dry,T_reinj,species_bg_output)

    if MPI_LOADED:
        comm.Barrier()
    # get the execution time
    if rank_0:
        et = time.time()
        elapsed_time = et - st

        mpiprint('\n Execution time: '+ str(elapsed_time)+' seconds',file=sys.stdout)
