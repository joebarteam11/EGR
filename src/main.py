import time
from mpi_func import *



if __name__ == '__main__':

    try: 
        from mpi4py import MPI
        # mpiprint("mpi4py properly installed, // available ",priority="info")
        MPI_LOADED=True
    except:
        MPI_LOADED=False
        mpiprint("mpi4py not installed, can only run on one CPU",priority="warning",file=sys.stdout)
        pass

    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

    # MPI Init.
    path = os.path.dirname(os.path.dirname(__file__))+'/'
    if MPI_LOADED:
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

        templistOx = [393] #[i for i in np.arange(290,305,10.0)]
        templistFuel = [393] #[i for i in np.arange(290,305,10.0)]
        templistEGR = templistOx
        presslist= [5.01325E5] #[i for i in np.arange(1E5,1.4E5,0.2E5)]
        phirange =  [i for i in np.arange(0.605,1.206,0.1)] # [0.85,0.1] # [0.6] [0.6,0.7,0.8,0.9,1.0,1.05,1.1005,1.2005,1.3005,1.4005]#
        # fuelblendrange = [0.0] if no fuel blend needed
        fuelblendrange = [0.0]#[i for i in np.arange(0.0,0.301,0.100)] # 
        # egrrange = [0.0] if no dilution needed
        egrrange = [0.2]#[i for i in np.arange(0.0,0.301,0.1)] 
        config = case(['CH4:1.0','H2:1.0'],       #fuel compo   e.g. ['CH4:1.0','H2:1.0'] with fuel blend = [0.0] is equivalent to pure CH4 as fuel, with fuel blend = [0.1] is equivalent to 90% CH4 and 10% H2
                    templistFuel,                 #tin fuel
                    presslist,                    #pin fuel
                    'O2:1.0 N2:3.76',             #ox compo
                    templistOx,                   #tin ox
                    presslist,                    #pin ox
                    'CO2:1.0',                    #egr compo
                    templistEGR,                  #tin egr
                    presslist,                    #pin egr
                    phirange,                     #[i for i in np.arange(0.60,2.51,0.1)],        #phi range
                    fuelblendrange,               #fuelblend 
                    egrrange,                     #[i for i in np.linspace(0.0,0.6,30)],#[0.0,0.1,0.15,0.2],            #egr range
                    'mole',                       #fuelblend rate unit mole / mass
                    'mole',                       #egr rate unit mole / mass
                    path+'schemes/CH4_15_256_9_AP.cti',  #fuel mechanism
                    #'gri30.cti',
                    'Mix', #transport model
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
    species_bg_output = ['CH4','O2','CO','CO2','H2O','N2']

    # tolerences for 1D flame solver
    tol_ss = [2.0e-5, 1.0e-9]  # tolerance [rtol atol] for steady-state problem
    tol_ts = [2.0e-5, 1.0e-9]  # tolerance [rtol atol] for time stepping

    real_egr = False
    dry=True
    T_reinj= 573

    restart_rate = None
    if (real_egr):
        save_file_name = path + "/results/" + dim + "M12_N2CO2H2OH2ReacOnly" + ".csv"
    else:
        save_file_name = path + "/results/" + dim + "TEST_optim_flame_resolution" + ".csv"

    if ncpu==1:
        if MPI_LOADED:
            PRINT_MONO_CPU_WARNING()
        else: 
            PRINT_MONO_CPU_WARNING_AND_MPI_NOT_LOADED()
            #for i in range(50):
            mpiprint('\n WARNING, mpi4py is not installed in your env. Please do not use mpirun -n ...',file=sys.stdout)
                

        MONO_CPU_CALCULATION(items,tol_ss,tol_ts,save_file_name,dim,restart_rate,real_egr,dry,T_reinj,species_bg_output)
    else:
        MPI_CALCULATION(rank_0,items,tol_ss,tol_ts,comm,ncpu,optimise_mpi_flame_order,save_file_name,dim,restart_rate,real_egr,dry,T_reinj,species_bg_output)

    if MPI_LOADED:
        comm.Barrier()
    # get the execution time
    if rank_0:
        et = time.time()
        elapsed_time = et - st

        mpiprint('\n Execution time: '+ str(elapsed_time)+' seconds',file=sys.stdout)
