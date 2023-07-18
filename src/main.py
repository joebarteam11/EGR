from lib_egr_260 import * #(compute_equilibrium,compute_solutions_0D,compute_solutions_1D)
import time
import pandas as pd
from mpi_func import *



# MPI Init.
comm,ncpu,myrank,rank_0=initialize_MPI()


sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

if ncpu == 1:
    for i in range(10):
        mpiprint('!!! WARNING !!!')
        mpiprint('I am built differently, I only know how to do MPI nproc >= 2 yet, sorry')


if __name__ == '__main__':

    if rank_0:
        path = os.getcwd()
        mpiprint('Current folder: '+path)
        mpiprint(f"Running Cantera version: {ct.__version__}")
        # get the start time
        st = time.time()
        temptlist = [i for i in np.arange(290,305,100.0)]
        presslist= [i for i in np.arange(1E5,1.4E5,0.2E5)]
        phirange = [i for i in np.arange(0.7,1.3,0.2)]
        egrrange = [i for i in np.arange(0.0,0.1,0.2)]
        config = case('CH4:1.',                     #fuel compo
                    temptlist,                    #tin fuel
                    presslist,                        #pin fuel
                    'O2:1. N2:3.76',              #ox compo
                    temptlist,                    #tin ox
                    presslist,                        #pin ox
                    'CO2:1.',                     #egr compo
                    temptlist,                    #tin egr
                    presslist,                        #pin egr
                    phirange,   #[i for i in np.arange(0.60,2.51,0.1)],        #phi range
                    egrrange,   #[i for i in np.linspace(0.0,0.6,30)],#[0.0,0.1,0.15,0.2],            #egr range
                    'mole',                       #egr rate unit mole / mass
                    'schemes/CH4_16_250_10_QC.cti',               #path to scheme
                    'Mix', #transport model
                    'ARC',  #is an ARC chemistry ? 'ARC' = yes, other = no
                    )
        

        Tins = [[t[i] for t in config.tin] for i in range(len(config.tin.fuel))]
        Pins = [[p[i] for p in config.pin] for i in range(len(config.pin.fuel))]
        optimise_mpi_flame_order = True
        items = [[config,phi,Tin,Pin,EGR] for phi in config.phi_range for EGR in config.egr_range for Tin in Tins for Pin in Pins ]
        #print(items)

        # Get start time for intermediate result saving
        time_file = time.strftime("%Y%m%d-%H%M%S")

        mpiprint('nCPU :'+str(ncpu))
        mpiprint('nItems :'+str(len(items)))

        #Progress bar declaration
        # pbar=tqdm(total=len(items),file=sys.stdout) 
        # def update(*a):
        #     pbar.update()
    else:
        items=[]

        
    dim='0D'

    if(dim=='equilibrate'):
        species = ['O2','CO','CO2']
        items = comm.bcast(items,root=0)
        results=[]
        n_calculs=(len(items)-1)//(ncpu)
        mpiprint(n_calculs)
        if myrank!=0:
            results = [compute_equilibrium(*item,species) for item in items[(myrank-1)*n_calculs:myrank*n_calculs]]
        else:
            if ncpu!=1:
                n_calculs_0=len(items)-n_calculs*(ncpu)
                mpiprint(str(n_calculs_0)+" : "+str( items[(ncpu-1)*n_calculs:]))
                results = [compute_equilibrium(*item,species) for item in items[(ncpu-1)*n_calculs:]]
            if ncpu==1:
                results = [compute_equilibrium(*item,species) for item in items]
                n_calculs=len(items)
                
        mpiprint('I am here,'+str(myrank))
        results = pd.concat(results[:],axis=0)

        results=comm.gather(results[:],root=0)

        if myrank==0:
            output=pd.concat(results,axis=0)
            output.to_csv(path+'/results'+'/plan_total_equilibrium'+'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv',index=False)
            print(output)

    elif(dim=='0D'):
        real_egr = False
        not_end = True
        species = ['O2','CO','CO2']
        res_time = 0.1
        dfs = []
        results = []

        if rank_0:
            if ncpu != 1:
                items_and_status, requests, itemtot, nb_of_started_flames, nb_of_finished_flames  = initialize_master_1D_flame(items,comm,ncpu)

                while itemtot!=nb_of_finished_flames:
                    intention,talking_to_cpu=receive_intention_from_slave(requests)

                    if intention=='data':
                        items_and_status = master_intention_is_data(comm,talking_to_cpu,items_and_status,results)
                        
                    if intention=='available':
                        items_and_status = master_intention_is_available(comm,items_and_status,talking_to_cpu,itemtot,nb_of_started_flames)

                    items_and_status, requests, nb_of_started_flames, nb_of_finished_flames = update_requests_and_nb_of_flammes(comm,items_and_status,requests,talking_to_cpu)
                    rank0_update_output_log(nb_of_started_flames,nb_of_finished_flames,itemtot)
                    if nb_of_finished_flames> 0 and  nb_of_finished_flames%10==0 : 
                        output.to_csv(path+'/BFER-For-ANSALDO'+'_'+time_file+'.csv',index=False)
                        print("Partial results has been saved")
                        sys.stdout.flush()


                # end of calc.
                output=pd.concat(results[:],axis=0)
                output.to_csv(path+'/BFER-For-ANSALDO'+'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv',index=False)
                print(output)
                print(items_and_status.to_string())
                sys.stdout.flush()

            else:
                mpiprint('I am built differently, I only know how to do MPI yet, sorry')

        else:
            proc0=int(0)
            while not_end:
                msg_to_send='available'
                sends_msg_to_proc(comm,msg_to_send,proc0,1) # Tells proc0 that available
                items=receive_msg_from_proc(comm,0,2) # receive from proc 0 items
                if items is not None:
                    slave_compute_and_communicate_0D(comm,items,real_egr,species)
                else:
                    not_end=False


    elif(dim=='1D'):
        not_end = True
        real_egr = False
        restart_rate = None # config.egr_range[0] #set to None if want to restart computation from the first egr value in egr_range
        results = []
        time_slower=0
        if rank_0:
            # If calc is in parrallel, then goes through // process for proc 0
            if ncpu != 1:
                items_and_status, requests, itemtot, nb_of_started_flames, nb_of_finished_flames  = initialize_master_1D_flame(items,comm,ncpu) 
                while itemtot!=nb_of_finished_flames:

                    if optimise_mpi_flame_order:
                        items_and_status,time_slower=update_priority(items_and_status,ncpu,nb_of_started_flames,time_slower)            

                    intention,talking_to_cpu=receive_intention_from_slave(requests)

                    if intention=='data':
                        items_and_status = master_intention_is_data(comm,talking_to_cpu,items_and_status,results)
                        
                        if nb_of_finished_flames> 0 and  nb_of_finished_flames%10==0 : 
                            output=pd.concat(results[:],axis=0)
                            output.to_csv(path+'/BFER-For-ANSALDO'+'_'+time_file+'.csv',index=False)
                            mpiprint("Partial results has been saved")
                            output=[]


                    elif intention=='available':
                        items_and_status = master_intention_is_available(comm,items_and_status,talking_to_cpu,itemtot,nb_of_started_flames)


                    items_and_status, requests, nb_of_started_flames, nb_of_finished_flames = update_requests_and_nb_of_flammes(comm,items_and_status,requests,talking_to_cpu)
                    rank0_update_output_log(nb_of_started_flames,nb_of_finished_flames,itemtot)
                    print(items_and_status.to_string())
                    sys.stdout.flush()


                # End of calculation for proc 0 
                output=pd.concat(results[:],axis=0)
                output.to_csv(path+'/BFER-For-ANSALDO'+'_'+time_file+'.csv',index=False)
                print(output)
                print(items_and_status.to_string())
                sys.stdout.flush()


            # If calc is not, then goes through serial process for proc 0
            else:
                mpiprint('I am built differently, I only know how to do MPI yet, sorry')

        else:
            proc0=int(0)
            while not_end:
                msg_to_send='available'
                sends_msg_to_proc(comm,msg_to_send,proc0,1) # Tells proc0 that available
                items=receive_msg_from_proc(comm,0,2) # receive from proc 0 items
                if items is not None:
                    slave_compute_and_communicate_1D(comm,items,restart_rate,real_egr,results=[])
                else:
                    not_end=False
    
    # get the execution time
    if rank_0:
        et = time.time()
        elapsed_time = et - st

        print('Execution time:', elapsed_time, 'seconds')
        sys.stdout.flush()

        MPI.COMM_WORLD.Abort(0)