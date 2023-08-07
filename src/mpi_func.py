import sys,os
import numpy as np
from mpi4py import MPI
import pandas as pd
import time
import logging
import cantera as ct
from lib_egr_260 import case,init_cases
from lib_egr_equilibrium import compute_equilibrium
from lib_egr_0D import compute_solutions_0D
from lib_egr_1D import compute_solutions_1D

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))


mpilog = logging.getLogger('mpi')
hl = logging.FileHandler(filename="0.log",mode='a')
format = logging.Formatter('%(message)s')
hl.setFormatter(format)
mpilog.setLevel(logging.INFO)
mpilog.addHandler(hl)


def mpiprint(message_to_log,priority="info",file=None):
    if file is not None:
        print(message_to_log)
        sys.stdout.flush()

    if priority == "info":
        mpilog.info(message_to_log)
    elif priority == "warning":
        mpilog.warning(message_to_log)

try:
    from tqdm import tqdm
except:
    mpiprint("tqdm not installed")
    pass


def initialize_MPI(comm):
    # Initialise MPI
    ncpu = comm.Get_size() # Number of CPU
    myrank = comm.Get_rank() # Gives the rank "rank" for each CPU
    if myrank==0:
        rank_0=True
    else:
        rank_0=False

    return ncpu,myrank,rank_0 # rank_0 is a True/False that tells if a proc is rank 0

def initialize_master_to_slave_communication(comm,ncpu):
    avmsg=[]                                            # List of open request communication
    for Master_is_talking_to_CPUn in range(1,ncpu):
        avmsg.append(comm.irecv(source=Master_is_talking_to_CPUn,tag=1)) # Start communication with proc Master_is_talking_to_CPUn 
    return avmsg

def receive_intention_from_slave(requests): 
    nproc_and_message=MPI.Request.waitany(requests)  
    message=nproc_and_message[1] # Extract the message from the request
    Master_is_talking_to_CPUn=nproc_and_message[0]+1 # Extract from which CPU the message is comming from the request
    return message,Master_is_talking_to_CPUn

def sends_msg_to_proc(comm,msg,Talking_to_CPU,tag): 
    send_msg_request=comm.isend(msg,dest=Talking_to_CPU,tag=tag) # sends message msg, to Talking-to-CPU with tag tag
    send_msg_request.wait() # Waits that talking to cpu receives data

def receive_msg_from_proc(comm,Talking_to_CPU,tag):
    received_msg_request=comm.irecv(source=Talking_to_CPU,tag=tag) # create a receive request for message received_msg_request, from Talking-to-CPU with tag tag
    received_msg = received_msg_request.wait() # waits that Talking_to_CPU send the message 
    return received_msg

def master_items_and_status_update_new_started_flamme(items_and_status,talking_to_cpu): # Function that finds the index of the next flame to start
    temp=True # Boolean to init. while
    i=0
    while temp:
        if items_and_status.loc[i,'started']==1:
            i +=1
        else:
            index_to_start=i
            temp=False
    msg=items_and_status.loc[index_to_start:index_to_start,'items'] # Fills msg with the next item to start 
    items_and_status.loc[index_to_start:index_to_start,'started']=1 # Declare flame as started
    items_and_status.loc[index_to_start:index_to_start,'by_cpu']=talking_to_cpu # Declare flame as calculated by proc talking_to_cpu
    return msg,items_and_status

def rank0_update_output_log(started,finished,itemtot,file=None):
    message_to_print=str('\n!------------------------------!\n!     '+str(finished)+' finished / '+str(itemtot)+'    !\n!     '+str(started)+' started / '+str(itemtot)+'    !\n!------------------------------!')
    mpiprint(message_to_print,file)

def update_priority(items_and_status,ncpu,nb_of_started_flames,time_slower):
    st=time.time()
    items = items_and_status['items']
    if nb_of_started_flames<1: # For the ncpu first flames, flames are started by priority according to pressure, temp, phi, egr ( max P, max T, max EGR, farther PHI from ONE)
        maxP=1
        maxT=1
        maxEGR=0
        fartherPHI=0

        for i in range(len(items)):
            T=items[i][2][1]
            P=items[i][3][1]
            Phi=items[i][1]
            EGR=items[i][4]
            if P>maxP:
                maxP=P
            if T>maxT:
                maxT=T
            if EGR>maxEGR:
                maxEGR=EGR
            if abs(Phi-1)>fartherPHI:
                fartherPHI=abs(Phi-1)
        
        for i in range(len(items)):
            T=items[i][2][1]
            P=items[i][3][1]
            EGR=items[i][4]
            Phi=items[i][1]
            phi_distance_to_1=abs(Phi-1)
            items_and_status.loc[i,'priority']=(items[i][3][1]/maxP + items[i][4]/(maxEGR+1e-16) + items[i][2][1]/maxT + phi_distance_to_1/(fartherPHI+1e-16) )/3

    elif nb_of_started_flames>ncpu: # For flames > ncpu priority is given to the flame that took the longest time to run and for flames that have similar global parameters
        idx_slower=items_and_status['time'].idxmax()
        new_time_slower=items_and_status.loc[idx_slower,'time']
        mpiprint(str(idx_slower)+"   "+str(items_and_status.loc[idx_slower,'time']))

        if new_time_slower!=time_slower:
            time_slower=new_time_slower
            temperature_slower=items[idx_slower][2][1]
            press_slower=items[idx_slower][3][1]
            phi_slower=items[idx_slower][1]
            egr_slower=items[idx_slower][4]
            dist=np.zeros(len(items))
            for i in range(len(items)):
                T=items[i][2][1]
                P=items[i][3][1]
                Phi=items[i][1]
                EGR=items[i][4]
                dist_T= 1 - T/temperature_slower
                dist_P= 1 - P/press_slower
                dist_phi= 1 - Phi/phi_slower
                dist_egr = 1 - EGR/(egr_slower+1e-16)
                dist[i]= (dist_P**2+dist_T**2+dist_phi**2+dist_egr**2)**(0.5)

            dist_max=max(dist)
            dist_min=min(dist)
            for i in range(len(items)):
                items_and_status.loc[i,'priority']=(dist[i]-dist_max)/(dist_min-dist_max+1e-16)

    items_and_status=items_and_status.sort_values(by="priority",ascending=False)
    items_and_status=items_and_status.reset_index(drop=True)
    et=time.time()
    elapsed_time=et-st
    mpiprint("Optimise priority took: "+str(elapsed_time)+" seconds")
    return items_and_status,time_slower

def initialize_master_1D_flame(items,comm,ncpu):
    itemtot=len(items)
    requests=initialize_master_to_slave_communication(comm,ncpu) # Initialise all request between proc 0 and slaves
    items_and_status = {'items': items,'started': np.zeros(len(items)),'finished': np.zeros(len(items)),'time':(np.zeros(itemtot)),'by_cpu':(np.zeros(itemtot)),'priority':(np.ones(itemtot))} # Create a list that contains all important information about ongoing calculations
    items_and_status = pd.DataFrame(data=items_and_status) # Goes into a dataframe
    nb_of_finished_flames=0 # initialise nb of finished & started flames
    nb_of_started_flames=0 # initialise nb of finished & started flames

    return items_and_status,requests, itemtot, nb_of_started_flames, nb_of_finished_flames

def master_intention_is_data(comm,talking_to_cpu,items_and_status,results):
    new_result = receive_msg_from_proc(comm,talking_to_cpu,3) # Received results from talking to cpu 
    results += new_result # Add new results to the lsit of results 
    runtime = receive_msg_from_proc(comm,talking_to_cpu,4) # Receive runtim from talking to cpu
    items_and_status['time'] = np.where( (items_and_status['by_cpu']==talking_to_cpu) & (items_and_status['started']==1) & (items_and_status['finished']==0) , runtime , items_and_status['time']) # Stores the run time of this flame in items and status
    items_and_status['finished'] = np.where( (items_and_status['by_cpu']==talking_to_cpu) & (items_and_status['started']==1) & (items_and_status['finished']==0) , 1 , items_and_status['finished']) # Declare flame as finished in items and status

    return items_and_status

def master_intention_is_available(comm,items_and_status,talking_to_cpu,itemtot,nb_of_started_flames):
    if itemtot==nb_of_started_flames : # Nore more calculation to do
        msg=None
    else:
        msg,items_and_status=master_items_and_status_update_new_started_flamme(items_and_status,talking_to_cpu) # Extract msg ( that contains an item), and updated items_and_status
    
    sends_msg_to_proc(comm,msg,talking_to_cpu,2) # sends items to proc talking to cpu

    return items_and_status

def update_requests_and_nb_of_flammes(comm,items_and_status,requests,talking_to_cpu):
    nb_of_started_flames=(items_and_status['started'].values == 1).sum()  # Count the nb of started flames
    nb_of_finished_flames=(items_and_status['finished'].values == 1).sum() # Count the nb of finished flames
    requests[talking_to_cpu-1]=comm.irecv(source=talking_to_cpu,tag=1) # Re-create a request for talking-to-cpu to receive future intentions
    return items_and_status, requests, nb_of_started_flames, nb_of_finished_flames

def slave_compute_and_communicate_1D(comm,items,restart_rate,real_egr,dry,T_reinj,results):
    proc0=int(0)
    st = time.time() # Beginning of solving flame
    results += [compute_solutions_1D(*item,restart_rate,real_egr,dry,T_reinj) for item in items.iloc[:1]]
    et = time.time() # End of solving flame
    elapsed_time = et - st 
    mpiprint("Elapsed time for solving this flame: "+str(elapsed_time)+" seconds")
    msg_to_send='data'      # New intention with proc 0 to send data        
    sends_msg_to_proc(comm,msg_to_send,proc0,1) # Tells proc0 that data is available to transfer
    sends_msg_to_proc(comm,results,0,3) # Sends results to proc0 
    sends_msg_to_proc(comm,elapsed_time,0,4) 
    
def slave_compute_and_communicate_0D(comm,items,real_egr,species):
    proc0=int(0)
    # msg_to_print='I am rank '+str(myrank)+' my item is',str(items)
    st = time.time() # Beginning of solving flame
    reactor_and_df = [compute_solutions_0D(*item,real_egr,species) for item in items.iloc[:1]]
    results = reactor_and_df[0][1]
    et = time.time() # End of solving flame
    elapsed_time = et - st
    msg_to_send='data'             
    sends_msg_to_proc(comm,msg_to_send,proc0,1) # Tells proc0 that data is available to transfer
    sends_msg_to_proc(comm,[results],0,3) # Sends results to proc0 
    sends_msg_to_proc(comm,elapsed_time,0,4) 

def slave_compute_and_communicate_equilibrate(comm,items,species,results):
    proc0=int(0)
    # msg_to_print='I am rank '+str(myrank)+' my item is',str(items)
    st = time.time() # Beginning of solving flame
    results += [compute_equilibrium(*item,species) for item in items.iloc[:1]]
    et = time.time() # End of solving flame
    elapsed_time = et - st
    msg_to_send='data'             
    sends_msg_to_proc(comm,msg_to_send,proc0,1) # Tells proc0 that data is available to transfer
    sends_msg_to_proc(comm,results,0,3) # Sends results to proc0 
    sends_msg_to_proc(comm,elapsed_time,0,4) 

def MPI_CALCULATION_MASTER(items,comm,ncpu,optimise_mpi_flame_order,save_file_name):
    time_slower = 0
    results = []
    pbar_started=tqdm(total=len(items),file=sys.stdout) # Progress bar declaration
    pbar_started.set_description("Started")
    pbar_finished=tqdm(total=len(items),file=sys.stdout) # Progress bar declaration
    pbar_finished.set_description("Finished")

    items_and_status, requests, itemtot, nb_of_started_flames, nb_of_finished_flames  = initialize_master_1D_flame(items,comm,ncpu)

    while True:  # While calculation is not finished : 
        intention,talking_to_cpu=receive_intention_from_slave(requests) # Receive intention from slaves

        if optimise_mpi_flame_order: # if use decided, 'optimise' calculation order
            items_and_status,time_slower=update_priority(items_and_status,ncpu,nb_of_started_flames,time_slower)     

        if intention=='data': # If slave has data to send to master,      
            items_and_status = master_intention_is_data(comm,talking_to_cpu,items_and_status,results) # Receive this data, 
            try:
                pbar_finished.update() # Update progress bar
                rank0_update_output_log(nb_of_started_flames,nb_of_finished_flames,itemtot) # log calculation status
            except:
                rank0_update_output_log(nb_of_started_flames,nb_of_finished_flames,itemtot,file=sys.stdout) # Prints calculation status
                #pass

        elif intention=='available': # If slave is ready for a calculation,
            items_and_status = master_intention_is_available(comm,items_and_status,talking_to_cpu,itemtot,nb_of_started_flames) # Send next calculation to do            
            try:
                if(nb_of_started_flames != pbar_started.n):
                    pbar_started.update() # Update progress bar
                rank0_update_output_log(nb_of_started_flames,nb_of_finished_flames,itemtot) # log calculation status
            except:
                rank0_update_output_log(nb_of_started_flames,nb_of_finished_flames,itemtot,file=sys.stdout) # Prints calculation status
                #pass

        else :
            rank0_update_output_log(nb_of_started_flames,nb_of_finished_flames,itemtot) # log calculation status

        if nb_of_finished_flames> 0 and  nb_of_finished_flames%10==0 : # Every 10 itération, saves a intermediate result file
            output=pd.concat(results[:],axis=0) 
            output.to_csv(save_file_name,index=False)
            mpiprint("Partial results has been saved")
            mpiprint(items_and_status.to_string())

        if(itemtot==nb_of_finished_flames): #do while loop emulation
            break

        items_and_status, requests, nb_of_started_flames, nb_of_finished_flames = update_requests_and_nb_of_flammes(comm,items_and_status,requests,talking_to_cpu) # Update calculations status

    # When calculation is finished, 
    output=pd.concat(results[:],axis=0) 
    output.to_csv(save_file_name,index=False) 
    pbar_started.close()
    pbar_finished.close()
    mpiprint(output, file=sys.stdout)
    mpiprint(items_and_status.to_string())
 
def MPI_CALCULATION_SLAVE(comm,species,dim,restart_rate,real_egr,dry,T_reinj):
    proc0=int(0)
    not_end = True
    while not_end:
        msg_to_send='available'
        sends_msg_to_proc(comm,msg_to_send,proc0,1) # Tells proc0 that available
        items=receive_msg_from_proc(comm,0,2) # receive from proc 0 items
        if items is not None:
            if dim=='equilibrate':
                slave_compute_and_communicate_equilibrate(comm,items,species,results=[])
            elif dim=='0D':
                slave_compute_and_communicate_0D(comm,items,real_egr,species)
            elif dim=='1D':
                slave_compute_and_communicate_1D(comm,items,restart_rate,real_egr,dry,T_reinj,results=[])
        else:
            not_end=False

def MPI_CALCULATION(rank_0,items,comm,ncpu,optimise_mpi_flame_order,save_file_name,species,dim,restart_rate,real_egr,dry,T_reinj):

    if rank_0:
        MPI_CALCULATION_MASTER(items,comm,ncpu,optimise_mpi_flame_order,save_file_name)
    else:
        MPI_CALCULATION_SLAVE(comm,species,dim,restart_rate,real_egr,dry,T_reinj)

def MONO_CPU_CALCULATION(items,species,save_file_name,dim,restart_rate,real_egr,dry,T_reinj):
    results = []
    if dim=='equilibrate':
        results = [compute_equilibrium(*item,species) for item in items]
    elif dim=='0D':
        reactor_and_df = [compute_solutions_0D(*item,real_egr,species) for item in items]
        for i in range(len(items)):
            results += [reactor_and_df[i][1]]
    elif dim=='1D':
        results = [compute_solutions_1D(*item,restart_rate,real_egr,dry,T_reinj) for item in items]

    output=pd.concat(results[:],axis=0) 
    output.to_csv(save_file_name,index=False) 
    mpiprint(output)

def PRINT_MONO_CPU_WARNING():
    mpiprint("--------------------------------------------------",file=sys.stdout)
    mpiprint("WARNING, I AM BETTER FOR PARALLEL MPI CALCULATIONS",file=sys.stdout)
    mpiprint("--------------------------------------------------",file=sys.stdout)
