from lib_egr_260 import * #(compute_equilibrium,compute_solutions_0D,compute_solutions_1D)
import time
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

if __name__ == '__main__':
    path = os.getcwd()
    print('Current folder: ',path)
    print(f"Running Cantera version: {ct.__version__}")

    # get the start time
    st = time.time()
    temptlist = [1700]#[i for i in np.arange(300,501,50)]
    presslist = [10e5,18e5]#[i for i in np.arange(1e5,5.1e5,5e4)]
    config = case('H2:1.',                     #fuel compo
                  temptlist,                    #tin fuel
                  presslist,                        #pin fuel
                  'O2:1. N2:3.76',              #ox compo
                  temptlist,                    #tin ox
                  presslist,                        #pin ox
                  'CO2:1.',                     #egr compo
                  temptlist,                    #tin egr
                  presslist,                        #pin egr
                  [0],#[i for i in np.arange(0.60,2.51,0.1)],        #phi range
                  [i for i in np.linspace(0.0,1.0,50)],#[0.0,0.1,0.15,0.2],            #egr range
                  'mole',                       #egr rate unit
                  'schemes/Aramco13.cti',               #scheme
                  'Mix', #transport model
                  'NA'  #is an ARC chemistry ? 'ARC' = yes, other = no
                 )
    # temptlist = [i for i in np.arange(300,500,5)]
    # config = case('CH4:1.',                     #fuel compo
    #               temptlist,                    #tin fuel
    #               [1e5,2e5,3e5],                        #pin fuel
    #               'O2:1. N2:3.76',              #ox compo
    #               temptlist,                    #tin ox
    #               [1e5,2e5,3e5],                        #pin ox
    #               'CO2:1.',                     #egr compo
    #               temptlist,                    #tin egr
    #               [1e5,2e5,3e5],                        #pin egr
    #               [0.8,1,1.2],#[i for i in np.arange(0.80,1.22,0.05)],        #phi range
    #               [0.1,0.15,0.2],            #egr range
    #               'mole',                       #egr rate unit
    #               'gri30.cti',               #scheme
    #               'Mix', #transport model
    #               'NA'  #is an ARC chemistry ? 'ARC' = yes, other = no
    #              )

    Tins = [[t[i] for t in config.tin] for i in range(len(config.tin.fuel))]
    Pins = [[p[i] for p in config.pin] for i in range(len(config.pin.fuel))]

    items = [[config,phi,Tin,Pin] for phi in config.phi_range for Tin in Tins for Pin in Pins]
    #print(items)

    ncpu = mp.cpu_count()
    print('nCPU :',ncpu)
    print('nItems :',len(items))

    #Progress bar declaration
    pbar=tqdm(total=len(items),file=sys.stdout) 
    def update(*a):
        pbar.update()

    dim='0D'
    #species = ['H2','O2','N2','H2O']
    species = ['O2','CO','CO2']

    if(dim=='equilibrate'):
        #species = ['O2','CO','CO2']

        #Computation pool 
        pool = mp.Pool(min(len(items),ncpu))
        results = [pool.apply_async(compute_equilibrium, args=item+[species],callback=update) for item in items]
        pool.close()
        # wait for all tasks to complete and processes to close
        pool.join()
        
        #get results & store them in csv
        unpacked=[res.get() for res in results]
        output=pd.concat(unpacked,axis=0)
        output.to_csv(path+'/results'+'/plan_total_equilibrium'+'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv',index=False)

        print(output)

    elif(dim=='0D'):
        #reactor,pdresult = compute_solutions_0D(config,phi,Tin,Pin,real_egr=False,species = ['CH4','H2','O2','CO','CO2','H2O'])
        real_egr = False
        #species = ['O2','CO','CO2']
        res_time = 3
        dfs = []
        for item in items:
            reactor,pdresult = compute_solutions_0D(*item,real_egr,species,res_time)
            dfs.append(pdresult)
            try:
                update()
            except:
                pass
        
        output=pd.concat(dfs,axis=0)
        output.to_csv(path+'/results'+'/plan_total_equilibrium_KP_fco2'+
                      #'_'+time.strftime("%Y%m%d-%H%M%S")+
                      '.csv',index=False)
        print(output)

    elif(dim=='1D'):
        real_egr = False
        restart_rate = None # config.egr_range[0] #set to None if want to restart computation from the first egr value in egr_range
        

        #Computation pool 
        pool = mp.Pool(min(len(items),ncpu))
        results = [pool.apply_async(compute_solutions_1D, args=item+[restart_rate,real_egr,species],callback=update) for item in items]
        pool.close()
        # wait for all tasks to complete and processes to close
        pool.join()
        
        #get results & store them in csv
        unpacked=[res.get() for res in results]
        output=pd.concat(unpacked,axis=0)
        output.to_csv(path+'/CRASHTEST'+'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv',index=False)

        print(output)

    
    # get the execution time
    et = time.time()
    elapsed_time = et - st

    print('Execution time:', elapsed_time, 'seconds')