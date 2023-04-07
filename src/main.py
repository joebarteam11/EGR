from lib_egr_260 import *
import time
import pandas as pd

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

if __name__ == '__main__':
    path = os.getcwd()
    print('Current folder: ',path)
    print(f"Running Cantera version: {ct.__version__}")
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
                  [i for i in np.arange(0.8,1.21,0.05)],        #phi range
                  [0.5,0.7],            #egr range
                  'mole',                       #egr rate unit
                  'schemes/Aramco13.cti'                   #scheme
                 )
    
    Tins = [[t[i] for t in config.tin] for i in range(len(config.tin.fuel))]
    Pins = [[p[i] for p in config.pin] for i in range(len(config.pin.fuel))]

    items = [[config,phi,Tin,Pin] for phi in config.phi_range for Tin in Tins for Pin in Pins]
    #print(items)

    ncpu = mp.cpu_count()
    print('nCPU :',ncpu)
    print('nItems :',len(items))

    #Progress bar declaration
    pbar=tqdm(total=len(items)*len(config.egr_range),file=sys.stdout) 
    def update(*a):
        pbar.update()

    if(True):
        previous_rate = 0.3
        for egr in config.egr_range:
            #Computation pool 
            pool = mp.Pool(min(len(items),ncpu))
            results = [pool.apply_async(compute_solutions_1D, args=item+[egr,previous_rate], callback=update) for item in items]
            pool.close()
            # wait for all tasks to complete and processes to close
            pool.join()
            previous_rate = egr
            #get results & store them in csv
            unpacked=[res.get() for res in results]
            output=pd.concat(unpacked,axis=0)
            output.to_csv(path+'/plan_partiel_dilution_'+str(round(egr,1))+'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv',index=False)

            print(output)

    
    # get the execution time
    et = time.time()
    elapsed_time = et - st

    print('Execution time:', elapsed_time, 'seconds')