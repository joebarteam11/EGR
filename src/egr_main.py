#import cantera as ct
from egr import *
import time

# get the start time
st = time.time()

config = case('CH4:1.',                     #fuel compo
              'O2:1., N2:3.76',             #ox compo
              'CO2:1.0',           #egr compo
              3960.0,                       #thermal output power NOT IMPLEMENTED YET
              0.0001,                       #egr rate
              'vol'                         #egr rate unit
             )

#set reservoirs thermo-state
config.res.fuel_gas, config.res.fuel = create_reservoir(config.compo.fuel,'gri30.xml', 300.0, 100000)
config.res.ox_gas, config.res.ox = create_reservoir(config.compo.ox,'air.xml', 300.0, 100000)
config.res.egr_gas, config.res.egr = create_reservoir(config.compo.egr,'gri30.xml', 300.0, 100000)

#range of computation
egr_percentages = np.arange(0.0,0.21,0.1)
for egr_rate in egr_percentages:
    config.egr_rate = egr_rate #override config.egr_rate set during object instanciation
    phi_range = np.arange(0.6,2.1,0.05)
    reactor,results = compute_solutions(config,phi_range,power_regulator=True)
    #print(reactor.volume)
    subplot_data(phi_range,results,'Phi',['T[K]','Pth[W]','Y_O2','Y_CO2'],'EGR rate (%):'+str(round(config.egr_rate*100,1)))
    #print_reactor(reactor)

config.pow = 100000
for egr_rate in egr_percentages:
    config.egr_rate = egr_rate #override config.egr_rate set during object instanciation
    phi_range = np.arange(0.6,2.1,0.05)
    reactor,results = compute_solutions(config,phi_range,power_regulator=True)
    #print(reactor.volume)
    subplot_data(phi_range,results,'Phi',['T[K]','Pth[W]','Y_O2','Y_CO2'],'EGR rate (%):'+str(round(config.egr_rate*100,1)),symbol='--o')
    #print_reactor(reactor)

# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
    
see_graphs('Mode:'+config.mode
          +' | fuel:'+str(round(config.res.fuel.thermo.P,1)) +' Pa, '+ str(round(config.res.fuel.thermo.T,1))+' K'
          +' | ox:'+str(round(config.res.ox.thermo.P,1)) +' Pa, '+ str(round(config.res.ox.thermo.T,1))+' K'
          +' | egr:'+str(round(config.res.egr.thermo.P,1)) +' Pa, '+ str(round(config.res.egr.thermo.T,1))+' K'
          )