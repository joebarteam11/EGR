#import cantera as ct
from egr import *

config = case('CH4:1.',                     #fuel compo
              'O2:1., N2:3.76',             #ox compo
              'CO2:1.0, H2O:1.0',           #egr compo
              3960.0,                       #thermal output power NOT IMPLEMENTED YET
              0.0001,                       #egr rate
              'vol'                         #egr rate unit
             )

#reservoirs def
config.res.fuel_gas, config.res.fuel = create_reservoir(config.compo.fuel,'gri30.xml', 300.0, 100000)
config.res.ox_gas, config.res.ox = create_reservoir(config.compo.ox,'air.xml', 300.0, 100000)
config.res.egr_gas, config.res.egr = create_reservoir(config.compo.egr,'gri30.xml', 800.0, 100000)

#range of computation
egr_percentages = np.arange(0.00,0.25,0.05)
for egr_rate in egr_percentages:
    config.egr_rate = egr_rate #override config.egr_rate set during object instanciation
    phi_range = np.arange(0.45,1.5,0.05)
    reactor,results = compute_solutions(config,phi_range,print_report=False,real_egr=True)
    #print_reactor(reactor)
    subplot_data(phi_range,results,'Phi',['T[K]','HRR[W/m3]','Y_O2','Y_CO2'],'EGR'+config.egr_unit+'%:'+str(round(egr_rate*100,1))+'%')

see_graphs('Mode:'+config.mode
          +' | fuel:'+str(round(config.res.fuel.thermo.P,1)) +' Pa, '+ str(round(config.res.fuel.thermo.T,1))+' K'
          +' | ox:'+str(round(config.res.ox.thermo.P,1)) +' Pa, '+ str(round(config.res.ox.thermo.T,1))+' K'
          +' | egr:'+str(round(config.res.egr.thermo.P,1)) +' Pa, '+ str(round(config.res.egr.thermo.T,1))+' K'
          )