import cantera as ct
import numpy as np
from matplotlib.pylab import *
import pandas as pd
import time
from cantera import ck2cti
import csv
print(f"Running Cantera version: {ct.__version__}")
# plt.style.use("ggplot")
# plt.style.use("seaborn-pastel")

# plt.rcParams["axes.labelsize"] = 18
# plt.rcParams["xtick.labelsize"] = 14
# plt.rcParams["ytick.labelsize"] = 14
# plt.rcParams["figure.autolayout"] = True
# plt.rcParams["figure.dpi"] = 120
# phitab=[0.1,1,9]
# Temp =np.linspace(500,1000,10)
residence_time = 0.43  # s
reactor_volume = 0.4*10**(-6)  # m3
max_simulation_time = 50
reactor_pressure =  100 * ct.one_atm
fuel = 'CH3OH: 1'
oxidizer = 'O2:1.0, N2:3.76'
nt=200
gas = ct.Solution('alcool.cti') 

temp_dependence = ct.SolutionArray(gas)
with open('stirred_reactor.csv', 'w') as outfile:
    writer = csv.writer(outfile)
   

    # We will use concentrations from the previous iteration to speed up convergence
    gas.TP = 600, reactor_pressure
    gas.set_equivalence_ratio(1, fuel, oxidizer)
    GB=ct.Solution('alcool.cti')
    GB.TP = 600, reactor_pressure
    GB.set_equivalence_ratio(1, fuel, oxidizer)
    GB.equilibrate("HP")

    stirred_reactor = ct.IdealGasReactor(GB, energy="on", volume=reactor_volume)
    exhaust = ct.Reservoir(GB)
    fuel_air_mixture_tank = ct.Reservoir(gas)
    mass_flow_controller = ct.MassFlowController(upstream=fuel_air_mixture_tank,downstream=stirred_reactor, mdot=stirred_reactor.mass/residence_time)
    pressure_regulator = ct.PressureController(upstream=stirred_reactor, downstream=exhaust, master=mass_flow_controller)
    reactor_network = ct.ReactorNet([stirred_reactor])
    reactor_network.set_initial_time(0.0)
    # Re-run the isothermal simulations
    tic = time.time()
    #reactor_network.rtol = 1e-4
    #reactor_network.atol = 1e-4
    reactor_network.reinitialize()
    #time=0.0
    dt=1e-5
    temp=[]
    t=0
    for k in range(nt):
        #time=time+dt
       
        t=k*dt
        print(t)
        reactor_network.advance(t)
        
        toc = time.time()
        concentrations = stirred_reactor.thermo.X

        temp_dependence.append(stirred_reactor.thermo.state)
    #writer.writerow(temp_dependence.T,temp_dependence("CH3OH").X)

#plt.figure()
plt.plot(temp_dependence.T, temp_dependence("CH3OH").X, "r-", label="$CH_{3}OH$")
plt.xlabel("Temperature (K)")
plt.ylabel("Mole Fractions of CH3OH")
plt.show()



