import sys
import time
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from matplotlib import *
import pandas as pd

class case:
    def __init__(self,fuel,ox,egr,phi_range,egr_range,egr_unit,scheme):
        self.compo = self.Compo(fuel,ox,egr)
        self.res = self.Reservoirs()
        self.phi_range = phi_range
        self.egr_range = egr_range
        self.egr_unit = egr_unit
        self.scheme = scheme

    class Compo:
        def __init__(self,fuel,ox,egr):
            self.fuel = fuel
            self.ox = ox
            self.egr = egr
    
    class Reservoirs:
        def __init__(self):
            self.fuel = None
            self.ox = None
            self.egr = None

def compute_mdots(config,egr_rate,phi,coef=50.0): 
    mw_fuel = config.res.fuel.thermo.mean_molecular_weight/1000 #kg/mol
    mw_oxidizer = config.res.ox.thermo.mean_molecular_weight/1000 #kg/mol
    mw_egr = config.res.egr.thermo.mean_molecular_weight/1000 #kg/mol

    Sx=2.0
    Sy=17.16
    n_mole_ox=4.76
    
    if(config.egr_unit=='mass'):
        fuel_mass=phi/(Sy+phi)
        air_mass=1-fuel_mass
        egr_mass=(egr_rate/(1-egr_rate))*(fuel_mass+air_mass)
        mdots=[fuel_mass*coef, air_mass*coef, egr_mass*coef]
        #print(sum(mdots))
        return mdots, sum(mdots)
    
    elif(config.egr_unit=='mole'):
        fuel_mol = phi/(n_mole_ox*Sx+phi)
        air_mol = 1-fuel_mol
        egr_mol=(egr_rate/(1-egr_rate))*(fuel_mol+air_mol)
        mdots=[fuel_mol*coef*mw_fuel, air_mol*coef*mw_oxidizer, egr_mol*coef*mw_egr]
        #print(sum(mdots))
        return mdots, sum(mdots)

    elif(config.egr_unit=='vol'):
        coef=coef/1000
        fuel_mol = phi/(n_mole_ox*Sx+phi)
        fuel_vol = fuel_mol * config.res.fuel.thermo.volume_mole
        air_vol = (1-fuel_mol) * config.res.ox.thermo.volume_mole
        egr_vol=(egr_rate/(1-egr_rate))*(fuel_vol+air_vol)
        mdots=[fuel_vol*coef*config.res.fuel.thermo.density_mass, air_vol*coef*config.res.ox.thermo.density_mass, egr_vol*coef*config.res.egr.thermo.density_mass]
        #print(sum(mdots))
        return mdots, sum(mdots)

    else:
        raise ValueError('egr_unit must be mass, mole or vol')
    
def create_reservoir(content, scheme, T, P):
    gas = ct.Solution(scheme)
    gas.TP = T, P
    gas.set_equivalence_ratio(1.0, content, 'O2:1.0, N2:3.76')
    return ct.Reservoir(gas)

def burned_gas(phi,config,egr_rate,scheme,ignition=True):
    gas = ct.Solution(scheme)
    gas.TP = config.res.ox.thermo.T,config.res.ox.thermo.P
    gas.set_equivalence_ratio(phi, config.compo.fuel, config.compo.ox)
        #additional params only compatible with cantera >= 2.5
        #basis="mole", diluent=config.compo.egr,fraction={"diluent":egr_rate})
    if(ignition):
        gas.equilibrate("HP")
    return gas

def build_reactor(mixture,volume):
    mix = ct.IdealGasReactor(mixture, energy='on',volume=volume) #energy eq. is on by default
    exhaust = ct.Reservoir(mixture)
    pressure_regulator = ct.Valve(mix, exhaust, K=10.0)
    return mix,pressure_regulator

def init_reservoirs_mdots(amont, mdots, aval):
    mfc=[]
    i=0
    for res in amont:
        mfc.append(ct.MassFlowController(res, aval, mdot=mdots[i]))
        i+=1
    return mfc

def edit_reservoirs_mdots(reactor,mdots):
    i=0
    for inlet in reactor.inlets:
        inlet.set_mass_flow_rate = mdots[i]
        i+=1

def reactor_0D(phi,config,egr_rate,real_egr):
    #create the reactor and fill it with burned gas to ignite the mixture
    gb = burned_gas(phi,config,egr_rate,scheme=config.scheme) 
    reactor,pressure_reg = build_reactor(gb,volume=1000.0)#high volume to ensure that the residence time is high enough
    reservoirs=[config.res.fuel, config.res.ox, config.res.egr]
    #create the reactor network
    sim=ct.ReactorNet([reactor])
    sim.reinitialize()
    #set the mass flow controllers according to phi and egr rate asked in the config
    mdots,mdot_tot = compute_mdots(config, egr_rate, phi)
    
    mfcs = None
    #create the mass flow controllers according to the number of reservoirs
    if(real_egr):
        mfcs = init_reservoirs_mdots([config.res.fuel, config.res.ox,reactor],mdots,reactor)## ajouter la detection du nombre de reservoir dans la config

    else:
        mfcs = init_reservoirs_mdots(reservoirs,mdots,reactor)## ajouter la detection du nombre de reservoir dans la config
    
    sim.rtol = 1e-11
    sim.atol = 1e-12
    #sim.max_steps = 1000000 #not in cantera < 2.5

    #compute steady state
    #sim.advance_to_steady_state()#only with cantera >= 2.5
    sim.set_initial_time(0.0)
    sim.advance(10000) #large number to ensure that the steady state is reached i.e. gas as travelled enough trough the reactor volume

    return mfcs, reactor

def compute_solutions_0D(config,real_egr=False,species = ['CH4','H2','O2','CO2','H2O']):
    #create a dataframe naming colums with 'phi', 'T' and all the species in the list
    #then fill it with the values of phi, T and mole fractions of species using the concatenation of two dataframes, for each phi
    df =  pd.DataFrame(columns=['EGR','phi','T','P']+species)
    print(('%10s %10s %10s %10s %10s %10s %10s' % ('phi','Xfuel', 'Xair', 'Xegr', 'HRR', 'T', 'P'))) 
    for egr in config.egr_range:
        for phi in config.phi_range:
            mfcs, reactor = reactor_0D(phi,config,egr,real_egr)
            #only with cantera >= 2.5
            #moldot_tot = mfcs[0].mass_flow_rate/(config.res.fuel.thermo.mean_molecular_weight/1000)+ mfcs[1].mass_flow_rate/(config.res.ox.thermo.mean_molecular_weight/1000)+mfcs[2].mass_flow_rate/(config.res.egr.thermo.mean_molecular_weight/1000)
            #print(('%10.3f %10.3f %10.3f %10.3e %10.3f %10.3f' % (mfcs[0].mass_flow_rate/mdot_tot, mfcs[1].mass_flow_rate/mdot_tot, (mfcs[2].mass_flow_rate/mdot_tot), reactor.thermo.heat_release_rate, reactor.T, reactor.thermo.P)))
            print((' %10.3f %10.3f %10.3f %10.3f %10s %10.3f %10.3f' % (phi, reactor.thermo['CH4'].X, reactor.thermo['O2'].X+reactor.thermo['N2'].X, reactor.thermo['CO2'].X,
                                                                    #only with cantera >= 2.5
                                                                    #(mfcs[0].mass_flow_rate/(config.res.fuel.thermo.mean_molecular_weight/1000))/moldot_tot,
                                                                    #(mfcs[1].mass_flow_rate/(config.res.ox.thermo.mean_molecular_weight/1000))/moldot_tot, 
                                                                    #(mfcs[2].mass_flow_rate/(config.res.egr.thermo.mean_molecular_weight/1000))/moldot_tot, 
                                                                    #reactor.thermo.heat_release_rate, #only with cantera >= 2.5
                                                                    'N/A',
                                                                    reactor.T, reactor.thermo.P)))

            df = pd.concat([df, pd.DataFrame([[egr, phi, reactor.T, reactor.thermo.P]+list(reactor.thermo[species].X)], columns=['EGR','phi','T','P']+species)]).astype(float)
        
    return reactor, df

def show_graphs(df,title,labels,xlabel,ylabel,style='-o'):
    fig, ax = plt.subplots(1,1)
    df.plot(ax=ax, style=style,title=title,xlabel=xlabel,ylabel=ylabel)
    plt.grid()
    plt.legend(loc='best')
    ax.legend(labels)
    plt.show()

if __name__ == '__main__':
    # get the start time
    st = time.time()

    config = case('CH30H:1.',                     #fuel compo
                  'O2:1. N2:3.76',              #ox compo
                  'CO2:1.',                     #egr compo
                  [0.8,0.9,1.0,1.1,1.2],        #phi range
                  [0.0],                    #egr range
                  'mole',                       #egr rate unit
                  'CRECK_new.cti'                   #scheme
                 )
    dfs=[]
    pressures = [10_000_000] #Pa
    for p in pressures:
        #set reservoirs thermo-state
        config.res.fuel = create_reservoir(config.compo.fuel,config.scheme, 600.0, p)
        config.res.ox = create_reservoir(config.compo.ox,'air.xml', 600.0, p)
        config.res.egr = create_reservoir(config.compo.egr,config.scheme, 600.0, p)

        reactor,pdresult = compute_solutions_0D(config,real_egr=False)
        
        dfs.append(pdresult)
    
    print(dfs)

    # get the end time
    et = time.time()
    # get the execution time
    elapsed_time = et - st

        
    #select the columns to plot (X,Xbis,Y)    
    dfsc=[df.pivot_table(index='phi',columns='EGR',values='T') for df in dfs]
    dfsc=pd.concat(dfsc, axis = 1, keys = pressures)
    
    print(dfsc)
    print('Execution time:', elapsed_time, 'seconds')

    title='(0D) Flame temperature vs equivalence ratio (Tin_EGR:'+str(config.res.egr.thermo.T)+'K)'
    human_labels = [str(round(p/100000,1))+' bar, '+str(round(e*100,1))+"%EGR"+config.egr_unit for p in pressures for e in config.egr_range]
    xlabel='Equivalence ratio'
    ylabel='Tad [K]'

    show_graphs(dfsc,title,human_labels,xlabel,ylabel)