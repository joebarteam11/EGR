import sys,os
import time
import numpy as np
import cantera as ct
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing as mp
from packaging import version
#from alive_progress import alive_bar
from tqdm import tqdm
import warnings

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
warnings.simplefilter("ignore", UserWarning)

class case:
    def __init__(self,fuel,Tinit_fuel,Pinit_fuel,ox,Tinit_ox,Pinit_ox,egr,Tinit_egr,Pinit_egr,phi_range,egr_range,egr_unit,scheme):
        self.compo = self.Compo(fuel,ox,egr)
        self.res = self.Reservoirs()
        self.gas = self.Gas(fuel,ox,egr)
        
        #try:
        if(len(Tinit_fuel) != len(Tinit_ox) or len(Tinit_fuel) != len(Tinit_egr)):
            raise ValueError('length of Tinit_fuel, Tinit_ox, Tinit_egr must be equal')
        else:
            self.tin = self.Tinit(Tinit_fuel,Tinit_ox,Tinit_egr)
        #except:
        #    raise ValueError('length of Tinit_fuel, Tinit_ox, Tinit_egr must be equal')
        #try:
        if(len(Pinit_fuel) != len(Pinit_ox) or len(Pinit_fuel) != len(Pinit_egr)):
            raise ValueError('length of Pinit_fuel, Pinit_ox, Pinit_egr must be equal')
        else:
            self.pin = self.Pinit(Pinit_fuel,Pinit_ox,Pinit_egr)
        #except:
        #    raise ValueError('length of Pinit_fuel, Pinit_ox, Pinit_egr must be equal')
        
        self.phi_range = phi_range
        self.egr_range = egr_range
        self.egr_unit = egr_unit
        self.scheme = scheme

    class Gas:
        def __init__(self,fuel,ox,egr):
            self.fuel = fuel
            self.ox = ox
            self.egr = egr

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
    
    class Tinit:
        def __init__(self,Tfuel,Tox,Tegr):
            self.fuel = Tfuel
            self.ox = Tox
            self.egr = Tegr

        def __iter__(self):
            yield self.fuel
            yield self.ox
            yield self.egr
    
    class Pinit:
        def __init__(self,Pfuel,Pox,Pegr):
            self.fuel = Pfuel
            self.ox = Pox
            self.egr = Pegr

        def __iter__(self):
            yield self.fuel
            yield self.ox
            yield self.egr

def compute_mdots(config,egr_rate,phi,coef=50.0,return_unit='mass'): 
    mw_fuel = config.gas.fuel.mean_molecular_weight/1000 #kg/mol
    mw_oxidizer = config.gas.ox.mean_molecular_weight/1000 #kg/mol
    mw_egr = config.gas.egr.mean_molecular_weight/1000 #kg/mol
    mol_weights=[mw_fuel,mw_oxidizer,mw_egr]

    warnings.warn('compute_mdots() is setted for methane-air mixtures only, change Sx, Sy and n_mole_ox for others reactants', UserWarning)
    Sx=2.0
    Sy=17.16
    n_mole_ox=4.76
    
    if(config.egr_unit=='mass'):
        fuel_mass=phi/(Sy+phi)
        air_mass=1-fuel_mass
        egr_mass=(egr_rate/(1-egr_rate))*(fuel_mass+air_mass)
        mdots=[fuel_mass*coef, air_mass*coef, egr_mass*coef]
        
        if(return_unit=='mass'):
            return mdots, sum(mdots)
        elif(return_unit=='mole'):
            warnings.warn('not tested yet, use mass instead', UserWarning)
            moldots=[m/n for m, n in zip(mdots,mol_weights)]
            return moldots, sum(moldots)
    
    elif(config.egr_unit=='mole'):
        fuel_mol = phi/(n_mole_ox*Sx+phi)
        air_mol = 1-fuel_mol
        egr_mol=(egr_rate/(1-egr_rate))*(fuel_mol)
        mdots=[fuel_mol*coef*mw_fuel, air_mol*coef*mw_oxidizer, egr_mol*coef*mw_egr]
        
        if(return_unit=='mass'):
            return mdots, sum(mdots)
        elif(return_unit=='mole'):
            warnings.warn('not tested yet, use mass instead', UserWarning)
            moldots=[m/n for m, n in zip(mdots,mol_weights)]
            return moldots, sum(moldots)

    elif(config.egr_unit=='vol'):
        coef=coef/1000
        fuel_mol = phi/(n_mole_ox*Sx+phi)
        fuel_vol = fuel_mol * config.gas.fuel.volume_mole
        air_vol = (1-fuel_mol) * config.gas.ox.volume_mole
        egr_vol=(egr_rate/(1-egr_rate))*(fuel_vol+air_vol)
        mdots=[fuel_vol*coef*config.gas.fuel.density_mass, air_vol*coef*config.gas.ox.density_mass, egr_vol*coef*config.gas.egr.density_mass]
        
        if(return_unit=='mass'):
            return mdots, sum(mdots)
        elif(return_unit=='mole'):
            warnings.warn('not tested yet, use mass instead', UserWarning)
            moldots=[m/n for m, n in zip(mdots,mol_weights)]
            return moldots, sum(moldots)

    else:
        raise ValueError('egr_unit must be mass, mole or vol')
    
def create_reservoir(content, scheme, T, P):
    gas = ct.Solution(scheme)
    gas.TPX = T, P, content
    return ct.Reservoir(gas), gas

def burned_gas(phi,config,egr_rate,ignition=True):
    gas = ct.Solution(config.scheme)
    gas.TP = config.gas.ox.T,config.gas.ox.P
    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        gas.set_equivalence_ratio(phi, config.compo.fuel, config.compo.ox,
            #additional params only compatible with cantera >= 2.5
            basis="mole", diluent=config.compo.egr,fraction={"diluent":egr_rate})
    else:
        gas.set_equivalence_ratio(phi, config.compo.fuel, config.compo.ox)
    
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

def reactor_0D(phi,config,egr_rate,real_egr,steady_state_only=True):
    #create the reactor and fill it with burned gas to ignite the mixture
    gb = burned_gas(phi,config,egr_rate) 
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
    sim.atol = 1e-11
    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        sim.max_steps = 100000 

    #compute steady state
    if(steady_state_only):
        sim.advance_to_steady_state()#only with cantera >= 2.5
    # else:
    #     residence_time = 1  # starting residence time
    #     while reactor.T > 500:
    #         sim.set_initial_time(0.0)  # reset the integrator
    #         sim.advance_to_steady_state()
    #         print('tres = {:.2e}; T = {:.1f}'.format(residence_time, reactor.T))
    #         #states.append(reactor.thermo.state, tres=residence_time)
    #         residence_time *= 0.9  # decrease the residence time for the next iteration
    #sim.set_initial_time(0.0)
    #sim.advance(10000) #large number to ensure that the steady state is reached i.e. gas as travelled enough trough the reactor volume

    return mfcs, reactor

def compute_solutions_0D(config,real_egr=False,species = ['CH4','H2','O2','CO','CO2','H2O']):
    #create a dataframe naming colums with 'phi', 'T' and all the species in the list
    #then fill it with the values of phi, T and mole fractions of species using the concatenation of two dataframes, for each phi
    df =  pd.DataFrame(columns=['EGR','phi','T','P']+species)
    if(version.parse(ct.__version__) >= version.parse("2.4.0")):
        print(('%10s %10s %10s %10s %10s %10s %10s' % ('phi','Xfuel', 'Xair', 'Xegr', 'HRR', 'T', 'P'))) 

    for egr in config.egr_range:
        for phi in config.phi_range:
            mfcs, reactor = reactor_0D(phi,config,egr,real_egr)
            #only with cantera >= 2.5
            moldot_tot = mfcs[0].mass_flow_rate/(config.gas.fuel.mean_molecular_weight/1000)+ mfcs[1].mass_flow_rate/(config.gas.ox.mean_molecular_weight/1000)+mfcs[2].mass_flow_rate/(config.gas.egr.mean_molecular_weight/1000)
            #print(('%10.3f %10.3f %10.3f %10.3e %10.3f %10.3f' % (mfcs[0].mass_flow_rate/mdot_tot, mfcs[1].mass_flow_rate/mdot_tot, (mfcs[2].mass_flow_rate/mdot_tot), reactor.thermo.heat_release_rate, reactor.T, reactor.thermo.P)))
            if(version.parse(ct.__version__) >= version.parse("2.4.0")):
                print((' %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f' % (phi, #reactor.thermo['CH4'].X, reactor.thermo['O2'].X+reactor.thermo['N2'].X, reactor.thermo['CO2'].X,
                                                                    #only with cantera >= 2.5
                                                                    (mfcs[0].mass_flow_rate/(config.gas.fuel.mean_molecular_weight/1000))/moldot_tot,
                                                                    (mfcs[1].mass_flow_rate/(config.gas.ox.mean_molecular_weight/1000))/moldot_tot, 
                                                                    (mfcs[2].mass_flow_rate/(config.gas.egr.mean_molecular_weight/1000))/moldot_tot, 
                                                                    reactor.thermo.heat_release_rate, #only with cantera >= 2.5
                                                                    #'N/A',
                                                                    reactor.T, reactor.thermo.P)))

            df = pd.concat([df, pd.DataFrame([[egr, phi, reactor.T, reactor.thermo.P]+list(reactor.thermo[species].X)], columns=['EGR','phi','T','P']+species)]).astype(float)
            try:
                update()
            except:
                pass

    return reactor, df

def show_graphs(df,title=None,labels=None,xlabel=None,ylabel=None,style='-o',subplot=1,ax=None,plot=True,save=False,path=None):
    warnings.warn("This function is still under developpement and may not work properly", UserWarning)
    if(ax==None):
        fig, ax = plt.subplots(1,subplot)
    if(version.parse(matplotlib.__version__) >= version.parse('3.5.0')):
        df.plot(ax=ax,style=style,title=title,xlabel=xlabel,ylabel=ylabel)
    else:
        df.plot(ax=ax,style=style,title=title)
        
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if(plot):
        plt.grid()
        plt.legend(loc='best')
        ax.legend(labels)
        plt.tight_layout()
    if(save):
        print('Saving plot')
        plt.savefig(path+ylabel+'_f('+xlabel+').png', dpi=300, bbox_inches='tight')
    #else:
    plt.show(block=True)

    return ax

if __name__ == '__main__':
    path = os.getcwd()
    print('Current folder: ',path)
    print(f"Running Cantera version: {ct.__version__}")

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
                  [0.0,0.1,0.3,0.5],            #egr range
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

    dim='0D'
    if(dim=='0D'):
        dfs=[]
        #pressures = [100_000,500_000] #Pa
        for p in Pins:
            #set reservoirs thermo-state
            config.res.fuel,config.gas.fuel = create_reservoir(config.compo.fuel,config.scheme, 300.0, p[0])
            config.res.ox,config.gas.ox = create_reservoir(config.compo.ox,'air.xml', 300.0, p[1])
            config.res.egr,config.gas.egr = create_reservoir(config.compo.egr,config.scheme, 300.0, p[2])

            reactor,pdresult = compute_solutions_0D(config,real_egr=False,species = ['CH4','H2','O2','CO','CO2','H2O'])
            
            dfs.append(pdresult)
        dfs=pd.concat(dfs,axis=0)
        print(dfs)
        dfs.to_csv(path+'/results'+'/plan_partiel_0D_dilution_'+time.strftime("%Y%m%d-%H%M%S")+'.csv',index=False)
    else:
        pass
        # previous_rate = 0.3
        # for egr in config.egr_range:
        #     #Computation pool 
        #     pool = mp.Pool(min(len(items),ncpu))
        #     results = [pool.apply_async(compute_solutions_1D, args=item+[egr,previous_rate], callback=update) for item in items]
        #     pool.close()
        #     # wait for all tasks to complete and processes to close
        #     pool.join()
        #     previous_rate = egr
        #     #get results & store them in csv
        #     unpacked=[res.get() for res in results]
        #     output=pd.concat(unpacked,axis=0)
        #     output.to_csv(path+'/plan_partiel_0D_dilution_'+str(round(egr,1))+'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv',index=False)

        #     print(output)
    
    # get the execution time
    et = time.time()
    elapsed_time = et - st

    print('Execution time:', elapsed_time, 'seconds')