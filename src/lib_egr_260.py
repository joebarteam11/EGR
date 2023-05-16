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

class FlameExtinguished(Exception):
    pass

class case:
    def __init__(self,fuel,Tinit_fuel,Pinit_fuel,ox,Tinit_ox,Pinit_ox,egr,Tinit_egr,Pinit_egr,phi_range,egr_range,egr_unit,scheme,isARC):
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
        self.isARC = True if isARC == 'ARC' else False

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
        egr_mass=(egr_rate/(1-egr_rate))*(air_mass) #+fuel_mass depending on the definition you want
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
        egr_mol=(egr_rate/(1-egr_rate))*(air_mol) #+fuel_mol depending on the definition you want
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
        egr_vol=(egr_rate/(1-egr_rate))*(air_vol) #+fuel_vol depending on the definition you want
        mdots=[fuel_vol*coef*config.gas.fuel.density_mass, air_vol*coef*config.gas.ox.density_mass, egr_vol*coef*config.gas.egr.density_mass]
        
        if(return_unit=='mass'):
            return mdots, sum(mdots)
        elif(return_unit=='mole'):
            warnings.warn('not tested yet, use mass instead', UserWarning)
            moldots=[m/n for m, n in zip(mdots,mol_weights)]
            return moldots, sum(moldots)

    else:
        raise ValueError('egr_unit must be mass, mole or vol')
    
def create_reservoir(config,content, T, P,scheme=None):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
    if(config.isARC):
        #remove extension from sheme variable
        ct.compile_fortran(config.scheme.split('.')[0]+'.f90')
    if(scheme is None):
        scheme = config.scheme
    gas = ct.Solution(scheme)
    gas.TPX = T, P, content
    return ct.Reservoir(gas), gas

def burned_gas(phi,config,egr_rate,ignition=True):
    if(config.isARC):
        #remove extension from sheme variable
        ct.compile_fortran(config.scheme.split('.')[0]+'.f90')
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
    
    sim.rtol = 1e-7
    sim.atol = 1e-8
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

def compute_equilibrium(config,phi,tin,pin,species = ['CH4','H2','O2','CO','CO2','H2O']):
    #create a dataframe naming colums with 'phi', 'T' and all the species in the list
    #then fill it with the values of phi, T and mole fractions of species using the concatenation of two dataframes, for each phi
    df =  pd.DataFrame(columns=['EGR','phi','T','P']+species)
    
    for egr in config.egr_range:
        for phi in config.phi_range:
            _,config.gas.fuel = create_reservoir(config,config.compo.fuel, tin[0], pin[0])
            _,config.gas.ox = create_reservoir(config,config.compo.ox, tin[1], pin[1],scheme='air.xml')
            _,config.gas.egr = create_reservoir(config,config.compo.egr, tin[2], pin[2])
            
            bg = burned_gas(phi,config,egr,ignition=True)
            
            df = pd.concat([df, pd.DataFrame([[egr, phi, bg.T, bg.P]+list(bg[species].X)], columns=['EGR','phi','T','P']+species)]).astype(float)
            try:
                update()
            except:
                pass

    return df


def compute_solutions_0D(config,phi,tin,pin,real_egr=False,species = ['CH4','H2','O2','CO','CO2','H2O']):
    #create a dataframe naming colums with 'phi', 'T' and all the species in the list
    #then fill it with the values of phi, T and mole fractions of species using the concatenation of two dataframes, for each phi
    df =  pd.DataFrame(columns=['EGR','phi','T','P']+species)
    if(version.parse(ct.__version__) >= version.parse("2.4.0")):
        print(('%10s %10s %10s %10s %10s %10s %10s' % ('phi','Xfuel', 'Xair', 'Xegr', 'HRR', 'T', 'P'))) 

    for egr in config.egr_range:
        for phi in config.phi_range:
            config.res.fuel,config.gas.fuel = create_reservoir(config,config.compo.fuel, tin[0], pin[0])
            config.res.ox,config.gas.ox = create_reservoir(config,config.compo.ox, tin[1], pin[1],scheme='air.xml')
            config.res.egr,config.gas.egr = create_reservoir(config,config.compo.egr, tin[2], pin[2])
            
            mfcs, reactor = reactor_0D(phi,config,egr,real_egr)
            #only with cantera >= 2.5
            if(version.parse(ct.__version__) >= version.parse("2.4.0")):
                moldot_tot = mfcs[0].mass_flow_rate/(config.gas.fuel.mean_molecular_weight/1000)+ mfcs[1].mass_flow_rate/(config.gas.ox.mean_molecular_weight/1000)+mfcs[2].mass_flow_rate/(config.gas.egr.mean_molecular_weight/1000)
                #print(('%10.3f %10.3f %10.3f %10.3e %10.3f %10.3f' % (mfcs[0].mass_flow_rate/mdot_tot, mfcs[1].mass_flow_rate/mdot_tot, (mfcs[2].mass_flow_rate/mdot_tot), reactor.thermo.heat_release_rate, reactor.T, reactor.thermo.P)))
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

def fresh_gas(phi,config,egr_rate,Tmix,Pmix):
    if(config.isARC):
        #remove extension from sheme variable
        ct.compile_fortran(config.scheme.split('.')[0]+'.f90')
    gas = ct.Solution(config.scheme)
    gas.TP = Tmix,Pmix

    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        gas.set_equivalence_ratio(phi, config.compo.fuel, config.compo.ox, 
                              basis="mole",
                              diluent=config.compo.egr,fraction={"diluent":egr_rate})
    else:
        gas.set_equivalence_ratio(phi, config.compo.fuel, config.compo.ox)

    return gas

def build_freeflame(mix,width=0.03):
    #mix = ct.Mixer(inlets=mfcs) #energy eq. is on by default
    flame = ct.FreeFlame(mix, width=width)
    return flame

def mixer(phi,config,egr,real_egr=False):
    if(config.isARC):
        #remove extension from sheme variable
        ct.compile_fortran(config.scheme.split('.')[0]+'.f90')
    gas=ct.Solution(config.scheme)
    mdots,mdot_tot = compute_mdots(config, egr, phi,return_unit='mass')

    fuel = ct.Quantity(gas,constant='HP')
    fuel.TPX = config.gas.fuel.T,config.gas.fuel.P,config.compo.fuel
    fuel.mass = mdots[0]/mdot_tot
    ox = ct.Quantity(gas,constant='HP')
    ox.TPX = config.gas.ox.T,config.gas.ox.P,config.compo.ox
    ox.mass = mdots[1]/mdot_tot
    
    if(real_egr):
        dilutent = ct.Quantity(config.gas.egr,constant='HP')
        dilutent.TP = config.gas.egr.T,config.gas.egr.P
    else:
        dilutent = ct.Quantity(gas,constant='HP')
        dilutent.TPX = config.gas.egr.T,config.gas.egr.P,config.compo.egr
    dilutent.mass = mdots[2]/mdot_tot

    mix = fuel+ox+dilutent
    #print(mix.report())
    return mix.T, mix.P, mix.X

def apply_egr_to_inlet(f,config,phi,egr):
    config.gas.egr.X = f.X[:,-1]
    _,_,X = mixer(phi,config,egr,real_egr=True)
    f.inlet.X = X
    return f

def flamme_thickness(f):
    #get the temperature profile
    T=f.T
    #get the grid
    X=f.grid
    #get the temperature at the flame base
    T0=T[0]
    #get the temperature at the flame tip
    T1=T[-1]
    #compute the maximal spatial temperature gradient along the flame
    dTdx=np.max(np.gradient(T,X))
    #compute the thickness
    thickness=(T1-T0)/dTdx
    #print('Thicknes (m)',thickness)
    return thickness

def compute_solutions_1D(config,phi,tin,pin,restart_rate,real_egr=False,vars=['EGR','phi','P','Tin','T','u','dF'],species = ['CH4','O2','CO2','H2O']):
    path = os.getcwd()+'/src'
    dfs=[]
    vartosave = vars+species
    df =  pd.DataFrame(columns=vartosave)

    for egr in config.egr_range:
        #create the gas object containing the mixture of fuel, ox and egr and all thermo data
        _, config.gas.fuel = create_reservoir(config,config.compo.fuel,tin[0], pin[0])
        _, config.gas.ox = create_reservoir(config,config.compo.ox, tin[1], pin[1],scheme='air.xml')
        _, config.gas.egr = create_reservoir(config,config.compo.egr, tin[2], pin[2])

        #create a dataframe naming colums with 'phi', 'T' and all the species in the list
        #then fill it with the values of phi, T and mole fractions of species using the concatenation of two dataframes, for each phi
        #vartosave = vars+species
        #df =  pd.DataFrame(columns=vartosave)

        #get the temperature and pressure of the mixture according to phi and egr rate
        T,P,X = mixer(phi, config, egr)
        f = build_freeflame(fresh_gas(phi,config,egr,T,P))
        #print(f.transport_model)
        tol_ss = [2.0e-7, 1.0e-7]  # tolerance [rtol atol] for steady-state problem
        tol_ts = [2.0e-7, 1.0e-7]  # tolerance [rtol atol] for time stepping

        f.flame.set_steady_tolerances(default=tol_ss)
        f.flame.set_transient_tolerances(default=tol_ts)
        f.transport_model = 'Mix'
        #print('flame X_CH4',f.inlet.thermo['CH4'].X)
        #print(('%10s %10s %10s %10s %10s %10s' % ('phi','Xfuel', 'Xair', 'Xegr', 'T', 'P'))) 
        #print(('%10.3f %10.3f %10.3f %10.3f %10.3f' % (phi, f['CH4'].X, f['O2'].X+f['N2'].X, f['CO2'].X, T, P)))
        flametitle=''
        if(restart_rate is None):
            flametitle = path+'/data/'+'egr'+str(round(egr,1))+'_phi'+str(round(phi,2))+'_T'+str(round(T,0))+'_P'+str(round(P/100000,0))+'_ARC_QC16.h5'
            f.set_initial_guess()
        else:
            flametitle = path+'/data/'+'egr'+str(round(restart_rate,1))+'_phi'+str(round(phi,2))+'_T'+str(round(T,0))+'_P'+str(round(P/100000,0))+'_ARC_QC16.h5'
            try:
                f.read_hdf(flametitle)
                #f.set_initial_guess(data=flametitle)
            except:
                raise Exception('Cannot restore flame from file '+flametitle)
            flametitle = path+'/data/'+'egr'+str(round(egr,1))+'_phi'+str(round(phi,2))+'_T'+str(round(T,0))+'_P'+str(round(P/100000,0))+'_ARC_QC16.h5'

        # print(f.X[:,-1])
        # config.gas.egr.X = f.X[:,-1]
        # _,_,X = mixer(phi,config,egr,real_egr=True)
        # print('Inlet composition',X)
        if(real_egr):
            f = apply_egr_to_inlet(f,config,phi,egr)
        else:
            f.inlet.T = T
            f.P = P
            f.inlet.X = X
        #print('Inlet composition',f.inlet.X)
        flame = solve_flame(f,flametitle,config,phi,egr,real_egr=real_egr)


        #restart_rate = egr
        if(version.parse(ct.__version__) >= version.parse('2.5.0')):
            SL0=f.velocity[0]
        else:
            SL0=f.u[0]

        index = [f.gas.species_index(specie) for specie in species]
        df = pd.concat([df, pd.DataFrame([[egr, phi, P, T, f.T[-1],SL0,flamme_thickness(f)]+list(f.X[index,-1])], columns=vartosave)]).astype(float) #+list(f.X[index][-1])] #
        #dfs = pd.concat([df],axis=0)
        try:
            update()
        except:
            print('Cannot update progress bar')
    return df

def solve_flame(f,flametitle,config,phi,egr,real_egr=False):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
    #################################################################
    # Iterations start here
    #################################################################
    verbose = 1
    loglevel  = 0                      # amount of diagnostic output (0 to 5)	    
    refine_grid = True                  # True to enable refinement, False to disable 	
    f.max_time_step_count=25000
    f.max_grid_points=500
    
    # first iteration
    ##################
    # No energy equilibrium activated (reduce the calculation)
    f.energy_enabled = False
    # Mesh refinement
    f.set_refine_criteria(ratio = 7.0, slope = 0.95, curve = 0.95)
    # Max number of times the Jacobian will be used before it must be re-evaluated
    f.set_max_jac_age(50, 50)
    #Set time steps whenever Newton convergence fails
    f.set_time_step(1e-06, [25, 40, 80, 140, 200, 350, 500, 700, 1000, 1300, 1700, 2000, 3000, 5000, 10000, 12000, 15000, 20000]) #s

    # Calculation
    if(verbose>0):
        print('1st iteration...')
    try:
        if(verbose>1):
            print('Inlet composition before f.solve',f.inlet.X)
        f.solve(loglevel, refine_grid)
        if(real_egr):
            f = apply_egr_to_inlet(f,config,phi,egr)
            if(verbose>0):
                print('EGR applied to inlet')
        #f.save(flametitle)
    except FlameExtinguished:
        print('Flame extinguished')
        
    except ct.CanteraError as e:
        print('Error occurred while solving:', e)
        
    # Second iteration
    #################
    # Energy equation activated
    f.energy_enabled = True
    # mesh refinement
    f.set_refine_criteria(ratio = 7.5, slope = 0.75, curve = 0.75)
    # Calculation
    if(verbose>0):
        print('2nd iteration...')
    try:
        if(verbose>1):
            print('Inlet composition before f.solve',f.inlet.X)
        f.solve(loglevel, refine_grid)
        if(real_egr):
            f = apply_egr_to_inlet(f,config,phi,egr)
            if(verbose>0):
                print('EGR applied to inlet')
        #f.save(flametitle)
    except FlameExtinguished:
        print('Flame extinguished')
        
    except ct.CanteraError as e:
        print('Error occurred while solving:', e)
        
    # Third iteration
    #################
    # On raffine le maillage
    f.set_refine_criteria(ratio = 5.0, slope = 0.4, curve = 0.4, prune = 0.03)
    # Calculation
    if(verbose>0):
        print('3rd iteration...')
    try:
        if(verbose>1):
            print('Inlet composition before f.solve',f.inlet.X)
        f.solve(loglevel, refine_grid)
        if(real_egr):
            f = apply_egr_to_inlet(f,config,phi,egr)
            if(verbose>0):
                print('EGR applied to inlet')
        #f.save(flametitle)
    except FlameExtinguished:
        print('Flame extinguished')
        
    except ct.CanteraError as e:
        print('Error occurred while solving:', e)
        
    #################
    # Fourth iteration
    # Mesh refinement
    f.set_refine_criteria(ratio = 5.0, slope = 0.1, curve = 0.1, prune = 0.01)
    # Calculation
    if(verbose>0):
        print('4th iteration...')
    try:
        if(verbose>1):
            print('Inlet composition before f.solve',f.inlet.X)
        f.solve(loglevel, refine_grid)
        if(real_egr):
            f = apply_egr_to_inlet(f,config,phi,egr)
            if(verbose>0):
                print('EGR applied to inlet')
        #f.save(flametitle)
    except FlameExtinguished:
        print('Flame extinguished')
        
    except ct.CanteraError as e:
        print('Error occurred while solving:', e)
        
    # Fifth iteration
    #################
    #Mesh refinement
    f.set_refine_criteria(ratio = 5.0, slope = 0.05, curve = 0.05, prune = 0.01)
    # Calculation
    if(verbose>0):
        print('5th iteration...')
    try:
        if(verbose>1):
            print('Inlet composition before f.solve',f.inlet.X)
        f.solve(loglevel, refine_grid)
        if(real_egr):
            f = apply_egr_to_inlet(f,config,phi,egr)
            if(verbose>0):
                print('EGR applied to inlet')
        try:
            f.write_hdf(flametitle)
        except:
            f.save(flametitle[:-3]+'.yaml')
    except FlameExtinguished:
        print('Flame extinguished')
        
    except ct.CanteraError as e:
        print('Error occurred while solving:', e)
    # # Sixth iteration
    #################
    #NO Mesh refinement
    #f.set_refine_criteria(ratio = 5.0, slope = 0.05, curve = 0.05, prune = 0.01)
    # Calculation
    if(verbose>0):
        print('6th iteration...')
    try:
        if(verbose>1):
            print('Inlet composition before f.solve',f.inlet.X)
        f.solve(loglevel, refine_grid=False)
        if(real_egr):
            f = apply_egr_to_inlet(f,config,phi,egr)
            if(verbose>0):
                print('EGR applied to inlet')
        #f.save(flametitle)
    except FlameExtinguished:
        print('Flame extinguished')
        
    except ct.CanteraError as e:
        print('Error occurred while solving:', e)
    
    # # Seventh iteration
    #################
    #NO Mesh refinement
    #f.set_refine_criteria(ratio = 5.0, slope = 0.05, curve = 0.05, prune = 0.01)
    # Calculation
    if(verbose>0):
        print('7th iteration...')
    try:
        if(verbose>1):
            print('Inlet composition before f.solve',f.inlet.X)
        f.solve(loglevel, refine_grid=False)
        if(real_egr):
            f = apply_egr_to_inlet(f,config,phi,egr)
            if(verbose>0):
                print('EGR applied to inlet')
        #f.save(flametitle)
        try:
            f.write_hdf(flametitle)
        except:
            f.save(flametitle[:-3]+'.yaml')
    except FlameExtinguished:
        print('Flame extinguished')
        
    except ct.CanteraError as e:
        print('Error occurred while solving:', e)
    
    return f

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
                  [300,500],                    #tin fuel
                  [1e5,5e5],                        #pin fuel
                  'O2:1. N2:3.76',              #ox compo
                  [300,500],                    #tin ox
                  [1e5,5e5],                        #pin ox
                  'CO2:1.',                     #egr compo
                  [300,500],                    #tin egr
                  [1e5,5e5],                        #pin egr
                  [i for i in np.arange(0.80,1.22,0.05)],        #phi range
                  [0.0,0.1,0.3,0.5],            #egr range
                  'mole',                       #egr rate unit
                  'schemes/CH4_16_250_10_QC.cti',               #scheme
                  'ARC', 
                 )
    

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

    dim='1D'

    if(dim=='0D'):

        dfs=[]
        for p in Pins:
            #set reservoirs thermo-state
            config.res.fuel,config.gas.fuel = create_reservoir(config,config.compo.fuel, 300.0, p[0])
            config.res.ox,config.gas.ox = create_reservoir(config,config.compo.ox, 300.0, p[1],scheme='air.xml')
            config.res.egr,config.gas.egr = create_reservoir(config,config.compo.egr, 300.0, p[2])

            reactor,pdresult = compute_solutions_0D(config,real_egr=False,species = ['CH4','H2','O2','CO','CO2','H2O'])
            
            dfs.append(pdresult)

        dfs=pd.concat(dfs,axis=0)
        dfs.to_csv(path+'/results'+'/plan_partiel_0D_dilution_'+time.strftime("%Y%m%d-%H%M%S")+'.csv',index=False)
        print(dfs)

    elif(dim=='1D'):
        real_egr = False
        restart_rate = None # config.egr_range[0] #set to None if want to compute from the first egr value in egr_range
        #for phi in config.phi_range:
            #Computation pool 
        pool = mp.Pool(min(len(items),ncpu))
        results = [pool.apply_async(compute_solutions_1D, args=item+[restart_rate,real_egr],callback=update) for item in items]
        pool.close()
        # wait for all tasks to complete and processes to close
        pool.join()
        
        #get results & store them in csv
        unpacked=[res.get() for res in results]
        output=pd.concat(unpacked,axis=0)
        output.to_csv(path+'/plan_total_QC_16_ARC_canavbp'+'_'+time.strftime("%Y%m%d-%H%M%S")+'.csv',index=False)

        print(output)

    
    # get the execution time
    et = time.time()
    elapsed_time = et - st

    print('Execution time:', elapsed_time, 'seconds')