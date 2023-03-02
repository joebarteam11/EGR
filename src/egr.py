import sys
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from matplotlib import *
import pandas as pd

def get_ox_mdot_v1(egr_percentage, ox_compo, egr_compo, mode='mdot_fixed', unit='vol', phi=1.0):
    try:
        if(phi<0.0 or phi>10.0):
            raise Exception('phi must be in range ]0.0 ; 10;0]')
        mw_air = ox_compo.mean_molecular_weight/1000 #kg/mol
        rho_air = ox_compo.density_mass #kg/m3
        mw_egr = egr_compo.mean_molecular_weight/1000 #kg/mol
        rho_egr = egr_compo.density_mass #kg/m3

        x_o2 = ox_compo.X[ox_compo.species_index('O2')]
        mdot_air = ((mw_air/x_o2)*2.0)/phi

        try:
            if(unit=='vol'):
                mdot_air_vol = mdot_air/rho_air
                mdot_egr_vol = mdot_air_vol * egr_percentage #volumetric
                mdot_egr = mdot_egr_vol * rho_egr
            elif(unit=='mol'):
                mdot_air_mol = mdot_air/mw_air
                mdot_egr_mol = mdot_air_mol * egr_percentage #molar
                mdot_egr = mdot_egr_mol * mw_egr
            elif(unit=='mass'): 
                mdot_egr = mdot_air * egr_percentage #mass
            else :
                raise Exception('Unknown unit')
            #print('Imposed egr massflow rate: '+ str(mdot_egr)+ '[kg/s]')
            try:
                if(mode=='mdot_fixed'):
                    mdot_air -= mdot_egr
                elif(mode=='mdot_free'):
                    pass
                else:
                    raise Exception('Unknown mode')
                #print('Imposed Air massflow rate: '+ str(mdot_air)+ '[kg/s]')
                return mdot_air, mdot_egr
            except:
                print("Unknown mode '"+ mode +"', please choose between ['mdot_fixed','mdot_free'] ")
        except:
            print("Unknown unit '"+ unit +"', please choose between ['mol','mass','vol']")

    except:
        print('The oxidizer provided does not contain pure O2 !')
    
get_ox_mdot_v1.__doc__="Calcule un débit massique a imposer en focntion du pourcentage de egr qu'on souhaite ajouter Ce pourcentage peut etre défini comme un pourcentage de 'remplacement' d'air : on va réduired'autant le débit d'air pour conserver un débit massique total constant défini par le débit d'air sans egr à la richesse souhaitée (mode mdot_ox_fixed) Il peut aussi être considéré comme un pourcentage de egr supplémentaire, défini comme un ajout de X% du débit d'air à la richesse considérée (mode mdot_ox_free)"

def compute_mdots(phi,config,reactor,coeff=1.0): #, egr_rate, mdot_out, mdot_o2_egr, y_co2_out, y_o2_ox
    coef=coeff*10
    mw_fuel = config.res.fuel.thermo.mean_molecular_weight/1000 #kg/mol
    mw_oxidizer = config.res.ox.thermo.mean_molecular_weight/1000 #kg/mol
    mw_egr = config.res.egr.thermo.mean_molecular_weight/1000 #kg/mol

    Sx=2.0
    n_mole_ox = 4.76
    
    if(config.egr_unit=='mass'):
        X_egr = config.egr_rate/mw_egr
    elif(config.egr_unit=='mol'):
        X_egr = config.egr_rate
    elif(config.egr_unit=='vol'):
        X_egr = config.egr_rate*config.res.egr.thermo.density_mole

    X_o2_egr = config.res.egr.thermo['O2'].X*X_egr
    X_egr_no_o2 = (1-config.res.egr.thermo['O2'].X)*X_egr
    X_o2_air = (Sx*(1-X_egr_no_o2)/(phi)-X_o2_egr*(1+Sx/phi))/(1+n_mole_ox*Sx/(phi)) 
    X_air = n_mole_ox*X_o2_air
    X_fuel = 1-X_air-X_egr
    #phi_cor = Sx*X_fuel/(X_o2_air+X_o2_egr)
    #print(phi,phi_cor)
    mdot_tot = X_fuel*coef*mw_fuel + X_air*coef*mw_oxidizer + X_egr*coef*mw_egr
    #print(config.res.egr.thermo['O2'].X,config.res.ox_gas['N2'].X)
    #print(('%10.4f %10.4f %10.4f %10.4f')%(moldot_fuel,moldot_o2,moldot_air,moldot_egr))
    return [X_fuel*coef*mw_fuel, X_air*coef*mw_oxidizer, X_egr*coef*mw_egr], mdot_tot#, phi_cor

class case:
    def __init__(self,fuel,ox,egr,Q,egr_rate,egr_unit,mode='mdot_fixed'):
        self.compo = self.Compo(fuel,ox,egr)
        self.res = self.Reservoirs()
        self.pow = Q
        self.egr_rate = egr_rate
        self.egr_unit = egr_unit
        self.mode = mode

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

def create_reservoir(content, scheme, T, P):
    gas = ct.Solution(scheme)
    gas.TPX = T, P, content
    return ct.Reservoir(gas)

def burned_gas(phi,p,compo='O2:1.,N2:3.76, CH4:0.5',scheme='gri30.xml',ignition=True):
    gas = ct.Solution(scheme)
    gas.TPX = 400.0,p,compo
    gas.set_equivalence_ratio(phi, 'CH4', 'O2:1.0, N2:3.76')
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
        inlet.mass_flow_rate = mdots[i]
        i+=1

def mixer(phi,config):
    gb = burned_gas(phi,config.res.ox.thermo.P,scheme=config.res.fuel.thermo.source,ignition=False)
    reactor,pressure_reg = build_reactor(gb,volume=100000.0)
    #create the reactor network
    
    sim=ct.ReactorNet([reactor])
    sim.rtol = 1e-10
    sim.atol = 1e-10
    sim.initialize()
    #set the mass flow controllers according to phi and egr rate asked in the config
    mdots,mdot_tot = compute_mdots(phi, config, reactor)
    mfcs = init_reservoirs_mdots([config.res.fuel, config.res.ox,config.res.egr],mdots,reactor)## ajouter la detection du nombre de reservoir dans la config
    #compute steady state
    sim.advance_to_steady_state()
    return reactor.T, reactor.thermo.P, reactor.thermo.X

def reactor_0D(phi,config,power_regulation):
    #create the reactor and fill it with burned gas to ignite the mixture
    gb = burned_gas(phi,config.res.ox.thermo.P,ignition=True,scheme=config.res.fuel.thermo.source)
    reactor,pressure_reg = build_reactor(gb,volume=100000.0)
    #create the reactor network
    sim=ct.ReactorNet([reactor])

    sim.initialize()
    #sim.setTolerances(sim,rtol=1e-10, atol=1e-10)
    #set the mass flow controllers according to phi and egr rate asked in the config
    mdots,mdot_tot = compute_mdots(phi, config, reactor)
    mfcs = init_reservoirs_mdots([config.res.fuel, config.res.ox,config.res.egr],mdots,reactor)## ajouter la detection du nombre de reservoir dans la config
    #compute steady state
    sim.advance_to_steady_state()

    #create a loop to compute the mdots until the power is reached within a certain margin and limit the number of iterations to 10
    if(power_regulation):
        i=0
        currentpower = reactor.thermo.heat_release_rate*reactor.volume
        while (abs(currentpower - config.pow) > 10 and i<10):
            power_regulator = config.pow/currentpower
            mdots = [mdot*power_regulator for mdot in mdots]
            #print(mdots)
            edit_reservoirs_mdots(reactor, mdots)
            sim.reinitialize()
            try:
                sim.advance_to_steady_state()
            except:
                print('unknown error trying to reach steady state (power regulation)')
                break
            currentpower = reactor.thermo.heat_release_rate*reactor.volume
            i+=1
            print(('%2i %10.3f %10.3f %10.3f %10.4f %10.4f' % (i, mfcs[0].mass_flow_rate, mfcs[1].mass_flow_rate, mfcs[2].mass_flow_rate, reactor.thermo.heat_release_rate*reactor.volume, reactor.T)))
    
    return mfcs, reactor

def fresh_gas(phi,p,compo='O2:1.,N2:3.76, CH4:0.5',scheme='gri30.xml'):
    gas = ct.Solution(scheme)
    gas.TPX = 300.0,p,compo
    gas.set_equivalence_ratio(phi, 'CH4', 'O2:1.0, N2:3.76')
    return gas

def build_freeflame(mix):
    #mix = ct.Mixer(inlets=mfcs) #energy eq. is on by default
    flame = ct.FreeFlame(mix, width=0.2)
    return flame

def solve_flame(f):
    #################################################################
    # Iterations start here
    #################################################################
    loglevel  = 0                       # amount of diagnostic output (0 to 5)	    
    refine_grid = True                  # True to enable refinement, False to disable 	

    # first iteration
    ##################
    # No energy equilibrium activated (reduce the calculation)
    f.energy_enabled = False
    # Mesh refinement
    f.set_refine_criteria(ratio = 7.0, slope = 1, curve = 1)
    # Max number of times the Jacobian will be used before it must be re-evaluated
    f.set_max_jac_age(50, 50)
    #Set time steps whenever Newton convergence fails
    f.set_time_step(0.1e-06, [2, 5,10, 20, 80]) #s
    # Calculation
    f.solve(loglevel, refine_grid)

    # Second iteration
    #################
    # Energy equation activated
    f.energy_enabled = True
    # mesh refinement
    f.set_refine_criteria(ratio = 5.0, slope = 0.5, curve = 0.5)
    # Calculation
    f.solve(loglevel, refine_grid)
    
    # Third iteration
    #################
    # On raffine le maillage
    f.set_refine_criteria(ratio = 5.0, slope = 0.3, curve = 0.3)
    # Calculation
    f.solve(loglevel, refine_grid)
    #################
    # Quatri�me it�ration
    # Mesh refinement
    f.set_refine_criteria(ratio = 5.0, slope = 0.1, curve = 0.1)
    # Calculation
    f.solve(loglevel, refine_grid)
    # Fifth iteration
    #################
    # Mesh refinement
    f.set_refine_criteria(ratio = 5.0, slope = 0.05, curve = 0.05, prune = 0.01)
    # Calculation
    f.solve(loglevel, refine_grid)
    # Six iteration
    #################
    # Mesh refinement
    f.set_refine_criteria(ratio = 5.0, slope = 0.02, curve = 0.02, prune = 0.01)
    # Calculation
    f.solve(loglevel, refine_grid)
    
    return f

def compute_solutions_0D(config,phi_range,power_regulation=False,species = ['CH4','H2','O2','CO2','H2O']):
    data=np.zeros((len(phi_range),8))
    #create a dataframe naming colums with 'phi', 'T' and all the species in the list
    #then fill it with the values of phi, T and mole fractions of species using the concatenation of two dataframes, for each phi
    df =  pd.DataFrame(columns=['EGR','phi','T']+species)
    print(('%10s %10s %10s %10s %10s %10s %10s' % ('phi','fuel', 'air', 'egr', 'HRR', 'T', 'P'))) 
    for phi in phi_range:
 
        mfcs, reactor = reactor_0D(phi,config,power_regulation)
        moldot_tot = mfcs[0].mass_flow_rate/(config.res.fuel.thermo.mean_molecular_weight/1000)+ mfcs[1].mass_flow_rate/(config.res.ox.thermo.mean_molecular_weight/1000)+mfcs[2].mass_flow_rate/(config.res.egr.thermo.mean_molecular_weight/1000)
        #print(('%10.3f %10.3f %10.3f %10.3e %10.3f %10.3f' % (mfcs[0].mass_flow_rate/mdot_tot, mfcs[1].mass_flow_rate/mdot_tot, (mfcs[2].mass_flow_rate/mdot_tot), reactor.thermo.heat_release_rate, reactor.T, reactor.thermo.P)))
        print((' %10.3f %10.3f %10.3f %10.3f %10.3e %10.3f %10.3f' % (phi, (mfcs[0].mass_flow_rate/(config.res.fuel.thermo.mean_molecular_weight/1000))/moldot_tot, 
                                                                  (mfcs[1].mass_flow_rate/(config.res.ox.thermo.mean_molecular_weight/1000))/moldot_tot, 
                                                                  (mfcs[2].mass_flow_rate/(config.res.egr.thermo.mean_molecular_weight/1000))/moldot_tot, 
                                                                  reactor.thermo.heat_release_rate, reactor.T, reactor.thermo.P)))

        df = pd.concat([df, pd.DataFrame([[config.egr_rate, phi, reactor.T]+list(reactor.thermo[species].X)], columns=['EGR','phi','T']+species)]).astype(float)
        data[np.where(phi_range==phi), 0] = phi
        data[np.where(phi_range==phi), 1] = reactor.T
        data[np.where(phi_range==phi), 2] = reactor.thermo.heat_release_rate
        data[np.where(phi_range==phi), 3:] = reactor.thermo['CH4','H2','O2','CO2','H2O'].X
    
    return reactor, data, df

def compute_solutions_1D(config,phi_range,power_regulation=False,species = ['CH4','H2','O2','CO2','H2O']):
    #mdots,mdot_tot = compute_mdots(1.0, config)
    data=np.zeros((len(phi_range),8))
    #create a dataframe naming colums with 'phi', 'T' and all the species in the list
    #then fill it with the values of phi, T and mole fractions of species using the concatenation of two dataframes, for each phi
    df =  pd.DataFrame(columns=['EGR','phi','P','Tin','T','u']+species)
    
    for phi in phi_range:
        f = build_freeflame(fresh_gas(phi,config.res.ox.thermo.P,scheme=config.res.fuel.thermo.source))

        tol_ss = [1.0e-5, 1.0e-8]  # tolerance [rtol atol] for steady-state problem
        tol_ts = [1.0e-5, 1.0e-8]  # tolerance [rtol atol] for time stepping

        f.flame.set_steady_tolerances(default=tol_ss)
        f.flame.set_transient_tolerances(default=tol_ts)

        T,P,X = mixer(phi, config)

        f.inlet.T = T
        f.P = P
        f.inlet.X = X
        print('flame',f.inlet.X[13])

        f = solve_flame(f)
        index = [f.gas.species_index(specie) for specie in species]
        #print(list(f.X[index,-1]))
        df = pd.concat([df, pd.DataFrame([[config.egr_rate, phi, P, T, f.T[-1], f.velocity[0]]+list(f.X[index,-1])], columns=['EGR','phi','P','Tin','T','u']+species)]).astype(float) #+list(f.X[index][-1])] #
        data[np.where(phi_range==phi), 0] = phi
        data[np.where(phi_range==phi), 1] = f.T[-1]
        data[np.where(phi_range==phi), 2] = f.velocity[0]
        data[np.where(phi_range==phi), 3:] = f.X[index,-1]
    return f, data, df

def print_reactor(df):
    print(df)

def subplot_data(x,y,xlabel,ylabel,legend=None,symbol='-x'):
    rcParams['figure.figsize'] = (14, 12)
    nData = np.shape(y)[1]
    for i in range(nData):
        plt.subplot(nData-2, 2, i+1)
        plt.plot(x, y[:,i],symbol,label=legend)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel[i])
        plt.grid()
        if legend != None:
            plt.legend()

def see_graphs(suptitle=''):
    plt.suptitle(suptitle)
    plt.show()