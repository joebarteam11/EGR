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
#from tqdm import tqdm
import warnings
import re

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

class FlameExtinguished(Exception):
    pass

class case:
    def __init__(self,fuels,Tinit_fuel,Pinit_fuel,ox,Tinit_ox,Pinit_ox,egr,Tinit_egr,Pinit_egr,phi_range,fuelblend_range,egr_range,fuelblend_unit,egr_unit,scheme,transport,isARC):
        self.compo = self.Compo(fuels,ox,egr)
        self.res = self.Reservoirs()
        self.gas = self.Gas(fuels,ox,egr)
        
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
        self.fuelblend_range = fuelblend_range
        self.egr_range = egr_range
        self.fuelblend_unit = fuelblend_unit
        self.egr_unit = egr_unit
        self.scheme = scheme
        self.transport = transport
        self.isARC = True if isARC == 'ARC' else False

    class Gas:
        def __init__(self,fuels,ox,egr):
            self.fuels = fuels
            self.ox = ox
            self.egr = egr

    class Compo:
        def __init__(self,fuels,ox,egr):
            self.fuels = fuels
            self.ox = ox
            self.egr = egr
    
    class Reservoirs:
        def __init__(self):
            self.fuels = None
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

def compute_mdots(config,phi,egr_rate,blend_ratio,coef=50.0,return_unit='mass',egr_def='%egr/%fuel'): 
    mw_fuel = config.gas.fuels.mean_molecular_weight/1000 #kg/mol
    mw_o2 = config.gas.ox['O2'].molecular_weights/1000.0#kg/mol
    mw_n2 = config.gas.ox['N2'].molecular_weights/1000.0#kg/mol
    mw_oxidizer = config.gas.ox.mean_molecular_weight/1000 #kg/mol
    mw_egr = config.gas.egr.mean_molecular_weight/1000 #kg/mol
    mol_weights=[mw_fuel,mw_oxidizer,mw_egr]

    warnings.warn('compute_mdots() is setted for methane-air mixtures only, change Sx, Sy and n_mole_ox for others reactants', UserWarning)
    sx=[]
    #check if there is carbon in the fuel, then find the number of carbon atoms in the fuel (the integer just after the letter C) etc...
    for fuel in config.compo.fuels:
        fuel=fuel.split(':')[0]
        nC=0
        nH=0
        nO=0
        if(fuel.find('C')!=-1):
            nC=re.findall(r'C(\d+)', fuel) if re.findall(r'C(\d)+', fuel) != [] else [1]
            nC = list(map(float, nC))[0]
        if(fuel.find('H')!=-1):
            nH=re.findall(r'H(\d+)', fuel) if re.findall(r'H(\d)+', fuel) != [] else [1]
            nH = list(map(float, nH))[0]
        if(fuel.find('O')!=-1):
            nO=re.findall(r'O(\d+)', fuel) if re.findall(r'O(\d)+', fuel) != [] else [1]
            nO = list(map(float, nO))[0]

        sx.append(nC+nH/4-nO/2)

    nO2=re.search(r'O2:(\d+\.\d+)', config.compo.ox).group(1)
    nO2=float(nO2)
    nN2=re.search(r'N2:(\d+\.\d+)', config.compo.ox).group(1)
    nN2=float(nN2)

    n_mole_ox=sum([nO2,nN2])
    ox_n2_ratio=nN2/nO2
    mw_ox=(nO2*mw_o2+ox_n2_ratio*mw_n2)[0]

    Sx=sx[0]
    Sy=sx[0]*mw_ox/mw_fuel

    if(blend_ratio is not None and blend_ratio != 0):
        if(config.fuelblend_unit=='mole'):
            Sx=(1-blend_ratio)*sx[0]+(blend_ratio)*sx[1]

        elif(config.fuelblend_unit=='mass'):
            Sy=((1-blend_ratio)*sx[0]+(blend_ratio)*sx[1])*mw_ox/mw_fuel

    print('mw_ox :',mw_ox)
    print('Sx :',Sx)
    print('Sy :',Sy)
    print('n_mole_ox :',n_mole_ox)
    
    if(config.egr_unit=='mass'):
        fuel_mass=phi/(Sy+phi)
        air_mass=1-fuel_mass

        if(egr_def=='%egr/%ox'):
            egr_comparison=air_mass
        elif(egr_def=='%egr/%fuel'):
            egr_comparison=fuel_mass
        elif(egr_def=='%egr/all'):
            egr_comparison=air_mass+fuel_mass
        elif(egr_def=='custom'):
            egr_comparison=None

        egr_mass=(egr_rate/(1-egr_rate))*(egr_comparison) #+fuel_mass depending on the definition you want
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

        if(egr_def=='%egr/%ox'):
            egr_comparison=air_mol
        elif(egr_def=='%egr/%fuel'):
            egr_comparison=fuel_mol
        elif(egr_def=='%egr/all'):
            egr_comparison=air_mol+fuel_mol

        egr_mol=(egr_rate/(1-egr_rate))*(egr_comparison) # depending on the definition you want
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
        fuel_vol = fuel_mol * config.gas.fuels.volume_mole
        air_vol = (1-fuel_mol) * config.gas.ox.volume_mole

        if(egr_def=='%egr/%ox'):
            egr_comparison=air_vol
        elif(egr_def=='%egr/%fuel'):
            egr_comparison=fuel_vol
        elif(egr_def=='%egr/all'):
            egr_comparison=air_vol+fuel_vol

        egr_vol=(egr_rate/(1-egr_rate))*(egr_comparison) #+fuel_vol depending on the definition you want
        mdots=[fuel_vol*coef*config.gas.fuels.density_mass, air_vol*coef*config.gas.ox.density_mass, egr_vol*coef*config.gas.egr.density_mass]
        
        if(return_unit=='mass'):
            return mdots, sum(mdots)
        elif(return_unit=='mole'):
            warnings.warn('not tested yet, use mass instead', UserWarning)
            moldots=[m/n for m, n in zip(mdots,mol_weights)]
            return moldots, sum(moldots)

    else:
        raise ValueError('egr_unit must be mass, mole or vol')
    
def create_reservoir(config,content, T, P,blend_ratio=None,scheme=None):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
    if(config.isARC):
        #replace extension from scheme variable
        ct.compile_fortran(config.scheme.split('.')[0]+'.f90')
        
    if(scheme is None):
        scheme = config.scheme
        
    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        gas=ct.Solution(scheme, transport_model=config.transport)
    else:
        gas=ct.Solution(scheme)

    if type(content) is not list :
        gas.TPX = T, P, content
        return ct.Reservoir(gas), gas
    
    elif type(content) is list and blend_ratio is None:
        raise ValueError('blend_ratio must be set if content is a list')
    
    # elif type(content) is not list and blend_ratio is not None:
    #     raise ValueError('content must be a list if blend_ratio is set')

    elif type(content) is list and blend_ratio < 0.00001 :
        config.compo.fuel = content[0]
        gas.TPX = T, P, config.compo.fuel
        return ct.Reservoir(gas), gas
    
    elif len(content) <= 2 and blend_ratio > 0.00001:
        if(config.fuelblend_unit == 'mole'):
            mixed_content = ""
            for i,cont in enumerate(content):
                splited_content = cont.split(':')
                if(i==0):
                    mixed_content += (splited_content[0]+':'+str(float(splited_content[1])*(1-blend_ratio))+' ')
                else:    
                    mixed_content += (splited_content[0]+':'+str(float(splited_content[1])*(blend_ratio)))
            config.compo.fuel = mixed_content
            gas.TPX = T, P, config.compo.fuel

        elif(config.fuelblend_unit == 'mass'):
            mixed_content = ""
            for i,cont in enumerate(content):
                splited_content = cont.split(':')
                if(i==0):
                    mixed_content += (splited_content[0]+':'+str(float(splited_content[1])*(1-blend_ratio)))
                else:    
                    mixed_content += (splited_content[0]+':'+str(float(splited_content[1])*(blend_ratio)))
            config.compo.fuel = mixed_content
            gas.TPY = T, P, config.compo.fuel


        else:
            raise ValueError('fuelblend_unit must be "mass" or "mole"')

        return ct.Reservoir(gas), gas
    
    else:
        raise ValueError('Only two streams can be mixed in varaying proportions (but each stream can be a mixture)')

def burned_gas(phi,config,egr_rate,ignition=True):
    if(config.isARC):
        #remove extension from sheme variable
        ct.compile_fortran(config.scheme.split('.')[0]+'.f90')
    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        gas=ct.Solution(config.scheme, transport_model=config.transport)
    else:
        gas=ct.Solution(config.scheme)
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
    if(version.parse(ct.__version__) >= version.parse("2.5.0")):
        for inlet in reactor.inlets:
            inlet.mass_flow_rate = mdots[i]
            i+=1
    else:
        for inlet in reactor.inlets:
            inlet.set_mass_flow_rate(mdots[i])
            i+=1

def reactor_0D(phi,config,egr_rate,real_egr,max_residence_time=1.0,steady_state_only=False):
    #create the reactor and fill it with burned gas to ignite the mixture
    gb = burned_gas(phi,config,egr_rate) 
    reactor,pressure_reg = build_reactor(gb,volume=10000.0)#high volume to ensure that the residence time is high enough
    reservoirs=[config.res.fuel, config.res.ox, config.res.egr]
    #create the reactor network
    sim=ct.ReactorNet([reactor])
    sim.reinitialize()
    sim.rtol = 1e-8
    sim.atol = 1e-16
    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        sim.max_steps = 100000 
    #set the mass flow controllers according to phi and egr rate asked in the config
    def mdot_tot_required(residence_time):
        return reactor.mass / residence_time
    
    mdots,mdot_tot = compute_mdots(config, egr_rate, phi)
    mfcs = None
        #create the mass flow controllers according to the number of reservoirs
    if(real_egr):
        mfcs = init_reservoirs_mdots([config.res.fuel, config.res.ox,reactor],mdots,reactor)## ajouter la detection du nombre de reservoir dans la config
    else:
        mfcs = init_reservoirs_mdots(reservoirs,mdots,reactor) ## ajouter la detection du nombre de reservoir dans la config

    states = ct.SolutionArray(gb, extra=['tres'])
    for t in np.linspace(1e-3,max_residence_time, 100):
        #mdots,mdot_tot = compute_mdots(config, egr_rate, phi)
        coef_correction = mdot_tot_required(t)/mdot_tot
        mdots_corrected = [mdot*coef_correction for mdot in mdots]
        new_mdot_tot = sum(mdots_corrected)
        print('new mdot tot :',new_mdot_tot)
        print('reactor mass :',reactor.mass)

        edit_reservoirs_mdots(reactor,mdots_corrected)
        
        print('Real residence time :',reactor.mass / new_mdot_tot,'Required residence time :',t)


        sim.reinitialize()
        sim.advance_to_steady_state()
        #time_mesh = np.linspace(0.0, residence_time, 100)
        #for t in np.linspace(0,residence_time, 100):
        #    sim.advance(t)
        states.append(reactor.thermo.state, tres=t)
    return mfcs, mdot_tot, reactor, states

def compute_equilibrium(config,phi,tin,pin,egr,fb,species = ['CH4','H2','O2','CO','CO2','H2O']):
    #create a dataframe naming colums with 'phi', 'T' and all the species in the list
    #then fill it with the values of phi, T and mole fractions of species using the concatenation of two dataframes, for each phi
    df =  pd.DataFrame(columns=['EGR','FB','phi','T','P']+species)
    
    # for egr in config.egr_range:
    #     for phi in config.phi_range:
    _,config.gas.fuels = create_reservoir(config,config.compo.fuels, tin[0], pin[0],blend_ratio=fb)
    _,config.gas.ox = create_reservoir(config,config.compo.ox, tin[1], pin[1],scheme='air.xml')
    _,config.gas.egr = create_reservoir(config,config.compo.egr, tin[2], pin[2])
            
    bg = burned_gas(phi,config,egr,ignition=True)
    
    df = pd.concat([df, pd.DataFrame([[egr,fb, phi, bg.T, bg.P]+list(bg[species].X)], columns=['EGR','FB','phi','T','P']+species)]).astype(float)
    try:
        update()
    except:
        pass

    return df

def compute_solutions_0D(config,phi,tin,pin,egr,fb,real_egr=False,species = ['CH4','H2','O2','CO','CO2','H2O'],res_time=1.0):
    #create a dataframe naming colums with 'phi', 'T' and all the species in the list
    #then fill it with the values of phi, T and mole fractions of species using the concatenation of two dataframes, for each phi
    df =  pd.DataFrame(columns=['tres','EGR','FB','phi','T','P']+species)
    if(version.parse(ct.__version__) >= version.parse("2.4.0")):
        print(('%10s %10s %10s %10s %10s %10s %10s %10s' % ('mdot_tot','phi','Xfuel', 'Xair', 'Xegr', 'HRR', 'T', 'P'))) 

    # for egr in config.egr_range:
    #     for phi in config.phi_range:
    config.res.fuel,config.gas.fuels = create_reservoir(config,config.compo.fuels, tin[0], pin[0],blend_ratio=fb)
    config.res.ox,config.gas.ox = create_reservoir(config,config.compo.ox, tin[1], pin[1],scheme='air.xml')
    config.res.egr,config.gas.egr = create_reservoir(config,config.compo.egr, tin[2], pin[2])
    
    mfcs, mdot_tot, reactor, states = reactor_0D(phi,config,egr,real_egr,res_time,steady_state_only=False)
    #only with cantera >= 2.5
    if(version.parse(ct.__version__) >= version.parse("2.4.0")):
        moldot_tot = mfcs[0].mass_flow_rate/(config.gas.fuels.mean_molecular_weight/1000)+ mfcs[1].mass_flow_rate/(config.gas.ox.mean_molecular_weight/1000)+mfcs[2].mass_flow_rate/(config.gas.egr.mean_molecular_weight/1000)
        #print(('%10.3f %10.3f %10.3f %10.3e %10.3f %10.3f' % (mfcs[0].mass_flow_rate/mdot_tot, mfcs[1].mass_flow_rate/mdot_tot, (mfcs[2].mass_flow_rate/mdot_tot), reactor.thermo.heat_release_rate, reactor.T, reactor.thermo.P)))
        print((' %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f' % (mdot_tot,phi, #reactor.thermo['CH4'].X, reactor.thermo['O2'].X+reactor.thermo['N2'].X, reactor.thermo['CO2'].X,
                                                            #only with cantera >= 2.5
                                                            (mfcs[0].mass_flow_rate/(config.gas.fuels.mean_molecular_weight/1000))/moldot_tot,
                                                            (mfcs[1].mass_flow_rate/(config.gas.ox.mean_molecular_weight/1000))/moldot_tot, 
                                                            (mfcs[2].mass_flow_rate/(config.gas.egr.mean_molecular_weight/1000))/moldot_tot, 
                                                            reactor.thermo.heat_release_rate, #only with cantera >= 2.5
                                                            #'N/A',
                                                            reactor.T, reactor.thermo.P)))

    #df = pd.concat([df, pd.DataFrame([[egr, phi, reactor.T, reactor.thermo.P]+list(reactor.thermo[species].X)], columns=['EGR','phi','T','P']+species)]).astype(float)
    if(states is not None):
        index = [reactor.thermo.species_index(specie) for specie in species]
        for i in range(len(states.tres)):
            df = pd.concat([df, pd.DataFrame([[states.tres[i]]+[egr,fb, phi, states.T[i], reactor.thermo.P]+list(states.X[i,index])], columns=['tres']+['EGR','FB','phi','T','P']+species)]).astype(float)
    else:
        df = pd.concat([df, pd.DataFrame([[egr,fb, phi, reactor.T, reactor.thermo.P]+list(reactor.thermo[species].X)], columns=['EGR','FB','phi','T','P']+species)]).astype(float)

    try:
        update()
    except:
        pass

    return reactor, df

def fresh_gas(phi,config,egr_rate,Tmix,Pmix):
    if(config.isARC):
        #remove extension from sheme variable
        ct.compile_fortran(config.scheme.split('.')[0]+'.f90')
    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        gas=ct.Solution(config.scheme, transport_model=config.transport)
    else:
        gas=ct.Solution(config.scheme)
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

def mixer(config,phi,egr,fb,real_egr=False,T_reinj=None):
    if(config.isARC):
        #remove extension from sheme variable
        ct.compile_fortran(config.scheme.split('.')[0]+'.f90')
    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        gas=ct.Solution(config.scheme, transport_model=config.transport)
    else:
        gas=ct.Solution(config.scheme)
    mdots,mdot_tot = compute_mdots(config, phi, egr, fb, return_unit='mass')

    fuel = ct.Quantity(gas,constant='HP')
    fuel.TPX = config.gas.fuels.T,config.gas.fuels.P,config.compo.fuel
    fuel.mass = mdots[0]/mdot_tot
    ox = ct.Quantity(gas,constant='HP')
    ox.TPX = config.gas.ox.T,config.gas.ox.P,config.compo.ox
    ox.mass = mdots[1]/mdot_tot
    
    if(real_egr):
        dilutent = ct.Quantity(config.gas.egr,constant='HP')
        if(T_reinj is not None):
            dilutent.TP = T_reinj,config.gas.egr.P
        else:
            dilutent.TP = config.gas.egr.T,config.gas.egr.P
    else:
        dilutent = ct.Quantity(gas,constant='HP')
        dilutent.TPX = config.gas.egr.T,config.gas.egr.P,config.compo.egr
    dilutent.mass = mdots[2]/mdot_tot

    mix = fuel+ox+dilutent
    #print(mix.report())
    return mix.T, mix.P, mix.X

def apply_egr_to_inlet(f,config,phi,egr,fb,dry=False,T_reinj=None):
    config.gas.egr.X = f.X[:,-1]
    #print('EGR composition',config.gas.egr.X)

    if(dry):
        index = f.gas.species_index('H2O')
        dry_egr = config.gas.egr.X
        dry_egr[index] = 0.0
        coef = 1/sum(dry_egr)
        dry_egr = [coef*i for i in dry_egr]
        config.gas.egr.X = dry_egr
        #print('DRY EGR composition',config.gas.egr.X)
        #print('sum dry egr',sum(config.gas.egr.X[:]))

    T,P,X = mixer(config,phi,egr,fb,real_egr=True,T_reinj=T_reinj)
    f.inlet.X = X
    f.inlet.T = T
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

def compute_solutions_1D(config,phi,tin,pin,egr,fb,restart_rate,real_egr,dry=False,T_reinj=None,species = ['CH4','O2','CO2','H2O'],vars=['EGR','FB','phi','AF','P','Tin','T','u','dF']):
    path = os.getcwd()+'/src'
    dfs=[]
    vartosave = vars+species
    df =  pd.DataFrame(columns=vartosave)

    #create the gas object containing the mixture of fuel, ox and egr and all thermo data
    _, config.gas.fuels = create_reservoir(config,config.compo.fuels,tin[0], pin[0],blend_ratio=fb)
    _, config.gas.ox = create_reservoir(config,config.compo.ox, tin[1], pin[1],scheme='air.xml')
    _, config.gas.egr = create_reservoir(config,config.compo.egr, tin[2], pin[2])

    #create a dataframe naming colums with 'phi', 'T' and all the species in the list
    #then fill it with the values of phi, T and mole fractions of species using the concatenation of two dataframes, for each phi
    #vartosave = vars+species
    #df =  pd.DataFrame(columns=vartosave)

    #get the temperature and pressure of the mixture according to phi and egr rate
    T,P,X = mixer(config,phi,egr,fb)
    f = build_freeflame(fresh_gas(phi,config,egr,T,P))
    #print(f.transport_model)
    tol_ss = [2.0e-7, 1.0e-7]  # tolerance [rtol atol] for steady-state problem
    tol_ts = [2.0e-7, 1.0e-7]  # tolerance [rtol atol] for time stepping

    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    f.transport_model = config.transport

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
            # f.read_hdf(flametitle)
            a=1+1
            #f.set_initial_guess(data=flametitle)
        except:
            raise Exception('Cannot restore flame from file '+flametitle)
        flametitle = path+'/data/'+'egr'+str(round(egr,1))+'_phi'+str(round(phi,2))+'_T'+str(round(T,0))+'_P'+str(round(P/100000,0))+'_ARC_QC16.h5'

    # print(f.X[:,-1])
    # config.gas.egr.X = f.X[:,-1]
    # _,_,X = mixer(phi,config,egr,real_egr=True)
    # print('Inlet composition',X)
    if(real_egr):
        f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
    else:
        f.inlet.T = T
        f.P = P
        f.inlet.X = X
    #print('Inlet composition',f.inlet.X)
    flame = solve_flame(f,flametitle,config,phi,egr,fb,real_egr=real_egr,dry=dry,T_reinj=T_reinj)
    if(version.parse(ct.__version__) == version.parse("2.3.0")):
        #flame.write_AVBP('solut_avbp_ARC.csv')
        species_names = flame.gas.species_names
        net_prod_rate = flame.net_production_rates
        for i,specie in enumerate(species_names):
            print(f"{specie}: {np.max(np.abs(net_prod_rate[i])):.6e} kmol/m^3/s")
        #print the maximum 'w_CH4' value
        #print('Maximum w_CH4',np.max(flame.w[0]))

    #restart_rate = egr
    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        SL0=f.velocity[0]
    else:
        SL0=f.u[0]

    index = [f.gas.species_index(specie) for specie in species]
    df = pd.concat([df, pd.DataFrame([[egr,fb, phi,1/phi, P, T, f.T[-1],SL0,flamme_thickness(f)]+list(f.X[index,-1])], columns=vartosave)]).astype(float) #+list(f.X[index][-1])] #
    #dfs = pd.concat([df],axis=0)
        
    return df

def solve_flame(f,flametitle,config,phi,egr,fb,real_egr=False,dry=False,T_reinj=None):
    #warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
    #################################################################
    # Iterations start here
    #################################################################
    
    residuals = []
    index = f.gas.species_index('CO2')
    last_residual = 1.0 
    i=0
    verbose = 0
    loglevel  = 0                      # amount of diagnostic output (0 to 5)	    
    refine_grid = True                  # True to enable refinement, False to disable 	
    f.max_time_step_count=25000
    f.max_grid_points=1000
    
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

    while(True or i<5):
        #print('While loop starts')
        # Calculation
        if(verbose>0):
            print('1st iteration...')
        try:
            if(verbose>1):
                print('Inlet composition before f.solve',f.inlet.X)
            f.solve(loglevel, refine_grid)
            
            if(real_egr):
                XCO2_1 = f.X[index,-1]
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
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
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
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
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
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
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
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
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    print('EGR applied to inlet')
            # try:
            #     f.write_hdf(flametitle)
            # except:
            #     f.save(flametitle[:-3]+'.yaml')
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
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
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
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    print('EGR applied to inlet')
            #f.save(flametitle)
            # try:
            #     f.write_hdf(flametitle)
            # except:
            #     f.save(flametitle[:-3]+'.yaml')
        except FlameExtinguished:
            print('Flame extinguished')
            
        except ct.CanteraError as e:
            print('Error occurred while solving:', e)
        if(real_egr):
            last_residual = abs(residuals[-1]-residuals[-2])
            #print('last residual',last_residual)
            #print('Residuals',residuals)
        else:
            break
        #print('last residual',last_residual)
        if(last_residual<1e-11):
            #print('BREAK HERE')
            break
        i+=1
        if(i>5):
            raise Exception('Too many iterations of solve_flame, cannot converge')
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