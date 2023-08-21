import sys,os
import warnings
import re
import logging
import hashlib  
import numpy as np
import cantera as ct
from mpi4py import MPI
import pandas as pd
from packaging import version
from math import floor

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

comm = MPI.COMM_WORLD
computelog = logging.getLogger(str(comm.Get_rank()))

hl = logging.FileHandler(filename=computelog.name+".log",mode='w')
format = logging.Formatter('%(message)s')
hl.setFormatter(format)
computelog.setLevel(logging.INFO)
computelog.addHandler(hl)


def logprint(message_to_log,priority="info",file=None):

    if file is not None:
        print(message_to_log)
        sys.stdout.flush()

    if priority == "info":
        computelog.info(message_to_log)
    elif priority == "warning":
        computelog.warning(message_to_log)

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

def init_cases(config):
    if(config.isARC):
        ct.compile_fortran(config.scheme.replace('.cti','.f90'))

def stoechiometric_ratios(config):
    mw_o2 = config.gas.ox['O2'].molecular_weights/1000.0 #kg/mol
    mw_n2 = config.gas.ox['N2'].molecular_weights/1000.0 #kg/mol
    sx = []
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
        try:
            nN2=re.search(r'N2:(\d+\.\d+)', config.compo.ox).group(1)
            nN2=float(nN2)
        except:
            nN2=0.0

        n_mole_ox=sum([nO2,nN2])
        ox_n2_ratio=nN2/nO2
        mw_ox=(nO2*mw_o2+ox_n2_ratio*mw_n2)[0]
        
    return sx,n_mole_ox,mw_ox

def compute_mdots(config,phi,egr_rate,blend_ratio,coef=50.0,return_unit='mass',egr_def='%egr/%fuel'): 
    mw_fuel = config.gas.fuels.mean_molecular_weight/1000 #kg/mol
    mw_oxidizer = config.gas.ox.mean_molecular_weight/1000 #kg/mol
    mw_egr = config.gas.egr.mean_molecular_weight/1000 #kg/mol
    mol_weights=[mw_fuel,mw_oxidizer,mw_egr]

    sx,n_mole_ox,mw_ox = stoechiometric_ratios(config) # get the molar stoechiometric ratio for each fuel in the blend

    Sx=sx[0]
    Sy=sx[0]*mw_ox/mw_fuel

    if(blend_ratio is not None and blend_ratio != 0):
        if(config.fuelblend_unit=='mole'):
            Sx=(1-blend_ratio)*sx[0]+(blend_ratio)*sx[1]

        elif(config.fuelblend_unit=='mass'):
            Sy=((1-blend_ratio)*sx[0]+(blend_ratio)*sx[1])*mw_ox/mw_fuel
    
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
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        gas=ct.Solution(config.scheme, transport_model=config.transport)
    else:
        gas=ct.Solution(config.scheme)

    gas.TP = config.gas.ox.T,config.gas.ox.P

    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        gas.set_equivalence_ratio(phi, config.compo.fuel, config.compo.ox,
            #additional params, only compatible with cantera >= 2.5
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

def mixer(config,phi,egr,fb,real_egr=False,T_reinj=None):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
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
    logprint("mixer composition")
    logprint(mix.report())

    return mix.T, mix.P, mix.X

def apply_egr_to_inlet(f,config,phi,egr,fb,dry=False,T_reinj=None):
    config.gas.egr.X = f.X[:,-1]
    logprint("-----------------------------------------")
    logprint('EGR composition BEFORE drying operation',)
    species_names = config.gas.egr.species_names
    for i,specie in enumerate(species_names):
            logprint(f"{specie}: {config.gas.egr.X[i]:.6e} X[mol]")


    if(dry):
        index = f.gas.species_index('H2O')
        dry_egr = config.gas.egr.X
        dry_egr[index] = 0.0
        coef = 1/sum(dry_egr)
        dry_egr = [coef*i for i in dry_egr]
        config.gas.egr.X = dry_egr
        
        logprint('EGR composition AFTER drying operation',)
        for i,specie in enumerate(species_names):
                logprint(f"{specie}: {config.gas.egr.X[i]:.6e} X[mol]")

    logprint("-----------------------------------------")
    T,P,X = mixer(config,phi,egr,fb,real_egr=True,T_reinj=T_reinj)
    f.inlet.X = X
    f.inlet.T = T

    return f

def generate_unique_filename(config,flame):
    # Convert parameters to strings for inclusion in the identifier
    compo_str = '_'.join([f"{specie}:{flame.inlet.X[i]:.2f}"  for i,specie in enumerate(flame.gas.species_names) if flame.inlet.X[i]>0])
    #oxidizer_str = '_'.join([f"{species}:{config.gas.ox.X[i]:.2f}" for i,species in enumerate(config.gas.ox.species_names) if config.gas.ox.X[i]>0])
    pressure_str = f"{floor(round(flame.P/1e5,3))}atm"
    temperature_str = f"{floor(round(flame.inlet.T/100,3))}K/100"

    # Combine the parameters into a single string
    parameter_str = f"{compo_str}_{pressure_str}_{temperature_str}_{config.scheme.split('/')[-1].split('.')[0]}"
    logprint(f"Flame name: {parameter_str}")
    # Generate a hash of the parameter string
    hash_object = hashlib.md5(parameter_str.encode())
    uid = hash_object.hexdigest()

    return parameter_str,uid