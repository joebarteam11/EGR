import sys
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from matplotlib import *
import colorama

colorama.init()

def get_ox_mdot(egr_percentage, ox_compo, egr_compo, mode='mdot_fixed', unit='vol', phi=1.0):
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
        print(colorama.Fore.RED + 'The oxidizer provided does not contain pure O2 !')
    



get_ox_mdot.__doc__="Calcule un débit massique a imposer en focntion du pourcentage de egr qu'on souhaite ajouter Ce pourcentage peut etre défini comme un pourcentage de 'remplacement' d'air : on va réduired'autant le débit d'air pour conserver un débit massique total constant défini par le débit d'air sans egr à la richesse souhaitée (mode mdot_ox_fixed) Il peut aussi être considéré comme un pourcentage de egr supplémentaire, défini comme un ajout de X% du débit d'air à la richesse considérée (mode mdot_ox_free)"
'''
def get_fuel_mdot(h2_percentage, fuel_compo, mode='mdot_fixed', unit='vol', phi=1.0):
    try:
        if(phi<0.0 or phi>10.0):
            raise Exception('phi must be in range ]0.0 ; 10;0]')

        mw_air = fuel_compo.mean_molecular_weight/1000 #kg/mol
        rho_air = fuel_compo.density_mass #kg/m3
        mw_egr = fuel_compo.mean_molecular_weight/1000 #kg/mol
        rho_egr = fuel_compo.density_mass #kg/m3

        x_o2 = fuel_compo.X[fuel_compo.species_index('O2')]
        mdot_air = ((mw_air/x_o2)*2.0)/phi
    except:
        print('The oxidizer provided does not contain pure O2 !')
        
    #calcule un débit massique a imposer en focntion du pourcentage de egr qu'on souhaite ajouter
    #Ce pourcentage peut etre défini comme un pourcentage de "remplacement" d'air : on va réduire
    #d'autant le débit d'air pour conserver un débit massique total constant défini par le débit d'air
    # sans egr à la richesse souhaitée (mode mdot_ox_fixed)
    #Il peut aussi être considéré comme un pourcentage de egr supplémentaire, défini comme un ajout de X% 
    #du débit d'air à la richesse considérée (mode mdot_ox_free)
    try:
        if(unit=='vol'):
            mdot_air_vol = mdot_air/rho_air
            mdot_egr_vol = mdot_air_vol * h2_percentage #volumetric
            mdot_egr = mdot_egr_vol * rho_egr
        elif(unit=='mol'):
            mdot_air_mol = mdot_air/mw_air
            mdot_egr_mol = mdot_air_mol * h2_percentage #molar
            mdot_egr = mdot_egr_mol * mw_egr
        elif(unit=='mass'): 
            mdot_egr = mdot_air * h2_percentage #mass
        else :
            raise Exception('Unknown unit')
        #print('Imposed egr massflow rate: '+ str(mdot_egr)+ '[kg/s]')
    except:
        print("Unknown unit '"+ unit +"', please choose between ['mol','mass','vol']")

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
'''

class case:
    def __init__(self,fuel,ox,egr,Q,egr_rate,egr_unit,mode='mdot_fixed'):
        self.compo = self.Compo(fuel,ox,egr)
        self.res = self.Reservoirs()
        self.Q = Q
        self.egr_rate = egr_rate
        self.egr_unit = egr_unit
        self.mode = mode
    #def __setattr__(self, res_fuel, res_ox, res_egr):
    #    self.res = self.Reservoirs(res_fuel, res_ox, res_egr)

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
            self.fuel_gas = self.Contents()
            self.ox_gas = self.Contents()
            self.egr_gas = self.Contents()

        class Contents:
            def __init__(self):
                self.contents = None


def compute_mdots(config,phi):
    coef=0.05
    mdot_fuel = 0.016
    mdot_air, mdot_egr=get_ox_mdot(config.egr_rate,
                                   config.res.ox_gas,
                                   config.res.egr_gas,
                                   unit=config.egr_unit,
                                   phi=phi
                                  )
    return mdot_fuel*coef, mdot_air*coef, mdot_egr*coef

def create_reservoir(content, scheme, T, P):
    gas = ct.Solution(scheme)
    gas.TPX = T, P, content
    stream = ct.Quantity(gas)
    return gas,ct.Reservoir(gas)#,stream

def burned_gas(compo='O2:1., CH4:0.5',scheme='gri30.xml'):
    gas = ct.Solution(scheme)
    gas.TPX = 300.0,1*100000,compo
    gas.equilibrate("HP")
    return gas

def build_reactor(mixture):
    mix = ct.IdealGasReactor(mixture, energy='on') #energy eq. is on by default
    exhaust = ct.Reservoir(mixture)
    pressure_regulator = ct.Valve(mix, exhaust, K=10.0)
    return mix

def set_reservoirs_mdots(amont, mdots, aval):
    mfc=[]
    i=0
    for res in amont:
        mfc.append(ct.MassFlowController(res, aval, mdot=mdots[i]))
        i+=1
    return mfc

def compute_solutions(config,phi_range,print_report=False,real_egr=False):
    data=np.zeros((len(phi_range),4))
    phi_bilger = np.zeros(len(phi_range))
    for phi in phi_range:
        gb = burned_gas()
        reactor = build_reactor(gb)
        mdots = compute_mdots(config,phi)
        if real_egr:
            mfcs = set_reservoirs_mdots([config.res.fuel, config.res.ox, reactor],mdots, reactor)## ajouter la detection du nombre de reservoir dans la config
        else:
            mfcs = set_reservoirs_mdots([config.res.fuel, config.res.ox,config.res.egr],mdots,reactor)## ajouter la detection du nombre de reservoir dans la config
        sim=ct.ReactorNet([reactor])
        sim.advance_to_steady_state()
        '''
        t=0.0
        for n in range(nit):
            tres = reactor.mass/sum(mfc.mass_flow_rate(t) for mfc in mfcs)
            t += tres
            sim.advance(t)
            if print_mdot :
                print(phi,reactor.T)
        '''
        phi_bilger[np.where(phi_range==phi)] = reactor.thermo.equivalence_ratio()
        data[np.where(phi_range==phi), 0] = reactor.T
        data[np.where(phi_range==phi), 1] = reactor.thermo.heat_release_rate
        data[np.where(phi_range==phi), 2:] = reactor.thermo['O2', 'CO2'].Y
    print(str(data[:, 0]))
    return phi_bilger,reactor, data

def print_reactor(reactor):
    print(reactor.thermo.report())

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