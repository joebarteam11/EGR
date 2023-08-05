from lib_egr_260 import *

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

    return df