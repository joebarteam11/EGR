from lib_egr_260 import *

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
    if(version.parse(ct.__version__) >= version.parse("2.5.0")):
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

    return reactor, df

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