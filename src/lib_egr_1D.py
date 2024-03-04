from lib_egr_260 import *


def flame_saver(f,loglevel,flamename,hash,config):
    try:
        # if hash is an existing file 
        if os.path.isfile(hash):
            logprint('Flame file '+hash+' is already existant, will not overwrite')
        else:
            logprint('Saving flame to file '+hash)
            f.save(hash, loglevel=loglevel, description=flamename + "\n Cantera version "+ct.__version__+' with '+config.transport+' transport model and mechanism '+config.scheme )
    except:
        logprint('Cannot save flame to file '+flamename)

def flame_saver_csv(f,config,phi,egr,fb):
    if config.saveCSV:
        if not os.path.isdir(config.saveCSVpath):
            os.makedirs(config.saveCSVpath)
        
        flamename = config.saveCSVpath+'/T'+str(config.gas.fuels.T)+'K_P'+str(round(config.gas.fuels.P/1e5,2))+'bar_phi'+str(round(phi, 2))+'_egr'+str(egr)+'_fb'+str(fb)
        f.write_AVBP(flamename+'.csv',quiet=False)
        logprint('Flame csv saved : '+flamename)

        


def solve_flame(f,hash,flamename,config,phi,egr,fb,real_egr=False,dry=False,T_reinj=None,restore_success=False):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...   
    
    #################################################################
    # Iterations start here
    #################################################################
    maxegrate_iter = 50
    residuals = []
    index = f.gas.species_index('H2')
    last_residual = 1.0 
    i=0
    verbose = 1
    loglevel  = 0                      # amount of diagnostic output (0 to 5)	    
    refine_grid = 'refine' #True                  # True to enable refinement, False to disable 	
    f.max_time_step_count=75000
    f.max_grid_points=1000
    auto = False

    # first iteration
    ##################
    # No energy equilibrium activated (reduce the calculation)
    # f.energy_enabled = False
    # Mesh refinement
    # f.set_refine_criteria(ratio = 7.0, slope = 0.95, curve = 0.95)
    # Max number of times the Jacobian will be used before it must be re-evaluated
    f.set_max_jac_age(10, 10)
    #Set time steps whenever Newton convergence fails
    f.set_time_step(1e-08, [25, 40, 80, 140, 200, 350, 500, 700, 1000, 1300, 1700, 2000, 3000, 5000, 10000, 12000, 15000, 20000]) #s

    # criteria_list determines the number of iterations that you want and tells : ratio, slope, curve, prune and energy enabled
    criteria_list = [[7.0, 0.95, 0.95,0, False]] # First iteration criterias ( 0 prune means disable)
    criteria_list += [[7.0, 0.75, 0.75,0, True]] # Second iteration criterias
    criteria_list += [[5.0, 0.4, 0.4,0.03, True]] # Third iteration criterias
    criteria_list += [[5.0, 0.1, 0.1,0.01, True]] # Fourth iteration criterias
    criteria_list += [[5.0, 0.05, 0.05,0.01, True]] # Fifth iteration criterias
    criteria_list += [[4.0, 0.05, 0.05,0.01, True]] # Sixth iteration criterias
    last_iteration = [4.0, 0.025, 0.025,0.005, True] # 7th iteration criterias

    
    while(True or i<maxegrate_iter):
        # Calculation
        if(not restore_success):
            auto_success = False
            auto = True 
            # One loop of auto refinement
            i = " auto iteration"
            f,auto_success,XCO2_1 = flame_iteration(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,i)
    
            if(not auto_success):
                auto = False
                boucle_over_flame_iterations(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,criteria_list)
            else:
                logprint('Successful auto refinement, going to last_iteration')
        else:
            logprint('Successful restore, going to last_iteration')
        #end if not restore_success

        # # Last iteration
        f.energy_enabled = True
        f.set_refine_criteria(ratio = last_iteration[0], slope = last_iteration[1], curve = last_iteration[2], prune = last_iteration[3])
        
        i = " last iteration"
        f,sucess,XCO2_1 = flame_iteration(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,i)
        T_flamme = f.T[-1]
        if sucess:
            flame_saver(f,loglevel,flamename,hash,config)
            flame_saver_csv(f,config,phi,egr,fb)
        else:
            pass


        if(real_egr):
            last_residual = abs(residuals[-1]-residuals[-2])
        else:
            break
        if(last_residual<1e-9):
            T_flamme = f.T[-1]
            break
        i+=1
        if(i>maxegrate_iter):
            T_flamme = np.Inf
            logprint('Too many iterations of solve_flame for real_egr loop, for flame ',hash,' cannot converge')
            break
    return f,T_flamme

def compute_solutions_1D(config,phi,tin,pin,egr,fb,restart_rate,real_egr,tol_ss=[2.0e-5, 1.0e-9],tol_ts=[2.0e-5, 1.0e-9],dry=True,T_reinj=None,species = ['CH4','O2','CO2','H2O'],vars=['EGR','FB','phi','AF','P','Tin','T','u','dF','rhoGF','rhoGB']):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
    path = os.getcwd()+'/src'


    # tol_ss = [2.0e-5, 1.0e-9]  # tolerance [rtol atol] for steady-state problem
    # tol_ts = [2.0e-5, 1.0e-9]  # tolerance [rtol atol] for time stepping

    #create the gas object containing the mixture of fuel, ox and egr and all thermo data
    _, config.gas.fuels = create_reservoir(config,config.compo.fuels,tin[0], pin[0],blend_ratio=fb)
    _, config.gas.ox = create_reservoir(config,config.compo.ox, tin[1], pin[1],scheme='air.yaml')
    _, config.gas.egr = create_reservoir(config,config.compo.egr, tin[2], pin[2])

    #get the temperature and pressure of the mixture according to phi and egr rate
    T,P,X = mixer(config,phi,egr,fb)
    gas=fresh_gas(phi,config,egr,T,P)
    f = build_freeflame(gas,width=0.02)

    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    f.transport_model = config.transport

    f.inlet.T = T
    f.P = P
    f.inlet.X = X

    flamename,hash = generate_unique_filename(config,f)
    logprint('Flame hash: '+hash)

    flamefile = path+'/data/'+hash+'.yaml'
    restore_success = False
    if os.path.isfile(flamefile):
        logprint('Flame file '+flamefile+' is found')
        try:
            f.restore(flamefile, loglevel=1)
            restore_success =True
            logprint('Flame restored from file '+flamefile)
        except:
            logprint('Cannot restore flame from file '+flamefile)
            f.set_initial_guess()
            logprint('Flame initial guess set')
    else:
        f.set_initial_guess()
        logprint('Flame initial guess set')

    if(real_egr):
        f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
    else:
        f.inlet.T = T
        f.P = P
        f.inlet.X = X

    f.inlet.T = T
    f.P = P
    f.inlet.X = X

    # alpha,beta=alpha_beta(f)
    # print('Alpha, beta',alpha_beta(f))
    # logprint('Inlet composition (mol)',f.inlet.X)
    # logprint('Inlet composition (mass)',f.inlet.Y)

    f,T_flamme = solve_flame(f,flamefile,flamename,config,phi,egr,fb,real_egr=real_egr,dry=dry,T_reinj=T_reinj,restore_success=restore_success)

    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        SL0=f.velocity[0]
    else:
        SL0=f.u[0]
    
    logprint('SL0 :',SL0)

    omega0=compute_omega0(f)
    #f.write_AVBP(path+'/'+'CH4_phi073'+'.csv')
    #print("Mean molecular weight: ",f.gas.mean_molecular_weight)

    if(real_egr):
        phi = get_equivalence_ratio(config,f,fb)
    
    # iCH4 = f.X[f.gas.species_index('CH4'),0]
    # iN2 = f.X[f.gas.species_index('N2'),0]
    # iCO2 = f.X[f.gas.species_index('CO2'),0]
    # iH2O = f.X[f.gas.species_index('H2O'),0]

    #for all non zero massfractions at inlet, build a list of species names and associated mass fractions
    i_species_names = ['Yi_'+specie for specie in f.gas.species_names if f.inlet.Y[f.gas.species_index(specie)] != 0.0]
    logprint('i_species_names',i_species_names)
    i_species_massfractions = [f.inlet.Y[f.gas.species_index(specie)] for specie in f.gas.species_names if f.inlet.Y[f.gas.species_index(specie)] != 0.0]
    logprint('i_species_massfractions',i_species_massfractions)

    vartosave = vars+i_species_names+species
    df =  pd.DataFrame(columns=vartosave)
    
    index = [f.gas.species_index(specie) for specie in species]
    df = pd.concat([df, pd.DataFrame([[egr,fb, phi,1/phi, P, T, T_flamme,SL0,flamme_thickness(f),f.density[0],f.density[-1],]+i_species_massfractions+list(f.X[index,-1])], columns=vartosave)]).astype(float) #+list(f.X[index][-1])] #

    return df

def fresh_gas(phi,config,egr_rate,Tmix,Pmix):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
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
    grid = np.linspace(0, width, 100)
    flame = ct.FreeFlame(mix, grid)
    return flame

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
    #print('Thickness (m)',thickness)
    return thickness

def compute_omega0(f):
    if(version.parse(ct.__version__) == version.parse("2.3.0")):
        species_names = f.gas.species_names
        net_prod_rate = f.net_production_rates
        omega0 = [(np.max(np.abs(net_prod_rate[i]))*f.gas[specie].molecular_weights).tolist()[0] for i,specie in enumerate(species_names)]

        logprint("-----------------------------------------")
        logprint("omega0 :")
        for i,specie in enumerate(species_names):
            logprint(f"{specie}: {omega0[i]:.6e} kg/m^3/s")
        
        return omega0
    if(version.parse(ct.__version__) == version.parse("3.0.0")):
        species_names = f.gas.species_names
        net_prod_rate = f.net_production_rates
        omega0 = [(np.max(np.abs(net_prod_rate[i]))*f.gas[specie].molecular_weights).tolist()[0] for i,specie in enumerate(species_names)]

        logprint("-----------------------------------------")
        logprint("omega0 :")
        for i,specie in enumerate(species_names):
            logprint(f"{specie}: {omega0[i]:.6e} kg/m^3/s")
        
        return omega0
    else:
        raise Exception('Cannot compute omega0 with cantera version '+ct.__version__+' (must be 2.3.0 or 3.0.0) or must be implemented for newer versions')
    
    
def alpha_beta(f):
    return f.X[f.gas.species_index('H2O'),0]/f.X[f.gas.species_index('O2'),0], f.X[f.gas.species_index('N2'),0]/f.X[f.gas.species_index('O2'),0]
    
def get_equivalence_ratio(config,f,blend_ratio):
    fuel_mol=0.0
    for fuel in config.compo.fuels:
        fuel_mol+=f.inlet.X[f.gas.species_index(fuel.split(':')[0])]
        #print('fuel_mol',fuel_mol)
    sx,_,_ = stoechiometric_ratios(config)

    Sx=sx[0]

    if(blend_ratio is not None and blend_ratio != 0):
        Sx=(1-blend_ratio)*sx[0]+(blend_ratio)*sx[1]

    phi = Sx*fuel_mol/f.inlet.X[f.gas.species_index('O2')]
    #print('ox',f.inlet.X[f.gas.species_index('O2')])
    logprint('Applied phi',phi)

    return phi

def apply_egr_to_inlet(f,config,phi,egr,fb,dry=False,T_reinj=None):
    config.gas.egr.X = f.X[:,-1]
    logprint("-----------------------------------------")
    logprint('EGR composition BEFORE drying operation',)
    species_names = config.gas.egr.species_names
    for i,specie in enumerate(species_names):
            logprint(f"{specie}: {config.gas.egr.X[i]:.6e} X[mol]")


    if(dry):
        # index = f.gas.species_index('H2O')
        # dry_egr = config.gas.egr.X
        # dry_egr[index] = 0.0
        # coef = 1/sum(dry_egr)
        # dry_egr = [coef*i for i in dry_egr]
        dry_egr = dryer(config,f,T_reinj)
        config.gas.egr.X = dry_egr

        logprint('EGR composition AFTER drying operation',)
        for i,specie in enumerate(species_names):
                logprint(f"{specie}: {config.gas.egr.X[i]:.6e} X[mol]")
    #keep only N2
    #config.gas.egr.X = [0.0 if specie != 'N2' else 1.0 for specie in species_names]
    
    # #keep N2 and CO2=1-X_N2
    # temp = [0.0 if specie != 'N2' else config.gas.egr.X[species_names.index('N2')] for specie in species_names]
    # temp = [1-temp[species_names.index('N2')] if (specie == 'CO2') else temp[species_names.index(specie)] for specie in species_names ]
    # config.gas.egr.X = temp
    # #keep N2 and CO2 and H2O=1-X_N2-X_CO2
    # temp = [0.0 if specie != 'N2' else config.gas.egr.X[species_names.index('N2')] for specie in species_names]
    # temp = [config.gas.egr.X[species_names.index('CO2')] if (specie == 'CO2') else temp[species_names.index(specie)] for specie in species_names ]
    # temp = [1-temp[species_names.index('N2')]-temp[species_names.index('CO2')] if (specie == 'H2O') else temp[species_names.index(specie)] for specie in species_names ]
    # config.gas.egr.X = temp
    
    # #keep N2 and CO2 and H2O and O2 or CH4
    # phi = get_equivalence_ratio(config,f,fb)
    # temp = [config.gas.egr.X[species_names.index(specie)] if specie in ['N2','CO2','H2O','H2'] else 0.0 for specie in species_names]
    # sum = np.sum(temp)
    # if phi<1:
    #     spec= 'O2'
    # else:
    #     spec= 'CH4'
    # temp = [1-sum if (specie == spec) else temp[species_names.index(specie)] for specie in species_names ]
   
    #temp = [1-sum if (specie == 'CH4') else temp[species_names.index(specie)] for specie in species_names ]
    # temp = [config.gas.egr.X[species_names.index('H2O')] if (specie == 'H2O') else temp[species_names.index(specie)] for specie in species_names ]
    # temp = [1-temp[species_names.index('N2')]-temp[species_names.index('CO2')-temp[species_names.index('H2O')]] if (specie == 'CH4') else temp[species_names.index(specie)] for specie in species_names ]
    # temp = [1-temp[species_names.index('N2')]-temp[species_names.index('CO2')-temp[species_names.index('H2O')]] if (specie == 'O2') else temp[species_names.index(specie)] for specie in species_names ]

    # config.gas.egr.X = temp


    print(config.gas.egr.X)

    logprint("-----------------------------------------")
    T,P,X = mixer(config,phi,egr,fb,real_egr=True,T_reinj=T_reinj)
    f.inlet.X = X
    f.inlet.T = T

    return f



def flame_iteration(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,i):
    XCO2_1 = 0
    if(verbose>0):
        if (verbose>1):
            logprint('Iteration n'+str(i),file=sys.stdout)
        else:
            logprint('Iteration n'+str(i),file=sys.stdout)
        
        try:
            if(verbose>1):
                logprint('Inlet composition before f.solve',f.inlet.X,file=sys.stdout)
            f.solve(loglevel, refine_grid, auto)
            success = True
            if (real_egr):
                index = f.gas.species_index('H2')
                XCO2_1 = f.X[index,-1]
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    logprint('EGR applied to inlet',file=sys.stdout)
            
        except FlameExtinguished:
            logprint('Flame '+flamename +': ',file=sys.stdout)
            logprint('Flame extinguished',file=sys.stdout)
            success = False

        except ct.CanteraError as e:
            logprint('Flame '+flamename +': ',file=sys.stdout)
            logprint(('Error occurred while solving: (ite '+str(i)+') \n', e),file=sys.stdout)
            success = False

        return f,success,XCO2_1
    
def boucle_over_flame_iterations(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,list_of_criteria):


    for i,criteria in enumerate(list_of_criteria):
        f.energy_enabled = criteria[4]
        f.set_refine_criteria(ratio = criteria[0], slope = criteria[1], curve = criteria[2], prune = criteria[3])
        f,success,XCO2 = flame_iteration(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,i) 
        if (success):
            pass
        else:
            break

    return f,success,XCO2




