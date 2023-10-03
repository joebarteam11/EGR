from lib_egr_260 import *

def solve_flame(f,hash,flamename,config,phi,egr,fb,real_egr=False,dry=False,T_reinj=None):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...   
    
    #################################################################
    # Iterations start here
    #################################################################
    maxegrate_iter = 50
    residuals = []
    index = f.gas.species_index('CO2')
    last_residual = 1.0 
    i=0
    verbose = 1
    loglevel  = 0                      # amount of diagnostic output (0 to 5)	    
    refine_grid = 'refine' #True                  # True to enable refinement, False to disable 	
    f.max_time_step_count=50000
    f.max_grid_points=1000
    auto = False

    # first iteration
    ##################
    # No energy equilibrium activated (reduce the calculation)
    f.energy_enabled = False
    # Mesh refinement
    f.set_refine_criteria(ratio = 7.0, slope = 0.95, curve = 0.95)
    # Max number of times the Jacobian will be used before it must be re-evaluated
    f.set_max_jac_age(50, 50)
    #Set time steps whenever Newton convergence fails
    f.set_time_step(1e-08, [25, 40, 80, 140, 200, 350, 500, 700, 1000, 1300, 1700, 2000, 3000, 5000, 10000, 12000, 15000, 20000]) #s

    
    while(True or i<maxegrate_iter):
        #print('While loop starts')
        # Calculation
        if(verbose>0):
            if(verbose>1):
                logprint('1st iteration...',file=sys.stdout)
            else:
                logprint('1st iteration...')
            logprint('on flame '+ flamename)
        try:
            if(verbose>1):
                logprint('Inlet composition before f.solve',f.inlet.X)
            f.solve(loglevel, refine_grid, auto=True)
            
            if(real_egr):
                XCO2_1 = f.X[index,-1]
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    logprint('EGR applied to inlet')
                
            #f.save(flamefile)
        except FlameExtinguished:
            logprint('Flame '+flamename +': ')
            logprint('Flame extinguished')
            
        except ct.CanteraError as e:
            logprint('Flame '+flamename +': ')
            logprint(('Error occurred while solving: (ite 1) \n', e))
            
        # Second iteration
        #################
        # Energy equation activated
        f.energy_enabled = True
        # mesh refinement
        f.set_refine_criteria(ratio = 7.5, slope = 0.85, curve = 0.85)
        # Calculation
        if(verbose>0):
            if(verbose>1):
                logprint('2nd iteration...',file=sys.stdout)
            else:
                logprint('2nd iteration...')
            logprint('on flame '+ flamename)
        try:
            if(verbose>1):
                logprint('Inlet composition before f.solve',f.inlet.X)
            f.solve(loglevel, refine_grid, auto=auto)
            
            if(real_egr):
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    logprint('EGR applied to inlet')
            #f.save(flamefile)
        except FlameExtinguished:
            logprint('Flame '+flamename +': ')
            logprint('Flame extinguished')
            
        except ct.CanteraError as e:
            logprint('Flame '+flamename +': ')
            logprint(('Error occurred while solving: (ite 2)', e))
            
        # Third iteration
        #################
        # On raffine le maillage
        f.set_refine_criteria(ratio = 5.0, slope = 0.4, curve = 0.4, prune = 0.03)
        # Calculation
        if(verbose>0):
            if(verbose>1):
                logprint('3rd iteration...',file=sys.stdout)
            else:
                logprint('3rd iteration...')
            logprint('on flame '+ flamename)
        try:
            if(verbose>1):
                logprint('Inlet composition before f.solve',f.inlet.X)
            f.solve(loglevel, refine_grid, auto=auto)
            if(real_egr):
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    logprint('EGR applied to inlet')
            #f.save(flamefile)
        except FlameExtinguished:
            logprint('Flame '+flamename +': ')
            logprint('Flame extinguished')
            
        except ct.CanteraError as e:
            logprint('Flame '+flamename +': ')
            logprint(('Error occurred while solving: (ite 3)', e))
            
        #################
        # Fourth iteration
        # Mesh refinement
        f.set_refine_criteria(ratio = 5.0, slope = 0.1, curve = 0.1, prune = 0.01)
        # Calculation
        if(verbose>0):
            if(verbose>1):
                logprint('4th iteration...',file=sys.stdout)
            else:  
                logprint('4th iteration...')
            logprint('on flame '+ flamename)
        try:
            if(verbose>1):
                logprint('Inlet composition before f.solve',f.inlet.X)
            f.solve(loglevel, refine_grid, auto=auto)
            if(real_egr):
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    logprint('EGR applied to inlet')
            #f.save(flamefile)
        except FlameExtinguished:
            logprint('Flame '+flamename +': ')
            logprint('Flame extinguished')
            
        except ct.CanteraError as e:
            logprint('Flame '+flamename +': ')
            logprint(('Error occurred while solving: (ite 4)', e))
            
        # Fifth iteration
        #################
        #Mesh refinement
        f.set_refine_criteria(ratio = 5.0, slope = 0.05, curve = 0.05, prune = 0.01)
        # Calculation
        if(verbose>0):
            if(verbose>1):
                logprint('5th iteration...',file=sys.stdout)
            else:
                logprint('5th iteration...')
            logprint('on flame '+ flamename)
        try:
            if(verbose>1):
                logprint('Inlet composition before f.solve',f.inlet.X)
            f.solve(loglevel, refine_grid)
            if(real_egr):
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    logprint('EGR applied to inlet')
            # try:
            #     f.write_hdf(flamefile)
            # except:
            #     f.save(flamefile[:-3]+'.yaml')
        except FlameExtinguished:
            logprint('Flame '+flamename +': ')
            logprint('Flame extinguished')
            
        except ct.CanteraError as e:
            logprint('Flame '+flamename +': ')
            logprint(('Error occurred while solving: (ite 5)', e))
        # # Sixth iteration
        #################
        #NO Mesh refinement
        f.set_refine_criteria(ratio = 5.0, slope = 0.05, curve = 0.05, prune = 0.01)
        # Calculation
        if(verbose>0):
            if(verbose>1):
                logprint('6th iteration...',file=sys.stdout)
            else:
                logprint('6th iteration...')
            logprint('on flame '+ flamename)
            
        try:
            if(verbose>1):
                logprint('Inlet composition before f.solve',f.inlet.X)
            f.solve(loglevel, refine_grid='disabled')
            if(real_egr):
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    logprint('EGR applied to inlet')
            #f.save(flamefile)
        except FlameExtinguished:
            logprint('Flame '+flamename +': ')
            logprint('Flame extinguished')
            
        except ct.CanteraError as e:
            logprint('Flame '+flamename +': ')
            logprint(('Error occurred while solving: (ite 6)', e))
        
        # # Seventh iteration
        #################
        #NO Mesh refinement
        f.set_refine_criteria(ratio = 5.0, slope = 0.05, curve = 0.05, prune = 0.01)
        # Calculation
        if(verbose>0):
            if(verbose>1):
                logprint('7th iteration...',file=sys.stdout)
            else:
                logprint('7th iteration...')
            logprint('on flame '+ flamename)
        try:
            if(verbose>1):
                logprint('Inlet composition before f.solve',f.inlet.X)
            f.solve(loglevel, refine_grid='disabled')
            T_flamme = f.T[-1]
            if(real_egr):
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    logprint('EGR applied to inlet')
            #f.save(flamefile)
            try:
                # if hash is an existing file 
                if os.path.isfile(hash):
                    logprint('Flame file '+hash+' is already existant, will not overwrite')
                else:
                    logprint('Saving flame to file '+hash)
                    f.save(hash, loglevel=loglevel, description=flamename + "\n Cantera version "+ct.__version__+' with '+config.transport+' transport model and mechanism '+config.scheme )
            except:
                logprint('Cannot save flame to file '+flamename)

        except FlameExtinguished:
            T_flamme= np.NaN
            logprint('Flame '+flamename +': ')
            logprint('Flame extinguished')
            
            break

        except ct.CanteraError as e:
            T_flamme = np.Inf
            logprint('Flame '+flamename +': ')
            logprint(('Error occurred while solving: (ite 7)', e))
            break

        if(real_egr):
            last_residual = abs(residuals[-1]-residuals[-2])
            #print('last residual',last_residual)
            #print('Residuals',residuals)
        else:
            break
        #print('last residual',last_residual)
        if(last_residual<1e-9):
            T_flamme = f.T[-1]
            #print('BREAK HERE')
            break
        i+=1
        if(i>maxegrate_iter):
            T_flamme = np.Inf
            logprint('Too many iterations of solve_flame for real_egr loop, for flame ',hash,' cannot converge')
            break
    return f,T_flamme

def compute_solutions_1D(config,phi,tin,pin,egr,fb,restart_rate,real_egr,dry=True,T_reinj=None,species = ['CH4','O2','CO2','H2O'],vars=['EGR','FB','phi','AF','P','Tin','T','u','dF','rhoGF','rhoGB']):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
    path = os.getcwd()+'/src'
    vartosave = vars+species
    df =  pd.DataFrame(columns=vartosave)

    tol_ss = [2.0e-9, 1.0e-9]  # tolerance [rtol atol] for steady-state problem
    tol_ts = [2.0e-9, 1.0e-9]  # tolerance [rtol atol] for time stepping

    #create the gas object containing the mixture of fuel, ox and egr and all thermo data
    _, config.gas.fuels = create_reservoir(config,config.compo.fuels,tin[0], pin[0],blend_ratio=fb)
    _, config.gas.ox = create_reservoir(config,config.compo.ox, tin[1], pin[1],scheme='air.xml')
    _, config.gas.egr = create_reservoir(config,config.compo.egr, tin[2], pin[2])

    #get the temperature and pressure of the mixture according to phi and egr rate
    T,P,X = mixer(config,phi,egr,fb)
    gas=fresh_gas(phi,config,egr,T,P)
    f = build_freeflame(gas)

    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
    f.transport_model = config.transport

    f.inlet.T = T
    f.P = P
    f.inlet.X = X

    flamename,hash = generate_unique_filename(config,f)
    logprint('Flame hash: '+hash)

    flamefile = path+'/data/'+hash+'.yaml'
    if os.path.isfile(flamefile):
        logprint('Flame file '+flamefile+' is found')
        try:
            f.restore(flamefile, loglevel=1)
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

    #print('Inlet composition (mol)',f.inlet.X)
    #print('Inlet composition (mass)',f.inlet.Y)

    f,T_flamme = solve_flame(f,flamefile,flamename,config,phi,egr,fb,real_egr=real_egr,dry=dry,T_reinj=T_reinj)

    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        SL0=f.velocity[0]
    else:
        SL0=f.u[0]

    omega0=compute_omega0(f)
    #f.write_AVBP(path+'/'+'CH4_phi073'+'.csv')
    #print("Mean molecular weight: ",f.gas.mean_molecular_weight)

    if(real_egr):
        phi = get_equivalence_ratio(config,f,fb)

    index = [f.gas.species_index(specie) for specie in species]
    df = pd.concat([df, pd.DataFrame([[egr,fb, phi,1/phi, P, T, T_flamme,SL0,flamme_thickness(f),f.density[0],f.density[-1]]+list(f.X[index,-1])], columns=vartosave)]).astype(float) #+list(f.X[index][-1])] #

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
    flame = ct.FreeFlame(mix, width=width)
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
    else:
        raise Exception('Cannot compute omega0 with cantera version '+ct.__version__+' (must be 2.3.0) or must be implemented for newer versions')
    
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
    #print('phi',phi)
    
    return phi