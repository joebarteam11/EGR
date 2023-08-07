from lib_egr_260 import *

def solve_flame(f,flametitle,config,phi,egr,fb,real_egr=False,dry=False,T_reinj=None):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...   
    
    #################################################################
    # Iterations start here
    #################################################################
    maxegrate_iter = 50
    residuals = []
    index = f.gas.species_index('CO2')
    last_residual = 1.0 
    i=0
    verbose = 0
    loglevel  = 0                      # amount of diagnostic output (0 to 5)	    
    refine_grid = 'refine' #True                  # True to enable refinement, False to disable 	
    f.max_time_step_count=50000
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
    f.set_time_step(1e-07, [25, 40, 80, 140, 200, 350, 500, 700, 1000, 1300, 1700, 2000, 3000, 5000, 10000, 12000, 15000, 20000]) #s

    while(True or i<maxegrate_iter):
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
        f.set_refine_criteria(ratio = 5.0, slope = 0.05, curve = 0.05, prune = 0.01)
        # Calculation
        if(verbose>0):
            print('6th iteration...')
        try:
            if(verbose>1):
                print('Inlet composition before f.solve',f.inlet.X)
            f.solve(loglevel, refine_grid='disabled')
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
        f.set_refine_criteria(ratio = 5.0, slope = 0.05, curve = 0.05, prune = 0.01)
        # Calculation
        if(verbose>0):
            print('7th iteration...')
        try:
            if(verbose>1):
                print('Inlet composition before f.solve',f.inlet.X)
            f.solve(loglevel, refine_grid='disabled')
            if(real_egr):
                XCO2_2 = f.X[index,-1]
                residuals.append(np.abs(XCO2_1-XCO2_2))
                XCO2_1 = XCO2_2
                f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
                if(verbose>0):
                    print('EGR applied to inlet')
            #f.save(flametitle)
            try:
                f.save(flametitle, loglevel=loglevel)
            except:
                logprint('Cannot save flame to file '+flametitle)

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
        if(last_residual<1e-9):
            #print('BREAK HERE')
            break
        i+=1
        if(i>maxegrate_iter):
            raise Exception('Too many iterations of solve_flame, cannot converge')
    return f

def compute_solutions_1D(config,phi,tin,pin,egr,fb,restart_rate,real_egr,dry=True,T_reinj=None,species = ['CH4','O2','CO2','H2O'],vars=['EGR','FB','phi','AF','P','Tin','T','u','dF']):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
    path = os.getcwd()+'/src'
    vartosave = vars+species
    df =  pd.DataFrame(columns=vartosave)

    tol_ss = [2.0e-7, 1.0e-9]  # tolerance [rtol atol] for steady-state problem
    tol_ts = [2.0e-7, 1.0e-9]  # tolerance [rtol atol] for time stepping

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

    if(real_egr):
        f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
    else:
        f.inlet.T = T
        f.P = P
        f.inlet.X = X

    flamename = generate_unique_filename(f)
    logprint('Flame hash: '+flamename)

    flametitle = path+'/data/'+flamename+'.xml'
    if os.path.isfile(flametitle):
        f.restore(flametitle, loglevel=0)
        logprint('Flame restored from file '+flametitle)
    else:
        f.set_initial_guess()
        logprint('Flame initial guess set')

    #print('Inlet composition (mol)',f.inlet.X)
    #print('Inlet composition (mass)',f.inlet.Y)

    solve_flame(f,flametitle,config,phi,egr,fb,real_egr=real_egr,dry=dry,T_reinj=T_reinj)

    if(version.parse(ct.__version__) >= version.parse('2.5.0')):
        SL0=f.velocity[0]
    else:
        SL0=f.u[0]

    omega0=compute_omega0(f)

    if(real_egr):
        phi = get_equivalence_ratio(config,f,fb)

    index = [f.gas.species_index(specie) for specie in species]
    df = pd.concat([df, pd.DataFrame([[egr,fb, phi,1/phi, P, T, f.T[-1],SL0,flamme_thickness(f)]+list(f.X[index,-1])], columns=vartosave)]).astype(float) #+list(f.X[index][-1])] #

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