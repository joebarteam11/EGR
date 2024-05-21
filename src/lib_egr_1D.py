from lib_egr_260 import *
import h5py
from write_table import write_table_hdf
import os

try_restore = False


def flame_saver(f,loglevel,flamename,hash,config):
    try:
        # if hash is an existing file 
        # if os.path.isfile(hash):
        #     # logprint('Flame file '+hash+' is already existant, will not overwrite')
        #     f.save(hash, loglevel=loglevel, description=flamename + "\n Cantera version "+ct.__version__+' with '+config.transport+' transport model and mechanism '+config.scheme )
        # else:
        logprint('Saving flame to file '+hash)
        f.save(hash, loglevel=loglevel, description=flamename + "\n Cantera version "+ct.__version__+' with '+config.transport+' transport model and mechanism '+config.scheme )
    except:
        logprint('Cannot save flame to file '+flamename)

def flame_saver_csv(f,config,phi,egr,fb):
    if config.saveCSV:
        if not os.path.isdir(config.saveCSVpath):
            os.makedirs(config.saveCSVpath)
        
        flamename = config.saveCSVpath+'/T'+str(config.gas.fuels.T)+'K_P'+str(round(config.gas.fuels.P/1e5,2))+'bar_phi'+str(round(phi, 2))+'_egr'+str(egr)+'_fb'+str(fb)
        f.write_csv(flamename+'.csv',quiet=False)
        logprint('Flame csv saved : '+flamename)



        
def solve_flame(f,hash,flamename,config,phi,egr,fb,real_egr=False,dry=False,T_reinj=None,restore_success=False):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...   
    
    #################################################################
    # Iterations start here
    #################################################################
    maxegrate_iter = 50
    residuals = [2.0]
    index = f.gas.species_index('CO2')
    last_residual = 1.0 
    XCO2_1 = 0.0
    i=0
    verbose = 0
    loglevel  = 0                      # amount of diagnostic output (0 to 5)	    
    refine_grid = 'refine' #True                  # True to enable refinement, False to disable 	
    f.max_time_step_count=10000
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
    f.set_time_step(1e-08, [25, 40, 80, 140, 200]) #s

    # criteria_list determines the number of iterations that you want and tells : ratio, slope, curve, prune and energy enabled
    criteria_list = [[7.0, 0.95, 0.95,0, False]] # First iteration criterias ( 0 prune means disable)
    criteria_list += [[5.0, 0.5, 0.5,0, True]] # Second iteration criterias
    criteria_list += [[5.0, 0.4, 0.4,0.03, True]] # Third iteration criterias
    criteria_list += [[5.0, 0.1, 0.1,0.01, True]] # Fourth iteration criterias
    criteria_list += [[5.0, 0.05, 0.05,0.01, True]] # Fifth iteration criterias
    criteria_list += [[4.0, 0.05, 0.05,0.01, True]] # Sixth iteration criterias
    last_iteration = [3.0, 0.05, 0.05, 0.01, True] # 7th iteration criterias


    while(True or i<maxegrate_iter):
        # Calculation
        if(not restore_success or real_egr):
            auto_success = False
            auto = True 
            # One loop of auto refinement
            itype = " auto iteration"
            if (not real_egr):
                f,auto_success,_,_ = flame_iteration(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,itype,XCO2_1)
            #residuals.append(residual)
            if(not auto_success or real_egr):
                auto = False
                f,sucess,XCO2_1,residu = boucle_over_flame_iterations(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,criteria_list)
                residuals.extend(residu)
            else:
                logprint('Successful auto refinement, going to last_iteration')
        else:
            logprint('Successful restore, going to last_iteration')
        #end if not restore_success

        # # Last iteration
        f.energy_enabled = True
        f.set_refine_criteria(ratio = last_iteration[0], slope = last_iteration[1], curve = last_iteration[2], prune = last_iteration[3])
        T_flamme = f.T[-1]
        itype = " last iteration"
        f,sucess,XCO2_1,residual = flame_iteration(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,itype,XCO2_1)
        residuals.append(residual)

        if sucess:
            flame_saver(f,loglevel,flamename,hash,config)
            flame_saver_csv(f,config,phi,egr,fb)
        else:
            pass


        if(real_egr):
            last_residual = abs(residuals[-1]-residuals[-2])
        else:
            break
        if(last_residual<1e-6):
            T_flamme = f.T[-1]
            break
        i+=1
        if(i>maxegrate_iter):
            T_flamme = f.T[-1]
            logprint('Too many iterations of solve_flame for real_egr loop, for flame ',hash,' cannot converge')
            break
    return f,T_flamme

def compute_solutions_1D(config,phi,tin,pin,egr,fb,restart_rate,real_egr,tol_ss=[2.0e-5, 1.0e-9],tol_ts=[2.0e-5, 1.0e-9],dry=True,T_reinj=None,species = ['CH4','O2','CO2','H2O'],vars=['EGR','FB','phi','AF','P','Tin','T','u','dF','hr','rhoGF','rhoGB','omega0']):
    warnings.simplefilter("ignore", UserWarning) #aramco speeks a lot...
    path = os.getcwd()+'/src'


    # tol_ss = [2.0e-5, 1.0e-9]  # tolerance [rtol atol] for steady-state problem
    # tol_ts = [2.0e-5, 1.0e-9]  # tolerance [rtol atol] for time stepping

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
    restore_success = False
    if os.path.isfile(flamefile) and try_restore:
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

    # max of omega0
    omega0=compute_omega0_max(f)
    #max of heat release rate
    hr_max = np.max(f.heat_release_rate)

    #f.write_AVBP(path+'/'+'CH4_AP_phi08'+'.csv')
    #print("Mean molecular weight: ",f.gas.mean_molecular_weight)

    if(real_egr):
        phi = get_equivalence_ratio(config,f,fb)
    
    # iCH4 = f.X[f.gas.species_index('CH4'),0]
    # iN2 = f.X[f.gas.species_index('N2'),0]
    # iCO2 = f.X[f.gas.species_index('CO2'),0]
    # iH2O = f.X[f.gas.species_index('H2O'),0]

    #for all non zero massfractions at inlet, build a list of species names and associated mass fractions
    # i_species_names = ['Yi_'+specie for specie in f.gas.species_names if f.inlet.Y[f.gas.species_index(specie)] != 0.0]
    # logprint('i_species_names',i_species_names)
    # i_species_massfractions = [f.inlet.Y[f.gas.species_index(specie)] for specie in f.gas.species_names if f.inlet.Y[f.gas.species_index(specie)] != 0.0]
    # logprint('i_species_massfractions',i_species_massfractions)

    vartosave = vars+species
                # i_species_names+
                # species
    df =  pd.DataFrame(columns=vartosave)

    #ambiguous definition with fuel blends
    max_omega0_fuel = np.max([omega0[f.gas.species_index(specie.split(':')[0])] for specie in config.compo.fuels])
    
    index = [f.gas.species_index(specie) for specie in species]
    df = pd.concat([df, pd.DataFrame([[egr,fb, phi,1/phi, P, T, T_flamme,SL0,flamme_thickness(f),hr_max,f.density[0],f.density[-1],max_omega0_fuel,]+list(f.X[index,-1])], columns=vartosave)]).astype(float) #+list(f.X[index][-1])] #
    # duplicate 10 times the first line
    #df = pd.concat([df]*5, ignore_index=True)
    #print(df)

    #table_generation(config,df,path)
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

def compute_omega0_max(f):
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


    #print(config.gas.egr.X)

    logprint("-----------------------------------------")
    T,P,X = mixer(config,phi,egr,fb,real_egr=True,T_reinj=T_reinj)
    f.inlet.X = X
    f.inlet.T = T

    return f

def table_generation(config,df,output_folder):
    
    data_names = {

        'u':'SL_0',
        'Tin':'TEMPERATURE',
        'T':'TEMPERATURE_BURNT',
        'P':'PRESSURE',
        'phi':'EQUIVALENCE_RATIO',
        'hr':'VOLUMETRIC_HEAT_RELEASE_MAX',
        'omega0':'W_FUEL_MAX',
        # 'omegaSpec':'W_YC_MAX',
        'dF':'DELTA_L_0',
    }
    # print("-THIS IS DF -"*400)
    # print(df)
    # print("--"*400)
    # sys.stdout.flush()
    # convert df columns names to data_names and only keep the columns that are in data_names
    tempdf = df.rename(columns=data_names)[data_names.values()]
    # remove ['EQUIVALENCE_RATIO','TEMPERATURE','PRESSURE',] from tempdf
    # tempdf = tempdf.drop(columns=['EQUIVALENCE_RATIO','TEMPERATURE','PRESSURE',])
    # print("-THIS IS tempdf -"*400)
    # print(tempdf)
    # print("--"*400)
    # sys.stdout.flush()

    #create a table object and fill it with the data
    table = Table(config,tempdf,output_folder)
    # print(table.data[table.variable_names[0]])
    # print(table.data['SL_0'])
    #write the table to hdf5 file
    write_table_hdf(table)
    # logprint('Table written to hdf format in provided output folder',sys.stdout)
    
#A class to store the data to be written in the hdf5 file
class Table:
    def __init__(self, config, dfdata, output_folder, phivitiated=False, tool='TFLES'):
        self.output_folder = output_folder
        self.variable_names = list(dfdata.columns.values)
        self.list_sorted_phi = sorted(set(dfdata['EQUIVALENCE_RATIO']))
        self.iPhi = [self.list_sorted_phi.index(phi) for phi in dfdata['EQUIVALENCE_RATIO']]
        self.list_sorted_temperature = sorted(set(dfdata['TEMPERATURE']))
        self.iTemperature = [self.list_sorted_temperature.index(T) for T in dfdata['TEMPERATURE']]
        self.list_sorted_pressure = sorted(set(dfdata['PRESSURE']))
        self.iPressure = [self.list_sorted_pressure.index(P) for P in dfdata['PRESSURE']]
        #remove equiv ratio, temp and pressure from data from self.variable_names
        self.variable_names = [var for var in self.variable_names if var not in ['EQUIVALENCE_RATIO','TEMPERATURE','PRESSURE']]  
        
        # print("-THIS IS lsitsorted -"*400)
        # print(self.list_sorted_phi)
        # print(self.list_sorted_temperature)
        # print(self.list_sorted_pressure)
        # print("--"*400)


        # print("-THIS IS dfdata -"*400)
        # print(dfdata)     
        # print("--"*400)
        # sys.stdout.flush()
        #data = data.drop(columns=['EQUIVALENCE_RATIO','TEMPERATURE','PRESSURE',])
        #define self.data is multiindex

        #convert dfdata to a 3D array of shape (len(self.list_sorted_phi),len(self.list_sorted_temperature),len(self.list_sorted_pressure)) and fill it with the data
        #conserve only the columns that are in self.variable_names and put all these arrays in self.data
        self.data = {var:np.zeros((len(self.list_sorted_phi),len(self.list_sorted_temperature),len(self.list_sorted_pressure))) for var in self.variable_names}

        for var in self.variable_names:
            for i,phi in enumerate(self.list_sorted_phi):
                for j,T in enumerate(self.list_sorted_temperature):
                    for k,P in enumerate(self.list_sorted_pressure):
                        # print("-THIS IS seld.data -"*400)
                        # print(self.data)
                        # print("--"*400)
                        # sys.stdout.flush()
                        self.data[var][i,j,k] = dfdata[(dfdata['EQUIVALENCE_RATIO']==phi) & (dfdata['TEMPERATURE']==T) & (dfdata['PRESSURE']==P)][var].values[0]
        
        self.mechanism = config.scheme.split('/')[-1].split('.')[0]
        self.phivitiated = phivitiated
        self.tool = tool


def flame_iteration(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,i,XCO2_1):

    if(verbose>0):
        if (verbose>1):
            logprint('Iteration n'+str(i),file=sys.stdout)
        else:
            logprint('Iteration n'+str(i))
        
    try:
        residual = 1.0
        if(verbose>1):
            logprint('Inlet composition before f.solve',f.inlet.X,file=sys.stdout)
        f.solve(loglevel, refine_grid, auto)
        success = True
        # print('before egr')
        if (real_egr and XCO2_1==0.0):
            # print('init egr loop')
            index = f.gas.species_index('CO2')
            XCO2_1 = f.X[index,-1]
            f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
            if(verbose>0):
                logprint('EGR applied to inlet',file=sys.stdout)
        elif(real_egr):
            # print('egr loop')
            index = f.gas.species_index('CO2')
            XCO2_2 = f.X[index,-1]
            residual = np.abs(XCO2_1-XCO2_2)
            #residuals.append(np.abs(XCO2_1-XCO2_2))
            XCO2_1 = XCO2_2
            f = apply_egr_to_inlet(f,config,phi,egr,fb,dry,T_reinj)
            if(verbose>0):
                logprint('EGR applied to inlet')
        else:
            XCO2_1 = 0.0

    except FlameExtinguished:
        logprint('Flame '+flamename +': ',file=sys.stdout)
        logprint('Flame extinguished',file=sys.stdout)
        success = False

    except ct.CanteraError as e:
        logprint('Flame '+flamename +': ',file=sys.stdout)
        logprint(('Error occurred while solving: (ite '+str(i)+') \n', e),file=sys.stdout)
        success = False

    # print('success',success)
    # print('XCO2_1',XCO2_1)
    # print('residual',residual)

    return f,success,XCO2_1,residual
    
def boucle_over_flame_iterations(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,list_of_criteria):

    XCO2 = 0.0
    residuals = []
    for i,criteria in enumerate(list_of_criteria):
        f.energy_enabled = criteria[4]
        f.set_refine_criteria(ratio = criteria[0], slope = criteria[1], curve = criteria[2], prune = criteria[3])
        f,success,XCO2,residual = flame_iteration(f,verbose,loglevel,refine_grid,auto,real_egr,flamename,config,phi,egr,fb,dry,T_reinj,i,XCO2)
        residuals.append(residual)
        if (success):
            pass
        else:
            break

    return f,success,XCO2,residuals




