#!/usr/bin/python
#-*- coding: latin-1 -*-
import os, sys
"""
Created on Mon Nov 14 14:11:59 2016

project 4: Exhaust gas recirculation (EGR) combustion technology in gas turbines
Combustion 2 course
Department of Energy and Propulsion
INSA Rouen-Normandie
2018-2019 term


Author: Pradip Xavier & Bruno Renou
"""

###############################################################################
#import 
###############################################################################
import cantera as ct
import numpy as np
import csv 
import matplotlib.pylab as plt 

###############################################################################
#program starts here
###############################################################################
# object creation with a specific chemical scheme
gas1 = ct.Solution('gri30.cti','gri30_mix')

# Initial conditions
Tin=293      # Temperature (K)
Pin=1*101325      # Pressure (Pa)
phi=1             # equivalence ratio

ifuel = gas1.species_index('CH4')
io2 = gas1.species_index('O2')
in2 = gas1.species_index('N2')
ico2 = gas1.species_index('CO2')
io = gas1.species_index('O')
ih = gas1.species_index('H')
ioh= gas1.species_index('OH')


air_N2_O2_ratio = 3.76
stoich_O2 = 2
#b=0.625
b=0

X = np.zeros(gas1.n_species)     # Molar fractions array of all species (size 'n_species')
X[ifuel] = phi
X[io2] = stoich_O2
X[in2] = stoich_O2*air_N2_O2_ratio
X[ico2] = stoich_O2*b
gas1.TPX = Tin, Pin, X           # Gas initialization with T, P and X input


# Initialization of the 1D flame
# Grid specfications with refinement at inlet and outlet, 6 points in x-direction :
initial_grid = 5*np.array([0.0, 0.001, 0.01, 0.02, 0.029, 0.03],'f')/3 # m


#Set tolerance properties
tol_ss    = [1.0e-5, 1.0e-8]        # [rtol atol] for steady-state problem
tol_ts    = [1.0e-5, 1.0e-8]        # [rtol atol] for time stepping

loglevel  = 1                       # amount of diagnostic output (0 to 5)	    
refine_grid = True                  # True to enable refinement, False to disable 				   


# Creation of the flame object
f = ct.FreeFlame(gas1, initial_grid)

f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

# Boundary conditions
f.inlet.X = X
f.inlet.T = Tin

#################################################################
# Iterations start here
#################################################################

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

##################################################################
# Flame speed

# Kinematic approach
print ('Kinematic approach: ', 'Flame speed = ',f.velocity,'m/s')

# Kinetic approach

M=16 #masse molaire du fuel
wk=M*np.trapz(f.net_production_rates[ifuel,:], f.flame.grid) #int�gral de tous les wk pour avoir le grand wk (en kg/m3/s)
print ('Taux de production wk= ', wk, ' kg/m3/s')

Y1= f.Y[ifuel,0]
Y2= f.Y[ifuel, len(f.flame.grid)-1]
print ('Y1= ', Y1, '  Y2= ', Y2)

rho=f.density_mass[0]
Sl=(-1.0/(rho*(Y1-Y2)))*wk  #0.665 = masse volumique du m�lange
print ('Kinetic approach: Flame speed = ', Sl, 'm/s')


#################################################################
# Result saving
#################################################################
f.write_csv('ch4_air_adiabatic_phi'+str(phi)+'-T'+str(Tin)+'-P'+str(Pin)+'.csv', quiet=False)


#################################################################
# Results plot
#################################################################
#Plot the velocity, temperature, density
z = f.flame.grid
T = f.T
u = f.velocity


fig, ax1 = plt.subplots()

ax1.plot(z,f.X[ifuel],'b--',label='$X_{CH4}$')
ax1.plot(z,f.X[io2],'b-',label='$X_{O2}$')
ax1.set_xlabel('axial position (m)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Mole fractions', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')
ax2 = ax1.twinx()
ax2.plot(z,T,'r-',label='Temperature')
ax2.set_ylabel('Temperature (K)', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r') 
    
legend = ax1.legend(loc='center right', shadow=True)

#plt.savefig('FPF_CH4_phi1_mole_fraction.eps', format='eps', dpi=1000)
#plt.savefig('FPF_CH4_phi1_mole_fraction.png', bbox_inches='tight')

#plt.xlim((0.0368,0.038))


fig, ax1 = plt.subplots()

ax1.plot(z,u,'b--')
ax1.set_xlabel('axial position (m)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Flow Velocity (m/s)', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')
ax2 = ax1.twinx()
ax2.plot(z,T,'r-',label='Temperature')
ax2.set_ylabel('Temperature (K)', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r') 
#plt.savefig('FPF_CH4_phi1_velocity.eps', format='eps', dpi=1000)
#plt.savefig('FPF_CH4_phi1_velocity.png', bbox_inches='tight')

#plt.xlim((0.0368,0.038))


fig, ax1 = plt.subplots()

ax1.plot(z,f.density_mass,'b--')

ax1.set_xlabel('axial position (m)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Flow density (kg/m3)', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')
ax2 = ax1.twinx()
ax2.plot(z,T,'r-',label='Temperature')
ax2.set_ylabel('Temperature (K)', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r') 
#plt.savefig('FPF_CH4_phi1_density.eps', format='eps', dpi=1000)
#plt.savefig('FPF_CH4_phi1_density.png', bbox_inches='tight')
#    legend = ax1.legend(loc='upper left', shadow=True)
    
#plt.xlim(0.0368,0.038)
plt.show()


############################################################################    
# plot avec CO2
print ('plot avec CO2')  

fig, ax1 = plt.subplots()

ax1.plot(z,f.Y[io],'b--',label='$Y_{O}$')
ax1.plot(z,f.X[ih],'c-',label='$Y_{H}$')
ax1.plot(z,f.X[ioh],'k-',label='$Y_{OH}$')
ax1.set_xlabel('axial position (m)')
# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Mass fractions', color='b')
for tl in ax1.get_yticklabels():
    tl.set_color('b')
ax2 = ax1.twinx()
ax2.plot(z,T,'r-',label='Temperature')
ax2.set_ylabel('Temperature (K)', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r') 
    
legend = ax1.legend(loc='center right', shadow=True)
 

plt.show()
