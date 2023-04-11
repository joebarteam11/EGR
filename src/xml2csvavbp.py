import cantera as ct
import os

path=os.getcwd()+"/src/data/"

f = ct.FreeFlame(ct.Solution('schemes/Aramco13.cti'), width=0.03)
file=path+'egr0.3_phi0.85_T300.0_P1.0.xml'
f.restore(file, loglevel=0)
f.write_AVBP(os.getcwd()+'/results/'+'Sol-CANTERA_P-100000-T-300.0-Phi-0.85.csv')
print('Done')
