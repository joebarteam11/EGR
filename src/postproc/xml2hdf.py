import cantera as ct
import os
from packaging import version
import pandas as pd 
from lib_egr_260 import flamme_thickness

#create h5 restore files from xml and save them in the subfolder 'data'
#also create a csv for plotting script


path=os.getcwd()+"/src/data/"
vars=['EGR','phi','P','Tin','T','u','dF']
species = ['CH4','H2','O2','CO','CO2','H2O']
vartosave = vars+species
df =  pd.DataFrame(columns=vartosave)

f = ct.FreeFlame(ct.Solution('schemes/Aramco13.cti'), width=0.03)

#create a for loop to convert all xml files in the subfolder 'data' to hdf5 files and remove '.xml' from the file name
for file in os.listdir(path):
    if file.endswith(".xml"):
        f.restore(file, loglevel=0)
        if(version.parse(ct.__version__) >= version.parse('2.5.0')):
            SL0=f.velocity[0]
        else:
            SL0=f.u[0]
        egr=float(file.split('_')[0].split('egr')[-1])
        print(egr)
        phi=float(file.split('_')[1].split('phi')[-1])
        print(phi)

        index = [f.gas.species_index(specie) for specie in species]
        df = pd.concat([df, pd.DataFrame([[egr, phi, f.P, f.inlet.T, f.T[-1],SL0,flamme_thickness(f)]+list(f.X[index,-1])], columns=vartosave)]).astype(float)
        print(df)
        #f.write_hdf(path+file[:-4]+'.h5',mode='w')
        #os.remove(file)
        #print(file)
df.to_csv(os.getcwd()+'/results/'+'all_data_with_flamthick.csv',index=False)

