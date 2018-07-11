from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import os
import sys
import glob
base = '/project2/rossby/group07'   # for saving netcdf converted data
sys.path.append(base+'/functions/') # load Momme's functions
import cannelmodel_hdf_to_netCDF as data_converter

import numpy as np
import m_tools as MT
import m_general as M
import matplotlib.pyplot as plt
#%matplotlib inline
#import xarray as xr
import imp


load_path = sys.argv[1]
key_name = load_path.split('/')[-1]
load_path += '/' # add slash to end
save_path=base+'/model_outputs_netCDF/'+key_name+'/'
print('Diagnostics for experiment: ',key_name)
MT.mkdirs_r(save_path)
# print('Current directory:', os.getcwd())
print('Retrieving data from: ',load_path)
print('Experiment name: ',key_name)
# print('Files present: ',glob.glob(load_path+'/*'))

# In[]
D=data_converter.experiment(load_path)

#time1=D.time_steps[-1]
eke=D.load_eke()
print('EKE = ', eke)


MT.pickle_save(path=save_path, name=key_name+'_EKE', data=eke, verbose=True)
MT.save_log_txt(path=save_path, name=key_name+'_EKE', hist=eke, verbose=True)


# In[]
#plt.contourf(u1_total)
D.load_data(steps=D.time_steps)

#D.load_data(steps=D.time_steps[~np.isnan(eke)])
# # In[]
#s=D.data#.sel(time=D.data.time[2])
#
#plt.contourf(s['u1_total'])
# s['v1'].mean(dim='x').plot()
# s['q1'].mean(dim='x').plot()
#
#
# plt.contourf(s['q1'].T)
# In[]
#D.data.to_netcdf(save_path+key_name+'netcdf4.nc', engine='netcdf4')
#D.data.to_netcdf(save_path+key_name+'h5netcdf.nc', engine='h5netcdf')
D.data.to_netcdf(save_path+key_name+'scipy.nc', engine='scipy')
print('Saved at: ',save_path+key_name+'scipy.nc')
