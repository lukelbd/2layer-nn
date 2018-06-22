from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import os
import sys
sys.path.append('/project2/rossby/group07/functions/')
import numpy as np
import m_tools as MT
import m_general as M
import os
import matplotlib.pyplot as plt
#%matplotlib inline
import imp

import xarray as xr
import cannelmodel_hdf_to_netCDF as data_converter


if len(sys.argv) <= 1:
    key_name='exp_control_naburu'
elif sys.argv[1] == '-f':
    key_name='exp_control_naburu'
else:
    key_name=sys.argv[1]

print('diagnose ' +key_name)
#key_name='2layer_long_run_2'
key_name='exp_forcing_naburu_tau10_output_test'
base='/scratch/midway2/holo/'
#base='/home/t-970c07/scratch-midway2/'
#base='/project2/rossby/group07/'
#base='/Projects/mount/'
load_path=base+key_name+'/'

save_path='/project2/rossby/group07/model_outputs_netCDF/'+key_name+'/'
MT.mkdirs_r(save_path)

#save_path=base+'/netCDF_files/'
print(os.getcwd())
print(load_path)
# In[]

# In[]

#u1_zonalmean=xr.open_mfdataset(load_path+'u1_zonalmean.df')

D=data_converter.experiment(load_path)
#time1=D.time_steps[-1]

eke=D.load_eke()
print('EKE = ' )
print(eke)

MT.pickle_save(path=save_path, name=key_name+'_EKE', data=eke, verbose=True)
MT.save_log_txt(path=save_path, name=key_name+'_EKE', hist=eke, verbose=True)
# In[]

var_list=['u1', 'v1', 'u2', 'v2']
var=var_list[0]
self=D
dx=self.dx
dy=self.dy
nx=self.nx
ny=self.ny

i=1
time=D.time_datasets_list# for i in steps],
self.time_real[steps]

data_all=xr.open_mfdataset(self.load_path+var+'.df')

for t in time:
    data_single=data_all[t].data
[time]

timearray=xr.DataArray(self.time_real, name= 'time')
ds=xr.DataArray(data_single, name= var, dims={'x': ('x', np.arange(0,(nx)*dx, dx )),
                            'y': ( 'y', np.arange(0,(ny)*dy, dy )),
                            })#'time': self.time_real[0]
#ds.to_dataset(name=var)
#ds=ds.expand_dims('time')

#ds['time']=              self.time_real[0]

ds

#help(xr.auto_combine)
dd=xr.concat([ds, ds, ds, ds], timearray[0:4])

dd.to_dataset(name=var)['time']
['u1']
#plt.contourf(u1_total)
#D.load_data(steps=D.time_steps[::10])

#D.data['u1_zm'].mean(dim='time').plot()
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

print('saved at:')
print(save_path+key_name+'scipy.nc')
