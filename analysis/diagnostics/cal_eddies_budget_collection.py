from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import sys
sys.path.append('/project2/rossby/group07/functions/')
import numpy as np
import m_tools as MT
import m_general as M
import matplotlib.pyplot as plt
%matplotlib inline
import xarray as xr
import os
import eddy_calculator as eddy_cal
import imp

if len(sys.argv) <= 1:
    key_name='exp_control_naburu'
elif sys.argv[1] == '-f':
    key_name='exp_control_naburu'
else:
    key_name=sys.argv[1]
key_name='runs2_famp_1e-9'
#'runs2_famp_1e-9'#'exp_forcing_naburu_tau10_u1_widejet'#, '2layer_highdamping'#'exp_forcing_naburu_tau10_u1_widejet'# 'exp_forcing_naburu_tau10_output_test'
print('diagnose ' +key_name)

#key_name='exp_both'
#base='/home/t-970c07/scratch-midway2/'
#base='/project2/rossby/group07/'
#base='/Projects/mount/'
#load_path=base+key_name+'/'


load_path='/project2/rossby/group07/model_outputs_netCDF/'+key_name+'/'
plot_path='/project2/rossby/group07/2layer_nn_plots/'+key_name+'/'
#MT.mkdirs_r(save_path)

#save_path=base+'/netCDF_files/'
print(os.getcwd())

# In[] load example

s_load=xr.open_mfdataset(load_path+key_name+'scipy.nc')
s=s_load.sel(time=s_load['time'])

params=dict()

params['L']=28000000.0
params['W']=72000000.0
params['Ld']=800000
params['dx']=params['L']/len(s['x'])
params['dy']=params['W']/len(s['y'])
#np.cumsum( s['u1'][0,:,:].data, axis=1)
params['x']=s['x']
params['y']=s['y']
params['X']=params['x']*params['dx']
params['Y']=params['y']*params['dy']

params['nx']=s['x'].shape[0]
params['ny']=s['y'].shape[0]

# In[]

u0s=[0, 0.5, 1, 2, 5, 10] # background shear
famps=[1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7] # forcing amplitude
tau_fs=[0.1, 0.2, 0.5, 1, 2, 5, 10] # friction damping


for famp in famps:
    print('runs2_'+'famp_'+str(famp))
    key_name='runs2_'+'famp_'+str(famp)

    load_path='/project2/rossby/group07/model_outputs_netCDF/'+key_name+'/'
    plot_path='/project2/rossby/group07/2layer_nn_plots/'+key_name+'/'

    s_load=xr.open_mfdataset(load_path+key_name+'scipy.nc')
    s=s_load.sel(time=s_load['time'])

    s=eddy_cal.cal_eddies(s, params)
    eddy_cal.plot_budgets(s, key_name, params,  plot_path=plot_path)

    seke=MT.pickle_load(path=load_path, name=key_name+'_EKE')
    eddy_cal.plot_eke(seke, key_name,  plot_path=plot_path)

for u0 in u0s:
    print('u0_'+str(u0))
    key_name='runs2_'+'u0_'+str(u0)

    load_path='/project2/rossby/group07/model_outputs_netCDF/'+key_name+'/'
    plot_path='/project2/rossby/group07/2layer_nn_plots/'+key_name+'/'

    s_load=xr.open_mfdataset(load_path+key_name+'scipy.nc')
    s=s_load.sel(time=s_load['time'])

    s=eddy_cal.cal_eddies(s, params)
    eddy_cal.plot_budgets(s, key_name, params,  plot_path=plot_path)

    seke=MT.pickle_load(path=load_path, name=key_name+'_EKE')
    eddy_cal.plot_eke(seke, key_name,  plot_path=plot_path)


for tau_f in tau_fs:
    print('tau_f_'+str(tau_f))
    key_name='runs2_'+'tau_f_'+str(tau_f)

    load_path='/project2/rossby/group07/model_outputs_netCDF/'+key_name+'/'
    plot_path='/project2/rossby/group07/2layer_nn_plots/'+key_name+'/'

    s_load=xr.open_mfdataset(load_path+key_name+'scipy.nc')
    s=s_load.sel(time=s_load['time'])

    s=eddy_cal.cal_eddies(s, params)
    eddy_cal.plot_budgets(s, key_name, params,  plot_path=plot_path)

    seke=MT.pickle_load(path=load_path, name=key_name+'_EKE')
    eddy_cal.plot_eke(seke, key_name,  plot_path=plot_path)
