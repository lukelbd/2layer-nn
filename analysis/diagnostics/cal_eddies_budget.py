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
key_name='runs2_noboru_standart3'#, '2layer_highdamping'#'exp_forcing_naburu_tau10_u1_widejet'# 'exp_forcing_naburu_tau10_output_test'
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

# In[]
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
s=eddy_cal.cal_eddies(s, params)


# In[]
yaxx=(s['y']- s['y'][int(params['ny']/2)])#*params['dy']/1e5
yaxx=yaxx/yaxx[-1]

limx=.3
# test
F=M.figure_axis_xy(x_size=6, y_size=7, container=True )
S1=plt.subplot(2, 1, 1)
S1=M.subplot_routines(S1)
S1.make_clear()

plt.title('Layer 1 | '+ key_name, loc='right')
plt.plot(yaxx,  s['pvflux1'].mean(dim='time'), c='red' ,  label='PV flux')
plt.plot(yaxx,  s['interact1'].mean(dim='time'), '--', c='blue', label='Heat Flux Term')
plt.plot(yaxx, -s['momflux_div_1'].mean(dim='time'), '--', c='darkblue', label='- ddy[u`v`]')

plt.plot(yaxx, (-s['momflux_div_1']+s['interact1']).mean(dim='time') , c='blue', linewidth=1.5,  label='RHS')
plt.plot(yaxx, (-s['momflux_div_1']+s['interact1'] -s['pvflux1'] ).mean(dim='time'),  c='grey',label='Residual')

plt.plot()
plt.legend()
plt.xlim(-limx, limx)
plt.xlabel('y/L/2')
plt.ylabel('$m/s^2$')

plt.legend(loc='upper right', frameon=False)
# In[
S2=plt.subplot(2, 1, 2)
S2=M.subplot_routines(S2)
S2.make_clear()
plt.title('Layer 2', loc='right')
plt.plot(yaxx,  s['pvflux2'].mean(dim='time'), c='red' , label='PV flux')
plt.plot(yaxx, -s['interact2'].mean(dim='time'), '--', c='blue',  label='Heat Flux Term')
plt.plot(yaxx, -s['momflux_div_2'].mean(dim='time'), '--', c='darkblue',  label='- ddy[u`v`]')

plt.plot(yaxx, (-s['momflux_div_2']-s['interact2']).mean(dim='time'), c='blue', linewidth=1.5,  label='RHS')
plt.plot(yaxx, (-s['momflux_div_2']-s['interact2'] -s['pvflux2'] ).mean(dim='time'), c='grey',  label='Residual')

plt.xlabel('y/L/2 ')
plt.ylabel('$m/s^2$')
plt.xlim(-limx, limx)
plt.legend(loc='lower right', frameon=False)

F.save_light(name=key_name+'_budget', path= plot_path )
F.save_pup(name=key_name+'_budget', path= plot_path )


# In[]
seke=MT.pickle_load(path=load_path, name=key_name+'_EKE')
# plot eke
F=M.figure_axis_xy(x_size=6, y_size=3.5)

F.make_clear()

plt.plot(seke[0:20])

plt.title('EKE  ' + key_name, loc='left')
plt.ylabel('EKE')
F.save_light(name=key_name+'_EKE_process', path= plot_path )


#plt.legend(loc='upper right', frameon=False)
# In[
