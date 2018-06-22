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
%matplotlib inline
import imp
import xarray as xr


if len(sys.argv) <= 1:
    key_name='exp_control_naburu'
elif sys.argv[1] == '-f':
    key_name='exp_control_naburu'
else:
    key_name=sys.argv[1]
key_name='exp_control_naburu2'
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


s=xr.open_mfdataset(load_path+key_name+'scipy.nc')

tstep=3
L=28000000.0
W=72000000.0
dx=L/len(s['x'])
dy=W/len(s['y'])
#np.cumsum( s['u1'][0,:,:].data, axis=1)
x=s['x']
y=s['y']
X=x*dx
Y=y*dy

xx, yy =np.meshgrid(np.sin(np.pi*X/L), np.sin(np.pi*Y/W))
tpsi=xx*yy
plt.figure()

plt.contourf(X, Y, tpsi)
plt.title('psi')
plt.colorbar()

tpsi=s['u1'][1,:,:]
tpsi.data=(xx*yy).T


tv, tu=np.gradient(tpsi, 1)
tu=- tu/dy
tv= tv/dx
plt.figure()
plt.subplot(2, 1, 1)
plt.contourf(X, Y,tu.T)
plt.colorbar()

plt.subplot(2, 1, 2)
plt.contourf(X, Y,tv.T)
plt.colorbar()

# In[] test call

tv_xr=tpsi
tv_xr.data=tv
tv_psi=tv_xr.cumsum(dim='x')*dx

_, tv_psi_u=np.gradient(tv_psi, 1)
tv_psi_u=- tv_psi_u/dy

plt.figure()
plt.subplot(2, 1, 1)
plt.contourf(X, Y,tv_psi.T)
plt.colorbar()
plt.title('psi from tv')

plt.subplot(2, 1, 2)
plt.contourf(X, Y,tv_psi_u.T)
plt.colorbar()

# In[]

def cal_psi_from_v(v, dx):
    """
    v is an xarray
    dx is the grid point distance (m)

    returns
    psi as an xarray
    """
    v_psi=v.cumsum(dim='x')*dx
    v_psi.name='psi'
    return v_psi

def cal_psi_from_v(u, dy):
    """
    u is an xarray
    dy is the grid point distance (m)

    returns
    psi as an xarray
    """
    u_psi=-u.cumsum(dim='y')*dy
    u_psi.name='psi'
    return u_psi

def cal_u_from_psi(psi, dy):
    """
    psi is an xarray
    dy is the grid point distance (m)

    returns
    u as an np.array
    """
    _, psi_u=np.gradient(psi, 1)
    psi_u=- psi_u/dy
    #tv_psi_u.name='u'
    return psi_u

s=xr.open_mfdataset(load_path+key_name+'scipy.nc')

v_empl=s['v1'][tstep,:,:]
psi_cal=cal_psi_from_v(v_empl, dx)

u_cal=cal_u_from_psi(psi_cal, dy)


# In[]
plt.figure()
plt.contourf(v_empl.T)
plt.title('v1')

plt.figure()
plt.contourf(psi_cal.T)
plt.title('psi')

plt.figure()
plt.contourf(u_cal.T)
plt.title('psiy')
plt.colorbar()

plt.figure()
plt.contourf(s['u1'][tstep,:,:].T)
plt.title('u')
plt.colorbar()

plt.figure()
plt.contourf(u_cal.T- s['u1'][tstep,:,:].T)
plt.title('psiy - u')
plt.colorbar()

# In[]

.trapz
s['u1']
