from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division


import sys
sys.path.append('/project2/rossby/group07/functions/')
import numpy as np
import xarray as xr


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

def cal_psi_from_u(u, dy):
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

# In[]
def cal_eddies(s, params):
    """ This script calcualates the Eddy budgets for the 2-Layer GCM"""

    for key,val in params.items():
            exec(key + '=val')

    s['psi1']=cal_psi_from_v(s['v1'], dx)
    s['psi2']=cal_psi_from_v(s['v2'], dx)
    #s['psi1']=cal_psi_from_u(s['u1_total'], dy)
    #s['psi2']=cal_psi_from_u(s['u2_total'], dy)


    # eddy fluxes
    #make primes
    pvars=['psi1', 'psi2']
    #s['psi1'].mean(dim='x').T.plot()
    for v in pvars:
        print(v+'prime')
        s[v+'prime']=s[v]-s[v].mean(dim='x')

    # eddy momentum flux
    s['up_vp_bar_1']=(s['u1']*s['v1']).mean(dim='x')
    s['up_vp_bar_2']=(s['u2']*s['v2']).mean(dim='x')

    #s['momflux_div_1']=s['up_vp_bar_1'].copy
    _, momfluxdiv=np.gradient(s['up_vp_bar_1'], 1)
    s['momflux_div_1']=xr.DataArray((momfluxdiv/dy).T, dims={'y': ( 'y', s['y']) ,'time': ('time', s['time']) })

    _, momfluxdiv=np.gradient(s['up_vp_bar_2'], 1)
    s['momflux_div_2']=xr.DataArray((momfluxdiv/dy).T, dims={'y': ( 'y', s['y']) ,'time': ('time', s['time']) })

    # interaction terms
    s['interact1']=( s['v1']*  (s['psi2prime']-s['psi1prime'])  ).mean(dim='x')/ Ld**2
    s['interact2']=( s['v2']*  (s['psi2prime']-s['psi1prime'])  ).mean(dim='x')/ Ld**2

    # PV fluxes
    s['pvflux1']=(s['v1']*s['q1']).mean(dim='x')
    s['pvflux2']=(s['v2']*s['q2']).mean(dim='x')

    return s

def plot_eke(data, key_name,  plot_path=None):
    import m_tools as MT
    import m_general as M
    import matplotlib.pyplot as plt
    F=M.figure_axis_xy(x_size=6, y_size=3.5)
    F.make_clear()

    plt.plot(data, color='black')

    plt.title('EKE  ' + key_name, loc='left')
    plt.ylabel('EKE')
    if plot_path is not None:
        F.save_light(name=key_name+'_EKE_process', path= plot_path )


def plot_budgets(s, key_name, params, plot_path=None):
    import m_tools as MT
    import m_general as M
    import matplotlib.pyplot as plt
    yaxx=(s['y']- s['y'][int(params['ny']/2)])#*params['dy']/1e5
    yaxx=yaxx/yaxx[-1]

    limx=.3
    # test
    F=M.figure_axis_xy(x_size=6, y_size=7, container=True )
    S1=plt.subplot(2, 1, 1)
    S1=M.subplot_routines(S1)
    S1.make_clear()

    plt.title('Layer 1 | '+ key_name, loc='right')
    plt.plot(yaxx,  s['pvflux2'].mean(dim='time')*0, c='black', linewidth=.5)

    plt.plot(yaxx,  s['pvflux1'].mean(dim='time'), c='red' ,  label='PV flux')
    plt.plot(yaxx,  s['interact1'].mean(dim='time'), '--', c='blue', label='Interaction Term')
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
    plt.plot(yaxx,  s['pvflux2'].mean(dim='time')*0, c='black', linewidth=.5)

    plt.plot(yaxx,  s['pvflux2'].mean(dim='time'), c='red' , label='PV flux')
    plt.plot(yaxx, -s['interact2'].mean(dim='time'), '--', c='blue',  label='Interaction Term')
    plt.plot(yaxx, -s['momflux_div_2'].mean(dim='time'), '--', c='darkblue',  label='- ddy[u`v`]')

    plt.plot(yaxx, (-s['momflux_div_2']-s['interact2']).mean(dim='time'), c='blue', linewidth=1.5,  label='RHS')
    plt.plot(yaxx, (-s['momflux_div_2']-s['interact2'] -s['pvflux2'] ).mean(dim='time'), c='grey',  label='Residual')

    plt.xlabel('y/L/2 ')
    plt.ylabel('$m/s^2$')
    plt.xlim(-limx, limx)
    plt.legend(loc='lower right', frameon=False)

    if plot_path is not None:
        F.save_light(name=key_name+'_budget', path= plot_path )
        F.save_pup(name=key_name+'_budget', path= plot_path )
