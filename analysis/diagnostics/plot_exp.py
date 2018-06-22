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


key_name='exp_control_naburu2_jet'
print('diagnose ' +key_name)

plot=True
#key_name='exp_both'
base='/home/t-970c07/scratch-midway2/'
#base='/project2/rossby/group07/'
#base='/Projects/mount/'
load_path=base+key_name+'/'

plot_path='/project2/rossby/group07/2layer_nn_plots/'+key_name+'/'

#save_path='/project2/rossby/group07/model_outputs_netCDF/'+key_name+'/'
#MT.mkdirs_r(save_path)
MT.mkdirs_r(plot_path)
#save_path=base+'/netCDF_files/'
print(os.getcwd())

# In[]
class experiment(object):
        def __init__(self, path, dt=None):
            """
            load experiment diagnostiscs
            returns:

            self.path
            self.diagnostics
            self.setup
            self.snapshots
            """
            self.path = load_path

            # Define basic lenghts:
            s=xr.open_mfdataset(load_path+'u1.df')

            dt=1600 if dt is None else dt

            self.dx=1
            self.dy=1
            self.nx=s['fakeDim1'].size
            self.ny=s['fakeDim0'].size


            self.time_datasets_list = [k for k in s.keys() if 'Data' in k]
            self.time_real=np.arange(0,(len(self.time_datasets_list))*dt, dt )
            self.time_steps=np.arange(0,(len(self.time_datasets_list)), 1 )


        def load_data(self, steps=None):

            steps=self.time_steps if steps is None else steps
            if type(steps) is not np.ndarray:
                steps =[steps]

            if len(steps) == 1:
                print('load single timestep')
                self.data=self.load_date_perstep(D.time_datasets_list[steps[0]], D.time_real[steps])

            else:
                #print('its a list')
                print('load all timesteps')
                dataset_list=list()
                for time1, time1_real in zip([D.time_datasets_list[i] for i in steps], D.time_real[steps]):
                    dataset_list.append(self.load_date_perstep(time1, time1_real))
                    print(str(round(100*time1_real/D.time_real[-1]))+ '%')

                self.data=xr.concat(dataset_list, dim='time')

        def load_date_perstep(self, time1, time1_real):

            dx=self.dx
            dy=self.dy
            nx=self.nx
            ny=self.ny

            u1=xr.open_mfdataset(load_path+'u1.df')[time1].data
            u2=xr.open_mfdataset(load_path+'u2.df')[time1].data
            v1=xr.open_mfdataset(load_path+'v1.df')[time1].data
            v2=xr.open_mfdataset(load_path+'v2.df')[time1].data

            #totalu1=xr.open_mfdataset(load_path+'u1_total.df')[time1].data
            #totalu2=xr.open_mfdataset(load_path+'u2_total.df')[time1].data

            #psi1=xr.open_mfdataset(load_path+'psi1.df')[time1].data
            #psi2=xr.open_mfdataset(load_path+'psi2.df')[time1].data

            q1_total=xr.open_mfdataset(load_path+'q1_total.df')[time1].data # q1
            q2_total=xr.open_mfdataset(load_path+'q2_total.df')[time1].data # q2

            zonalmean_q1=xr.open_mfdataset(load_path+'q1_zonalmean.df')[time1].data #\bar q 1
            zonalmean_q2=xr.open_mfdataset(load_path+'q2_zonalmean.df')[time1].data #\bar q 2

            # print(np.arange(0,(nx)*dx, dx ).shape, np.arange(0,(ny)*dy, dy ).shape )
            #
            # print(time1_real.shape)
            # print('total u1' +str(totalu1.T.shape) )
            # print('total u2' +str(totalu2.T.shape) )
            # print('q1 total' +str(q1_total.T.shape) )
            #
            # print('psi1' +str(psi1.T.shape) )

            return xr.Dataset({'u1': (['x', 'y'],  u1.T),
                             'u2': (['x', 'y'],  u2.T)  ,

                             'v1': (['x', 'y'],  v1.T)  ,
                             'v2': (['x', 'y'],  v2.T)  ,

                              #'totalu1': (['x', 'y'],  totalu1.T)  ,
                              #'totalu2': (['x', 'y'],  totalu2.T)  ,

                             #'psi1': (['x', 'y'],  psi1.T) ,
                             #'psi2': (['x', 'y'],  psi2.T) ,

                              'totalq1': (['x', 'y'],   q1_total.T)  ,
                              'totalq2': (['x', 'y'],  q2_total.T)  ,

                              'zonalmeanq1': (['y'],  zonalmean_q1)  ,
                              'zonalmeanq2': (['y'],  zonalmean_q2),
                               } ,
                                coords={'x': ('x', np.arange(0,(nx)*dx, dx )),
                                        'y': ( 'y', np.arange(0,(ny)*dy, dy )),
                                        'time': time1_real })
# #In[]
D=experiment(load_path)
#time1=D.time_steps[-1]
#D.time_datasets_list
D.load_data(steps=D.time_steps[-1])


# In[]
#D.data.to_netcdf(save_path+key_name+'netcdf4.nc', engine='netcdf4')
#D.data.to_netcdf(save_path+key_name+'h5netcdf.nc', engine='h5netcdf')
#D.data.to_netcdf(save_path+key_name+'scipy.nc', engine='scipy')

# In[]

if plot:
    for key, data in D.data.sel(time=D.data.time[-1]).data_vars.items():
        pdata=data

        if len(pdata.shape) ==1:
            F=M.figure_axis_xy(view_scale=.3)
            plt.semilogy(D.data.y, pdata)
            plt.title(key)
            F.save_light(path=plot_path, name='last_timestep_'+key)
        else:
            F=M.figure_axis_xy(x_size=3.2, y_size=3, view_scale=.4)
            plt.contourf(D.data.x, D.data.y,  pdata.T - pdata.T.mean(axis=1))
            plt.colorbar()
            plt.title(key + ' Zonal Mean Anomaly')
            F.save_light(path=plot_path, name='last_timestep_anomaly'+key)

            F=M.figure_axis_xy(x_size=3.2, y_size=3, view_scale=.4)
            plt.contourf(D.data.x, D.data.y,  pdata.T)
            plt.colorbar()
            plt.title(key + ' Total')
            F.save_light(path=plot_path, name='last_timestep_'+key)

#
# # In[]
# plt.figure()
# plt.plot(t, ke)
# plt.xlabel('time'); plt.ylabel('KE'); plt.title('KE')
# plt.savefig(plot_path+key_name+'_KE'+'.png', bbox_inches='tight')
# # In[]
#
#
# plt.figure()
# q_upper = m.q[0] + m.Qy[0]*m.y
# plt.contourf(m.x, m.y, q_upper, 12, cmap='RdBu_r')
# plt.xlabel('x'); plt.ylabel('y'); plt.title('Upper Layer PV')
# plt.colorbar();
# plt.savefig(plot_path+key_name+'_UpperPV'+'.png', bbox_inches='tight')
#
# plt.figure()
# q_lower = m.q[1] + m.Qy[1]*m.y
# plt.contourf(m.x, m.y, q_upper, 12, cmap='RdBu_r')
# plt.xlabel('x'); plt.ylabel('y'); plt.title('Lower Layer PV')
# plt.colorbar();
# plt.savefig(plot_path+key_name+'_LowerPV'+'.png', bbox_inches='tight')
#
# plt.figure()
# speed_upper = np.sqrt(m.u[0]**2, m.v[0]**2)
# plt.contourf(m.x, m.y, speed_upper, 12, cmap='RdBu_r')
# plt.xlabel('x'); plt.ylabel('y'); plt.title('Upper Layer Speed')
# plt.colorbar();
# plt.savefig(plot_path+key_name+'_UpperSpeed'+'.png', bbox_inches='tight')
