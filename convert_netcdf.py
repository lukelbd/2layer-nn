from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import os
import sys
import glob
base = '/project2/rossby/group07'   # for saving netcdf converted data
sys.path.append(base+'/functions/') # load Momme's functions
import numpy as np
import m_tools as MT
import m_general as M
import matplotlib.pyplot as plt
#%matplotlib inline
import imp
import xarray as xr
load_path = sys.argv[1]
key_name = load_path.split('/')[-1]
save_path=base+'/model_outputs_netCDF/'+key_name+'/'
print('Diagnostics for experiment: ',key_name)
MT.mkdirs_r(save_path)
# print('Current directory:', os.getcwd())
print('Retrieving data from: ',load_path)
print('Experiment name: ',key_name)
print('Files present: ',glob.glob(load_path+'/*'))

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
                self.data=self.load_date_perstep(self.time_datasets_list[steps[0]], self.time_real[steps])

            else:
                #print('its a list')
                print('load all timesteps that are not EKE nan')
                dataset_list=list()
                for time1, time1_real in zip([self.time_datasets_list[i] for i in steps], self.time_real[steps]):
                    dataset_list.append(self.load_date_perstep(time1, time1_real))
                    print(str(round(100*time1_real/D.time_real[-1]))+ '%')

                self.data=xr.concat(dataset_list, dim='time')

        def load_eke(self):
            ekedata=list()
            edata=xr.open_mfdataset(load_path+'eke.df')

            for time1 in self.time_datasets_list:
                ekedata.append(np.array(edata[time1].data )[0])
            return ekedata

        def load_date_perstep(self, time1, time1_real):

            dx=self.dx
            dy=self.dy
            nx=self.nx
            ny=self.ny

            # all field varaibles are zonal deviations, if not named otherwise
            u1=xr.open_mfdataset(load_path+'u1.df')[time1].data
            u2=xr.open_mfdataset(load_path+'u2.df')[time1].data
            v1=xr.open_mfdataset(load_path+'v1.df')[time1].data
            v2=xr.open_mfdataset(load_path+'v2.df')[time1].data

            q1=xr.open_mfdataset(load_path+'q1.df')[time1].data # q1
            q2=xr.open_mfdataset(load_path+'q2.df')[time1].data # q2

            u1_total=xr.open_mfdataset(load_path+'u1_total.df')[time1].data
            u2_total=xr.open_mfdataset(load_path+'u2_total.df')[time1].data

            #psi1=xr.open_mfdataset(load_path+'psi1.df')[time1].data
            #psi2=xr.open_mfdataset(load_path+'psi2.df')[time1].data

            #q1_total=xr.open_mfdataset(load_path+'q1_total.df')[time1].data # q1
            #q2_total=xr.open_mfdataset(load_path+'q2_total.df')[time1].data # q2


            q1_zm=xr.open_mfdataset(load_path+'q1_zonalmean.df')[time1].data #\bar q 1
            q2_zm=xr.open_mfdataset(load_path+'q2_zonalmean.df')[time1].data #\bar q 2

            eke=xr.open_mfdataset(load_path+'eke.df')[time1].data #\bar q 2
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

                              'q1': (['x', 'y'],   q1.T)  ,
                              'q2': (['x', 'y'],  q2.T)  ,

                              'u1_total': (['x', 'y'],  u1_total.T)  ,
                              'u2_total': (['x', 'y'],  u2_total.T)  ,

                             #'psi1': (['x', 'y'],  psi1.T) ,
                             #'psi2': (['x', 'y'],  psi2.T) ,

                              #'totalq1': (['x', 'y'],   q1_total.T)  ,
                              #'totalq2': (['x', 'y'],  q2_total.T)  ,

                              'q1_zm': (['y'], q1_zm)  ,
                              'q2_zm': (['y'], q2_zm),

                              #'eke' : (eke)
                               } ,
                                coords={'x': ('x', np.arange(0,(nx)*dx, dx )),
                                        'y': ( 'y', np.arange(0,(ny)*dy, dy )),
                                        'time': time1_real })

# In[]
D=experiment(load_path)
#time1=D.time_steps[-1]
eke=D.load_eke()
print('EKE = ', eke)

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
