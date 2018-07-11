from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import numpy as np
import xarray as xr


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
            self.path = path
            self.load_path=path

            # Define basic lenghts:
            s=xr.open_mfdataset(self.path+'u1.df')

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
                    print(str(round(100*time1_real/self.time_real[-1]))+ '%')

                self.data=xr.concat(dataset_list, dim='time')

        def load_eke(self):
            ekedata=list()
            edata=xr.open_mfdataset(self.path+'eke.df')

            for time1 in self.time_datasets_list:
                ekedata.append(np.array(edata[time1].data )[0])
            return ekedata

        def load_date_perstep(self, time1, time1_real):

            dx=self.dx
            dy=self.dy
            nx=self.nx
            ny=self.ny

            # all field varaibles are zonal deviations, if not named otherwise
            u1=xr.open_mfdataset(self.load_path+'u1.df')[time1].data
            u2=xr.open_mfdataset(self.load_path+'u2.df')[time1].data
            v1=xr.open_mfdataset(self.load_path+'v1.df')[time1].data
            v2=xr.open_mfdataset(self.load_path+'v2.df')[time1].data

            q1=xr.open_mfdataset(self.load_path+'q1.df')[time1].data # q1
            q2=xr.open_mfdataset(self.load_path+'q2.df')[time1].data # q2

            u1_zonalmean=xr.open_mfdataset(self.load_path+'u1_zonalmean.df')[time1].data
            u2_zonalmean=xr.open_mfdataset(self.load_path+'u2_zonalmean.df')[time1].data

            #psi1=xr.open_mfdataset(load_path+'psi1.df')[time1].data
            #psi2=xr.open_mfdataset(load_path+'psi2.df')[time1].data

            #q1_total=xr.open_mfdataset(load_path+'q1_total.df')[time1].data # q1
            #q2_total=xr.open_mfdataset(load_path+'q2_total.df')[time1].data # q2


            q1_zm=xr.open_mfdataset(self.load_path+'q1_zonalmean.df')[time1].data #\bar q 1
            q2_zm=xr.open_mfdataset(self.load_path+'q2_zonalmean.df')[time1].data #\bar q 2

            eke=xr.open_mfdataset(self.load_path+'eke.df')[time1].data #\bar q 2
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

                              'u1_zm': (['y'],  u1_zonalmean.T)  ,
                              'u2_zm': (['y'],  u1_zonalmean.T)  ,

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

class experiment_fast(object):
        def __init__(self, path, dt=None):
            """
            load experiment diagnostiscs
            returns:

            self.path
            self.diagnostics
            self.setup
            self.snapshots
            """
            self.path = path
            self.load_path=path

            # Define basic lenghts:
            s=xr.open_mfdataset(self.path+'q1.df')

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
                var_list=['u1', 'v1', 'u2', 'v2']

                print('load all timesteps that are not EKE nan')
                #for var in var_list:
                var=var_list[0]

                dataset_list=list()
                for time1, time1_real in zip([self.time_datasets_list[i] for i in steps], self.time_real[steps]):
                    dataset_list.append(self.load_date_perstep(time1, time1_real))
                    print(str(round(100*time1_real/self.time_real[-1]))+ '%')

                self.data=xr.concat(dataset_list, dim='time')

        def load_eke(self):
            ekedata=list()
            edata=xr.open_mfdataset(self.path+'eke.df')

            for time1 in self.time_datasets_list:
                ekedata.append(np.array(edata[time1].data )[0])
            return ekedata

        def load_data_pervar(self, var, time, time_real):
            dx=self.dx
            dy=self.dy
            nx=self.nx
            ny=self.ny

            data=xr.open_mfdataset(self.load_path+var+'.df')[time].data

            return xr.DataArray(data, name= var, dims={'x': ('x', np.arange(0,(nx)*dx, dx )),
                                        'y': ( 'y', np.arange(0,(ny)*dy, dy )),
                                        'time': time_real })


        def load_date_perstep(self, time1, time1_real):

            dx=self.dx
            dy=self.dy
            nx=self.nx
            ny=self.ny

            # all field varaibles are zonal deviations, if not named otherwise
            u1=xr.open_mfdataset(self.load_path+'u1.df')[time1].data
            u2=xr.open_mfdataset(self.load_path+'u2.df')[time1].data
            v1=xr.open_mfdataset(self.load_path+'v1.df')[time1].data
            v2=xr.open_mfdataset(self.load_path+'v2.df')[time1].data

            q1=xr.open_mfdataset(self.load_path+'q1.df')[time1].data # q1
            q2=xr.open_mfdataset(self.load_path+'q2.df')[time1].data # q2

            u1_zonalmean=xr.open_mfdataset(self.load_path+'u1_zonalmean.df')[time1].data
            u2_zonalmean=xr.open_mfdataset(self.load_path+'u2_zonalmean.df')[time1].data

            q1_zm=xr.open_mfdataset(self.load_path+'q1_zonalmean.df')[time1].data #\bar q 1
            q2_zm=xr.open_mfdataset(self.load_path+'q2_zonalmean.df')[time1].data #\bar q 2

            eke=xr.open_mfdataset(self.load_path+'eke.df')[time1].data #\bar q 2
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

                              'u1_zm': (['y'],  u1_zonalmean.T)  ,
                              'u2_zm': (['y'],  u1_zonalmean.T)  ,

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
