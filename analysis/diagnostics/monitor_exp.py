from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import sys
sys.path.append('/project2/rossby/group07/functions/')
import numpy as np
import matplotlib.pyplot as plt


import xarray as xr

import time
import os.path as path

key_name='exp_both'
base='/home/t-970c07/scratch-midway2/'
#base='/project2/rossby/group07/'
#base='/Projects/mount/'
load_path=base+key_name+'/'

plot_path='/project2/rossby/group07/2layer_nn_plots/'+key_name+'/'

import os
print(os.getcwd())
import imp

# In[]

def monitor_eke(file3= None, plot=False, last_digets=10):
    file2=np.copy(file3)
    if file2 is None:
        raise Warning('file input not defined')
        file2 =1
    dt=1600

    s=xr.open_mfdataset(file2)
    time_datasets_list = [k for k in s.keys() if 'Data' in k]
    time_real=np.arange(0,(len(time_datasets_list))*dt, dt )
    time_steps=np.arange(0,(len(time_datasets_list)), 1 )

    eke=list()
    for time1 in time_datasets_list:
        eke.append(np.array( xr.open_mfdataset(file2)[time1].data )[0])

    if plot:
        plt.plot(time_real, np.array(eke))
    else:
        print('timestep: '+ str(time_real[-1]) +'  EKE:')
        print([round(i,  3) for i in eke[-last_digets::]])

    s.close()
    del s, file2

    return time_real, eke


# In[]


ffile=load_path+'eke.df'

count = 0
mondtime_old=0
eke_diff=10
while (eke_diff > 1e-4):
    mondtime=path.getmtime(ffile)
    #print(mondtime - mondtime_old)
    if (mondtime - mondtime_old) > .1:
        time_real, eke =monitor_eke(file3=ffile, plot=False, last_digets=10)
        eke_diff=eke[-1] - eke[-2]
        print('EKE diff  '+ str(eke_diff) )
        print('--')
        mondtime_old=mondtime
    time.sleep(2)
    count = count + 1

print("Good bye!")
