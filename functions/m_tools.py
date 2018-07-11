from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

# import matplotlib
# #matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import numpy as np
#from scipy.io import netcdf
#from mpl_toolkits.basemap import Basemap,cm
import os
import pickle
#from matplotlib.dates import DateFormatter, MinuteLocator
#from matplotlib import dates
#import datetime as DT


def dt_form_timestamp(timestamp, unit=None):
    unit='h' if unit is None else unit
    return (timestamp[1]-timestamp[0]).astype('m8['+unit+']')


def tick_formatter(a, interval=2, rounder=2, expt_flag=True):

    O=int(np.log10(a.max()))
    fact=10**(O-1)
    b=np.round(a/fact, rounder+1)*fact
    ticklabels=[' ' for i in range(len(b))]
    N=int(np.ceil(len(b)/interval))

    tt=np.arange(0,len(b),interval)

    for t in tt:
        if expt_flag:
            ticklabels[int(t)]='{:.2e}'.format(b[t])
        else:

            ticklabels[int(t)]=format(b[t], '.2f').rstrip('0').rstrip('.')#'{:.{2}f}'.format(b[t])

    #ticks=a
    return ticklabels, b


def freq_lim_string(low, high):
	a='%2.1e' % low
	b='%2.1e' % high

	return a[0:3] +'-'+ b +' Hz'

def mkdirs_r(path):
    if not os.path.exists(path):
                os.makedirs(path)

def check_year(inputstr, yearstring):
    a=np.datetime64(inputstr).astype(object).year
    ref=np.datetime64(yearstring).astype(object).year
    if a == ref:
        return True
    else:
        return False

def datetime64_to_sec(d):
    return d.astype('M8[s]').astype('float')
def datetime64_to_day(d):
    return d.astype('M8[D]').astype('float')

def float_plot_time_to_sec(pp):
    return np.datetime64(dates.num2date(pp)).astype('M8[s]').astype('float')
def float_plot_time_to_dt64(pp):
    return np.datetime64(dates.num2date(pp)).astype('M8[s]')

def sec_to_dt64(pp):
    return pp.astype('M8[s]')

def sec_to_float_plot(pp):
    from matplotlib import dates
    import datetime as DT
    return dates.date2num(pp.astype('M8[s]').astype(DT.datetime))
def sec_to_float_plot_single(pp):
    from matplotlib import dates
    import datetime as DT
    return dates.date2num(np.datetime64(int(pp), 's').astype('M8[s]').astype(DT.datetime))
def fake_2d_data(verbose=True, timeaxis=False):
    x=np.arange(0,100,1)
    y=np.arange(0,40,1)
    XX, YY= np.meshgrid(x,y)

    mu=x.size/2
    sigma=x.size/5
    z2= 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (XX - mu)**2 / (2 * sigma**2) )
    z2=z2/z2.max()

    mu=y.size/2
    sigma=y.size/5
    z3= 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (YY - mu)**2 / (2 * sigma**2) )
    z3=z3/z3.max()
    if verbose:
        print('x' , x.shape)
        print('y' , y.shape)
        print('z' , z3.shape)

        plt.contourf(x, y,z2/2+z3/2)
        plt.colorbar()
        plt.axis('scaled')
        plt.show()


    return x, y, z3

def pickle_save(name, path, data, verbose=True):
    if not os.path.exists(path):
        os.makedirs(path)
    full_name= (os.path.join(path,name+ '.npy'))


    with open(full_name, 'wb') as f2:
        pickle.dump(data, f2)
    if verbose:
        print('save at: ',full_name)

def pickle_load(name, path, verbose=True):
    #if not os.path.exists(path):
    #    os.makedirs(path)
    full_name= (os.path.join(path,name+ '.npy'))

    with open(full_name, 'r') as f:
        data=pickle.load(f)

    if verbose:
        print('load from: ',full_name)
    return data


def json_save(name, path, data, verbose=False, return_name=False):
    import json
    if not os.path.exists(path):
        os.makedirs(path)
    full_name_root=os.path.join(path,name)
    full_name= (os.path.join(full_name_root+ '.json'))
    with open(full_name, 'w') as outfile:
        json.dump(data, outfile)
    if verbose:
        print('save at: ',full_name)
    if return_name:
        return full_name_root
    else:
        return


def json_load(name, path, verbose=False):
    import json
    full_name= (os.path.join(path,name+ '.json'))

    with open(full_name, 'r') as ifile:
        data=json.load(ifile)
    if verbose:
        print('loaded from: ',full_name)
    return data

def h5_load(name, path, verbose=False):
    import pandas as pd
    full_name= (os.path.join(path,name+ '.h5'))
    data=pd.read_hdf(full_name)
    #with pd.HDFStore(full_name) as data:
    #    data = pd.HDFStore(path+self.ID+'.h5')
    #    data.close()
    return data

def write_log(hist, string, verbose=False, short=True , date=True):
    import datetime as datetime
    if short:
        now = datetime.datetime.now().strftime("%Y%m%d")
    else:
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    if date:
        message='\n'+now+' '+string
    else:
        message='\n '.ljust(16)+' '+string

    if verbose== True:
        print(message)
    elif verbose == 'all':
        print(hist+message)
    return hist+message

def save_log_txt(name, path, hist,verbose=False):
    if not os.path.exists(path):
        os.makedirs(path)
    full_name= (os.path.join(path,name+ '.hist.txt'))
    with open(full_name, 'w') as ifile:
        ifile.write(str(hist))
    if verbose:
        print('saved at: ',full_name)

def load_log_txt(name, path):
    import glob
    hist_file=name#'DR01.LHN.stormdetect.A02_geometry_cut_storm.hist.txt' #ID.string+'.A02**.txt'
    f=[]
    for h in glob.glob(os.path.join(path,hist_file)):
        f.append(open(h, 'r').read())
    return '\n'.join(f)

def shape(a):
	for i in a:
		print(i.shape)

def find_O(a, case='round'):
    if case=='round':
        for k in np.logspace(0,24,25):
            if np.ceil(a/k) == 1:
                return k
                break
    elif case=='floor':
        for k in np.logspace(0,24,25):
            if np.ceil(a/k) == 1:
                return k
                break

    elif case=='ceil':
        for k in np.logspace(0,24,25):
            if np.ceil(a/k) == 1:
                return k
                break
    else:
        raise Warning('no propper case')

def stats(a):
	print('shape' , a.shape)
	print('Nans',np.sum(np.isnan(a)))
	print('max' , np.nanmax(a))
	print('min' ,np.nanmin(a))
	print('mean' ,np.nanmean(a))

def stats_format(a, name=None):
	print('Name:', str(name),'   Shape:' , a.shape ,'   NaNs:',np.sum(np.isnan(a)),' max:', np.nanmax(a),' min', np.nanmin(a),' mean:', np.nanmean(a))
