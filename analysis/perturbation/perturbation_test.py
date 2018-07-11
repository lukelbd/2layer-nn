import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline



x=np.arange(0, 265)
y=np.arange(0, 265)
xx, yy= np.meshgrid(x,y )

y_p=y- y[-1]/2

# Parameters
A0=  1          # amplitude
kappa0 = .2     # center of narrow banded Process
sigma_kappa = .5  # width of narrow banded Process

y_m = 0         # center of the mask
sigma_y = 20       # width of the mask

tau_i = 5 # days       #  injection rate

# mask
mask= np.exp( - abs(y_p - y_m)/ sigma_y )
plt.plot(mask)
mask_rep=np.repeat(np.array([mask]), x.size, axis=0).T
#plt.contourf(mask_rep)
# amplitude
def amp(A0, kappa0, sigma_kappa, k,l):

    kappa= np.sqrt(k**2 + l**2)
    return A0 *np.exp( - ( kappa- kappa0)**2 / sigma_kappa**2 )


li=.1
ki=.2
phi_l_i=np.pi/2
phi_k_i=np.pi/3

A_kappa=amp(A0, kappa0, sigma_kappa , .1, .1)

def field(li, ki, phi_l_i,phi_k_i ):
    return np.sin(li *yy  + phi_l_i) * np.real(  np.exp(  xx*0 +1j* ki *xx +phi_k_i))

F= mask_rep *field(li, ki, phi_l_i,phi_k_i )* A_kappa
plt.contourf(F)

# randomness
n_wave=10
wave_min=.2
wave_max=.7

np.random.uniform(low=wave_min, high=wave_max, size=n_wave)
