# A 2-Layer Forced Quasi-Geostrophic Channel Model
This model simulates a simple 2-layer system governed by quasi-geostrophic physics by solving two prognostic equations: 1) tendency of the mean background potential vorticity gradient "dq_bar/dy" and 2) tendency of the anomalous "eddy" potential vorticity "q". 
The current `make.sh` and `Makefile` compile the source code for the University of Chicago's "Midway" supercomputer.

## Usage
In a nutshell, to run the model, call `run.sh expname [param1=value1] [param2=value2]` (or write an `sbatch` script that calls it). The experiment name is the name of the directory in `scratch` where output will be saved (the directory will be created if it does not exist). If you call `run.sh` with just an experiment name, `run.sh` will look for the namelist modifications associated with that experiment in the file `experiments.txt`. If you call `run.sh` with any `param=value` pairs, the those assignments will be added to the namelist instead. Check out the source code for `run.sh` for more information. The run script also calls `make`, to compile any modificiations to the source code.

## Forcing
The model is forced in five ways:
  1. "Radiation": Relaxation of the upper and lower-layer winds toward an "equilibrium", spatially uniform shear (the value "shear" provided in the namelist).
  2. "Friction": Relaxation of the lower-layer relative vorticity toward zero.
  3. "Diffusion": Damping of potential vorticity proportional to its n'th spatial derivative; defaults to the 6th derivative. 
  4. Sponge: Relaxation of the lower and upper-layer relative vorticity toward zero, in the exact same way as friction.
  5. PV Injection: Application of PV *tendency* anomalies in narrow spectral band, localized to center of channel and with e-folding timescale of the autocorrelation specified by an "injection timescale".

## Namelist
The namelist file is meant to tune the  model timing, background state, and forcing schemes. See the default namelist file `input.nml` and the global variables file `global_variables.f90` for details. Here are a few explanations:

  * `y_sp` controls the proportion of each top/bottom "half" of the channel covered by sponge damping.
  * `y_i` controls the width of the pv injection region, same units as above; weighting of the injection is specified by the Gaussian weight `exp(-(y-0.5*width)^2/(y_i)^2)`.
  * `tau_r`, `tau_f`, `tau_sp` are the timescales for radiation, friction, and sponge damping in days, respectively.
  * `tau_i` is the e-folding time, in *seconds*, for the autocorrelation function underlying the pv injections.
  * `wmin_i` and `wmax_i` are the minimum and maximum integer wavenumbers for the pv injections.
  * `amp_i` is the maximum amplitude of the pv injections, in units 1/s^2 (remember we inject pv *tendencies*).
  * `shear` and `beta` control the background state.
  * `rd` is the Rossby radius of deformation; it is a function of `beta`, gravity, and layer height.
