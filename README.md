# A 2-Layer Forced Quasi-Geostrophic Channel Model
This model simulates a simple 2-layer system governed by quasi-geostrophic physics by solving two prognostic equations: 1) tendency of the mean background potential vorticity gradient $$\partial \bar{q}/\partial y$$ and 2) tendency of the anomalous "eddy" potential vorticity $q'$. 
<!-- integrates the equations for a simple 2-layer quasi-geostrophic equations -->
The current `make.sh` and `Makefile` compile the source code for the University of Chicago's "Midway" supercomputer.
* You can vary its width with "y_sp" (the units are "percentage of the top and bottom halves of channel covered by sponge" -- that is, if y_sp is 0.3, then both the bottom and top 15% of the y-domain are covered with sponge.
