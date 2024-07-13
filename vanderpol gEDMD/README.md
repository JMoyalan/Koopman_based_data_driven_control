## Running the code
1. Please ensure that both the folders i.e., `SOS Toolbox` and `vanderpol gEDMD` are in the same location on your PC to solve the convex optimization using `sosopt` toolbox.
2. Run the file `vanderpol gEDMD/sampling/main_VDP_sampling.m` to generate the time-series data from the Van der Pol dynamics.
3. Run the file `vanderpol gEDMD/gEDMD/main_gedmdc.m` to obtain the Koopman generator matrix for the Van der Pol dynamics.
4. Run the file `vanderpol gEDMD/L1OCP/main_VDP_L1OCP.m` to solve the convex optimization to find the optimal control to stabilize the dynamics with respect to $`\mathcal{L}_1`$-norm in the cost function. Alternatively, run the file `vanderpol gEDMD/L2OCP/main_VDP_L2OCP.m` to solve the convex optimization with respect to $`\mathcal{L}_2`$-norm in the cost function.

## Results

<p align="center">
<img src="figures/VDP_results.PNG" width="1000">
</p>

