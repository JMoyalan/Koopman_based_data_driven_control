## Running the code
1. Please ensure that both the folders i.e., `SOS Toolbox` and `Lorentz attractor` are in the same location on your PC to solve the convex optimization using `sosopt` toolbox.
2. Run the file `Lorentz attractor/sampling/main_Lorentz_sampling.m` to generate the time-series data from the Lorentz attractor dynamics.
3. Run the file `Lorentz attractor/Lorentz_main_L1OCP.m` to solve the convex optimization to find the optimal control to stabilize the dynamics with respect to $`\mathcal{L}_1`$-norm in the cost function. Alternatively, run the file `Lorentz attractor/Lorentz_main_L2OCP.m` to solve the convex optimization with respect to $`\mathcal{L}_2`$-norm in the cost function.

## Results
The following figure showcases the optimal trajectories of the Lorentz attractor converging to $`\boldsymbol{X}_T=[0,0,0]`$. This is obtained by solving for the optimal control in the space of density with $`\mathcal{L}_1`$-norm in the cost function. 
<p align="center">
<img src="../figures/Lorentz_results.PNG" width="500">
</p>
