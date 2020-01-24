# Characterizing_noise

## This repository contains codes and data used for analysis in the article "Noise-Induced Effects in Collective Dynamics and Inferring Interactions from Data".

# Codes for Figure 2

1. **Generate polarisation time-series for a given set of parameters:** A matlab code (__/Figure_2/Gillespie_stochastic_process.m__) can be used to compute the polarization/order parameter.
	1. The important parameters that need to be set are N (system size), r1, r2, r3, r4 (which are for different reaction rates).
		1. N = 50, 100, 200; r1 = 0.01; r2 = 1; r3 = 0 (for pairwise model) & 0.08 (for ternay model)
	2. Other parameters include Tint, Tend.
		1. Tint = 50; 
		2. Tend = Tint*Number of iterations. Set number of iterations in this formula to 1000000.
	3. Parameter rel is for number of realizations/replicates of the simulations. 
		1. Set rel >= 1. For 1 realization with 1 million iterations may take around 10-15 mins. **If the system is out of memory reduce the number of iterations.**
	4. The output is a time series array (S) of size (number of iterations*rel), and an array (tSample) storing system time.
	5. Autocorrelation time of the time series (est_tau) can be calculated. exp_tau is the expected autocorrealtion time for the given reaction rate and system size parameters.
	
2. **Reconstruct Deterministic and Stochastic part from time series:** A matlab code (__/Figure_2/SDE_different_Dt.m__) can be used to reconstruct the underlying functions from the time series for different time scales. This code can generate all plots in Figure 2 of the main text of the paper.
	1. Run one column of the time series array (for e.g., S(:,1)) at a time on this code.
	2. In line 8: Set different time scales over which you want to calculate the deterministic and the stochastic part.
	3. Uncomment line 33 and set the exact parameters used to generate the time series using code __/Figure_2/Gillespie_stochastic_process.m__
	4. The code uses tSample from the code __/Figure_2/Gillespie_stochastic_process.m__ (or from real data) for calculaing Tint in line 35. One can also set Tint manually here.
	5. Plotting begins from line 76.
	
# Codes for Figure 3

1. **Figure 3A and B:** Run the previous code (__/Figure_2/SDE_different_Dt.m__) to generate plots for Distance between the expected and the derived Deterministic function for pairwise model (r3 = 0) and the ternary model (r3 = 0.08), as a function of a more continuous $\Delta t$.
2. **Figure 3C:** Run __/Figure_3/optDt_changing_r1.m__ to generate the data for the figure (r1 = s = flip(1./(2*(10:10:200))), r2 = 1, r3 = 0.08, N = size = 50) and use plot or scatter command in matlab to plot optimal Dt vs Correlation time.
3. **Figure 3D:** Run __/Figure_3optDt_changingN.m__ to generate the data for the figure (r1 = 0.01, r2 = 1, r3 = 0.08, N = size = 50:1:200) and use plot or scatter command in matlab to plot optimal Dt vs Correlation time.

# Codes for Figure 4

1. **Figure 4A:** 
	1. Use code __/Figure_2/Gillespie_stochastic_process.m__ to generate three time series for the pairwise model (r1 = 0.01, r2 = 1, r3 = 0) for N = 50, 100 and 200.
	2. Run each of these time series on code __/Figure_2/SDE_different_Dt.m__ to calculate Distance between the expected and the derived Stochastic function as a function of $\delta t$.
2. **Figure 4B:**
	1. Use code __/Figure_4/varying_resolution.m__ to generate time series with different resolution. Set r1 = 0.01, r2 = 1, r3 = 0 for the pairwise model.
	2. Perform the above step for different N (50, 100, 200). The outputs will be stored in the folder __Characterizing_noise/pairwise/varying_resolution/N_15__, for example. 
3. **Figure 4C**
	1. Perform the steps to generate Figure 4A for the ternary models case (r1 = 0.01, r2 = 1, r3 = 0.08).
4. **Figure 4D**
	1. Perform the steps in Figure 4B for the ternary models case (r1 = 0.01, r2 = 1, r3 = 0.08).
5. **Figure 4E and 4F**
	1. Generate a time series using code __/Figure_2/Gillespie_stochastic_process.m__ with Tint = 10 and 50 for the pairwise model (r1 = 0.01, r2 = 1, r3 = 0) for N = 50.
	2. Run these time series on code __/Figure4/noise_analysis.m__ to get the plots 4E and 4F.
6. **Figure 4G and 4H**
	1. Perform the steps used to generate Figure 4E and 4F for the ternary model case ((r1 = 0.01, r2 = 1, r3 = 0.08)) for N = 50.
