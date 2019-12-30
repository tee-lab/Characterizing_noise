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
	1. Run one column of the time series array (for e.g., S(:,1)) at a time with this code.
	2. In line 8: Set different time scales over which you want to calculate the deterministic and the stochastic part.
	3. Uncomment line 33 and set the exact parameters used to generate the time series using code __/Figure_2/Gillespie_stochastic_process.m__
	4. tSample from the code __/Figure_2/Gillespie_stochastic_process.m__ (or from real data) is important for calculaing Tint in line 35.
	5. Plotting begins from line 76.
