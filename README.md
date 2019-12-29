# Characterizing_noise

## This repository contains codes and data used for analysis in the article "Noise-Induced Effects in Collective Dynamics and Inferring Interactions from Data".

# Codes for Figure 2

1. **Generate polarisation time-series for a given set of parameters:** A matlab code (/Figure_2/Gillespie_stochastic_process.m) can be used to compute the polarization/order parameter.
	1. The important parameters that need to be set are N (system size), r1, r2, r3, r4 (which are for different reaction rates).
	2. Other parameters include Tint (sampling time which is 50 in the paper), Tend (which is equal to Tint*number of iterations). Set number of iterations in this formula.
	3. Set parameter rel for number of realizations/replicates of the simulations.
	4. The output is a time series array (S) of size (number of iterations*rel), and an array (tSample) storing system time.
	5. Autocorrelation time of the time series (est_tau) can be calculated. exp_tau is the expected autocorrealtion time for the given reaction rate and system size parameters.
