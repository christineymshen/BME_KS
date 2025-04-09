This repo contains codes for simulation study results for the manuscript "Asymptotic well-calibration of the posterior predictive $p$-value under the modified Kolmogorov-Smirnov test".

These codes were run on a computing cluster using slurm array jobs. There are a few differences to the codes I used:

1. I have removed the folder names in all R scripts.
2. For `rstan` runs,  the codes use a file named "gamma_mod.rds". This file can be prepared using `stan_model("gamma.stan")`. As the compiler runs slightly differently in different environment, I did not upload the compiled file used on the server.

The codes are organized as follows:
1. run number 1-4 are for Figure 1 which is simulation results on the level of the $p$-values. Each of these run have the same setup except for using different sample sizes.
2. run number 5-8 are for Figure 2 which is simulation results on the power of the $p$-values. How each run corresponds to each alternative model can be found out from the "summary.R" file.
3. the `functions.R` file contains all the helper functions for simulation runs, aggregating results, and for plotting.
4. the `summary.R` file contains codes to aggregate run results as well as plotting.
