Install HDfair using HDfair_1.0.0.tar.gz

1. lowD_dense.R is the code to reproduce the simulation experiments in the low-dimensional dense scenario.

2. sparse_scenario.R is the code to reproduce all other simulation experiments

  2.1 Modify the following settings in setting.R:
	M = 1 or 3 sources
	P = 100 or 1000 variables
	sd.th.m = 0.1, 0.2, 0.3 between-source heterogeneity
	sd.th.a = 0.1, 0.2 between group heterogeneity

  2.2 Run the sim() function in sparse_scenario.R, the only input is seed number
	Line 397-399 may need to be modified depending on where and how the function is executed.
	Download ARMUL from https://github.com/kw2934/ARMUL into the root folder

  2.3 Run summary.R to generate the result tables