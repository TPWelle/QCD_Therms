import numpy as np
import therm

#instantiating models with their parameters
th_h = therm.EXI(1.44991e10)
th_q = therm.pQCD(11.57418, 18.02198)
th_SF = therm.Standard(198.31430, 4)
#Crossover takes two equations of state and a switching function
th = therm.Crossover(th_h, th_q, th_SF)

#arrays of temperatures and chemical potentials
Ts = np.linspace(10.0,800.0,200)
mus = 0.0 * Ts

#include these arrays and compute the EOS at each point
th.set_Tmu(Ts, mus)
th.compute()

#results
Pressures = th.Y
EntropyDensities = th.Y_T
BaryonDensities = th.Y_m
Temperatures = th.Temps #same as Ts
ChemicalPotentials = th.mu_Bs #same as mu
