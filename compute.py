
#=============================================================================
# The following contains examples to show how to setup the various models
#=============================================================================

import numpy as np
import therm

#=============================================================================
#First, instantiate a hadronic and a quark-gluon model and set their parameters
#=============================================================================

#====Hadronic models====

# th_h = therm.PT() #PT model
th_h = therm.EXI(1.44991e10) #EXI model (takes eps_0 as a parameter)
# th_h = therm.EXII(1.06307e10) #EXII model (takes eps_0 as a parameter)

# --- RMF model ---
# th_h = therm.RMF()

# set the couplings

#sigma, omega, and rho couplings to nucleons
#other coupling determined from these by symmetries
# th_h.g_sig, th_h.g_om, th_h.g_rho = 7.932718, 8.608836, 8.615788

#2nd and 3rd order sigma couplings
# th_h.b, th_h.c = 0.008872, -0.001628


# - which baryons/mesons to use can be changed using the following - :

# th_h.N_bar, th_h.N_mes = 8, 4

# these are the default values corresponding to N,Lambda,Sigma,Xi in order of decreasing charge
# th_h.m_baryons=[938.27,939.57,1115.7,1189.4,1192.6,1197.4,1314.9,1321.7];

# these are the default values corresponding to sigma, omega, phi, rho
# th_h.m_mesons=[500.0, 783.0, 1019.5, 775.5];

# all baryons are assumed to be spin 1/2, while the first meson species is assumed to be
# a scalar, and all others vectors.

# - some things computed in the RMF model -
# th_h.Meson_RMFs  (Values for the Meson mean fields for each T and mu)
# th_h.ns_Bar (density of each baryon species for each T and mu)
# th_h.ns_scl_Bar (scalar density of each baryon species for each T and mu)



#====Quark-gluon models====

th_q = therm.pQCD(3.684, 4.345) #perturbative QCD (takes C_E and C_S)
#Taken to order g^6 by default. can be changed by changing th_q.K

# th_q = therm.ColdQCD() # Low-T expansion for perturbative QCD

#=============================================================================
#Next, if you are using a Crossover model, you need to set the switching function.
#Then collect these into the Crossover model
#=============================================================================

th_SF = therm.Standard(198.31430, 4) #Standard switching function (params are T0 and r)

#Continuous switching function  with critical point(params are T0, Tc, muc, A, and r)
# th_SF = therm.Continuous(198.31430, 130.0, 450.0, .8, 4)

#mu_0 set to 3*pi*T0 by default in both of these
#it can be changed using th_SF.mu0


# Crossover model
th_CO = therm.Crossover(th_h, th_q, th_SF)

#=============================================================================
# If you are using the Embedded or Scaling models, these will modify the Standard function
# to include a critical point. These are set up as follows
#=============================================================================

# --- Embedded ---
th_emb = therm.Embedded(th_CO, 130.0, 450.0) #(params are th_CO, Tc, muc)

# These are the critical exponents. By default, the are set to the 3d-Ising values
# th_emb.alpha, th_emb.beta, th_emb.gamma, th_emb.delta

# the exponent in the function r_x(T,mu) is 4 by default.
# this can be changed using th_emb.m_exponent. (should be even)

#set parameters a0, Ta, Td
th_emb.a0, th_emb.Ta, th_emb.Td = .2, 80.0, 200.0

#set parameters d_- and d_+
th_emb.d_m = 50.0
th_emb.d_p = th_emb.d_m/3.0


# - some things computed in the Embedded model -

# th_emb.nc #critical baryon density
# th_emb.Pc #critical Pressure

# These are dictionaries containing the values of mu_x(T) and its first 2
# derivatives defined at each value of T which is computed, i.e. calling th_emb.mux[100.0]
# gives mu_x(100.0 MeV) if that is one the temperatures in th_emb.Temps.
# th_emb.mux, th_emb.dmux, th_emb.ddmux

# this is an array containing the function R at each T, mu and also its first two derivatives
# with respect to T and mu. This is given in the order {R,R_T,R_m,R_TT,R_Tm,R_mm}
# th_emb.Rs


# --- Scaling ---

th_scl=therm.Scaling(th_CO,120.0,750.0)  #(params are th_CO, Tc, muc)

# These are the critical exponents. By default, the are set to the 3d-Ising values
# th_scl.alpha, th_scl.beta, th_scl.gamma, th_scl.delta

# takes the ratios n0/n_c and P0/Pc and from them sets n0, nc, P0, and Pc
th_scl.set_n0P0(0.1, 0.05)

# set other parameters
th_scl.h0, th_scl.m0, th_scl.c_star=.2, .5, .7

# sets the values of coefficients h_3 and h_5, and from them
# determines the coefficients g0, g1, g2, g3 from g(theta)
th_scl.set_hs(-.762, .008)

# - some things computed in the Scaling model -

# th_scl.theta0 #the maximum value of theta

# Gives the functions theta, and R at each value of T and mu_B
# th_scl.thetas, th_scl.Rs

# These are dictionaries containing the values of mu_x(T) and its first 2
# derivatives defined at each value of T which is computed, i.e. calling th_emb.mux[100.0]
# gives mu_x(100.0 MeV) if that is one the temperatures in th_scl.Temps.
# th_scl.mux, th_scl.dmux, th_scl.ddmux


# The window function and its T and mu derivative
# th_scl.Wfs, th_scl.Wf_Ts, th_scl.Wf_ms

#=============================================================================
# Enter the temperatures and chemical potentials you want to compute at
#=============================================================================

# You need two arrays with the same shape, listing all points where you want to
# compute thermodynamic variables.

# -- For 100 temperatures between 10 MeV and 800 MeV at mu_B = 0 MeV --
Ts = np.linspace(10.0,800.0,100)
mus = 0.0*Ts

# -- For 100 mus between 0 MeV and 1000 MeV at T = 100 MeV --
# mus = np.linspace(10.0,1000.0,100)
# Ts = 0.0*mus + 100.0

# -- For a grid of 50 by 50 Ts and mu_Bs --
# Ts = np.linspace(100.0,800.0,50)
# mus = np.linspace(0.0,1500.0,50)
# Ts, mus = np.meshgrid(Ts, mus)


th = th_CO
# th = th_emb
# th = th_scl

#Include these arrays (also initializes the output variables among other things)
th.set_Tmu(Ts, mus)

#=============================================================================
# Finally, perform the computation on get the results
#=============================================================================

# Perform the calculations
th.compute()

# === results ===

Pressures = th.Y
EntropyDensities = th.Y_T
BaryonDensities = th.Y_m
Temperatures = th.Temps #same as Ts
ChemicalPotentials = th.mu_Bs #same as mu



# Individual models can be obtained from composite models
# in the crossover models, this looks like this
th_1 = th.pP1 #the hadronic model
th_2 = th.pP2 #the quark-gluon model
th_s = th.pSF #the switching function

# Information can then be obtained from each model individually

P_hadr = th_1.Y   # the partial pressure of the hadronic phase
n_qgp  = th_2.Y_m # baryon density of the quark phase
SF     = th_s.Y   # The switching function
SF_T   = th_s.Y_T # The T-derivative of the switching function
T_0    = th_s.T0  # The T0 parameter from the switching function model
eps_0  = th_1.eps0# The eps0 parameter from the hadronic model

#in the Embedded/Scaling model, the background model can be obtained via
# th_BG = th_emb.pP
# th_BG = th_scl.pP
