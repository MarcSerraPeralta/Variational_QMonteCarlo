import numpy as np
from sympy import symbols, diff, exp, pi, Heaviside, sqrt

import lib as lib
from inputs import *

#############################################

# HARMONIC OSCILLATOR
E_local_f = E_local_Harmonic_Oscillator
prob_density = prob_density_Harmonic_Oscillator
dim = 1 # dimension of configuration space
opt_args = {"method":"scan1D", "init_alpha":np.array([0.25]), "step":0.05, "final":np.array([0.75])}

# Monte Carlo integration params
N_steps = 25000
N_walkers = 250
N_skip = 1000
L_start = 5

# Results saving params
file_name = "output.csv"

#############################################

data_list = []
alpha_list = []

optimizer = lib.Optimizer(opt_args)

while not optimizer.converged:
	data = lib.MC_integration(E_local_f, prob_density, optimizer.alpha, dim,
							N_steps=N_steps, N_walkers=N_walkers, N_skip=N_skip, L_start=L_start)
	alpha = optimizer.alpha

	data_list += [data]
	alpha_list += [alpha]

	optimizer.update_alpha(data)

	if opt_args["method"] == "scan1D":
		print("E({:0.5f}) = {:0.15f} | var(E) = {:0.15f}".format(*alpha, *data))


lib.save(file_name, alpha_list, data_list)