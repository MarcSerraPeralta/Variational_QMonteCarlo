import numpy as np

import lib as lib
from inputs import *

#############################################

# HARMONIC OSCILLATOR
E_local_f = E_local_Harmonic_Oscillator
prob_density = prob_density_Harmonic_Oscillator
dim = 1 # dimension of configuration space
opt_args = {"method":"scan1D", "init_alpha":np.array([0.25]), "step":0.05, "final":np.array([0.75])}

# Monte Carlo integration params
N_steps = 30000
N_walkers = 400
N_skip = 4000
L_start = 5
N_cores = -1 # if set to -1, it uses all cores available

# Results saving params
file_name = "output.csv"

#############################################

data_list = []
alpha_list = []
trial_move = None

optimizer = lib.Optimizer(opt_args)

while not optimizer.converged:
	data = lib.MC_integration(E_local_f, prob_density, optimizer.alpha, dim, N_steps=N_steps, N_walkers=N_walkers, 
								N_skip=N_skip, L_start=L_start, trial_move=trial_move, N_cores=N_cores)

	if abs(data[2] - 0.5) >= 0.05: #if acceptance ratio deviates from 0.5 by 0.05 we change the trial move
		trial_move = None
		print("The trial move will be updated")
	else:
		trial_move = data[3]

	alpha = optimizer.alpha

	data_list += [data]
	alpha_list += [alpha]
	

	# UPDATE ALPHA

	if opt_args["method"] == "scan1D":
		optimizer.update_alpha()
		print("E({:0.7f}) = {:0.15f} | var(E) = {:0.15f}".format(*alpha, *data))

	if opt_args["method"] == "steepest_descent1D":
		optimizer.update_alpha(data)
		print("E({:0.7f}) = {:0.15f} | var(E) = {:0.15f}".format(*alpha, *data))

	if opt_args["method"] == "steepest_descent_ANY_D":

		# calculate gradient
		if np.all(optimizer.alpha_old != None):
			E = np.zeros(optimizer.dim_alpha)
			for j in range(optimizer.dim_alpha):
				e_j = np.zeros(optimizer.dim_alpha)
				e_j[j] = 1
				alpha_j = optimizer.alpha_old+optimizer.step*e_j
				data_j = lib.MC_integration(E_local_f, prob_density, alpha_j, dim,
								N_steps=N_steps, N_walkers=N_walkers, N_skip=N_skip, L_start=L_start, trial_move=trial_move)
				E[j] = data_j[0]
		else:
			E = data[0]

		optimizer.update_alpha(E)
		print("E(", alpha,") = {:0.15f} | var(E) = {:0.15f}".format(data[0], data[1]))
		

lib.save(file_name, alpha_list, data_list)