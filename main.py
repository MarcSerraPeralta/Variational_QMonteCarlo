import numpy as np

import lib as lib

#############################################

H = ""
psi_t = ""
init_alpha = np.array([])

N_steps = 5000
N_walkers = 250
N_skip = 100
L_start = 1

opt_args = {"method":None, "init_alpha": init_alpha}

file_name = "output.csv"

#############################################

E_local_f = lib.get_E_local_f(H, psi_t)
prob_density = lib.get_prob_density(psi_t)

data_list = []
alpha_list = []

optimizer = lib.Optimizer(opt_args)

while not optimizer.converged:
	data = lib.MC_integration(E_local_f, prob_density, optimizer.alpha, N_steps=N_steps, N_walkers=N_walkers, N_skip=N_skip, L_start=L_start)
	alpha = optimizer.alpha

	data_list += [data]
	alpha_list += [alpha]

	optimizer.update_alpha(data)

lib.save(file_name, alpha_list, data_list)