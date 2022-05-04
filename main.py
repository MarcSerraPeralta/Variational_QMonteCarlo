import numpy as np

import lib as lib
from inputs import *

#############################################



# HARMONIC OSCILLATOR
E_local_f = E_local_Harmonic_Oscillator
prob_density = prob_density_Harmonic_Oscillator
dim = 1 # dimension of configuration space
opt_args = {"method":"scan1D", "init_alpha":np.array([0.25]), "step":0.05, "final":np.array([0.75])}

# HYDROGEN ATOM
E_local_f = E_local_Hydrogen_atom
prob_density = prob_density_Hydrogen_atom
dim = 3 # dimension of configuration space
opt_args = {"method":"scan1D", "init_alpha":np.array([0.75]), "step":0.05, "final":np.array([1.25])}
opt_args = {"method":"steepest_descent1D", "init_alpha":np.array([0.75]), "init_step":0.05, "gamma":0.25, "precision":1E-3}
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([0.75]), "init_step":np.array([0.05]), "gamma":0.25, "precision":1E-3}

# HELIUM ATOM
E_local_f = E_local_Helium_atom_GS_1param_numeric
E_local_f = E_local_Helium_atom_GS_1param_analytic
prob_density = prob_density_Helium_atom_GS_1param
dim = 6
opt_args = {"method":"steepest_descent1D", "init_alpha":np.array([0.15]), "init_step":0.01, "gamma":0.03, "precision":1E-5}
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([0.15]), "init_step":np.array([0.01]), "gamma":0.03, "precision":1E-5}
opt_args = {"method":"scan1D", "init_alpha":np.array([0.15]), "step":0.005, "final":np.array([0.17])}

# HELIUM ATOM GROUND STATE 2 PARAMETERS
E_local_f = E_local_Helium_atom_GS
prob_density = prob_density_Helium_atom_GS
dim = 6
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([2,0.5]), "init_step":np.array([0.1,-0.1]), "gamma":0.5, "precision":1E-5}

# HELIUM ATOM FIRST EXCITED 3 PARAMETERS
E_local_f = E_local_Helium_atom_1E
prob_density = prob_density_Helium_atom_1E
dim = 6
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([2,2,0.5]), "init_step":np.array([-0.1,-0.1,-0.1]), "gamma":0.5, "precision":1E-5}

# # HELIUM ATOM FIRST EXCITED 3 PARAMETERS
# E_local_f = E_local_Helium_atom_2E
# prob_density = prob_density_Helium_atom_2E
# dim = 6
# opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([2,2,0.5]), "init_step":np.array([-0.1,-0.1,-0.1]), "gamma":0.5, "precision":1E-5}

# # HELIUM ATOM FIRST EXCITED 3 PARAMETERS
# E_local_f = E_local_Helium_atom_3E
# prob_density = prob_density_Helium_atom_3E
# dim = 6
# opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([2,2,0.5]), "init_step":np.array([-0.1,-0.1,-0.1]), "gamma":0.5, "precision":1E-5}

# # HELIUM ATOM FIRST EXCITED 3 PARAMETERS
# E_local_f = E_local_Helium_atom_4E
# prob_density = prob_density_Helium_atom_4E
# dim = 6
# opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([2,2,0.5]), "init_step":np.array([-0.1,-0.1,-0.1]), "gamma":0.5, "precision":1E-5}

#############################################

# Monte Carlo integration params
N_steps = 30000
N_walkers = 400
N_skip = 4000
L_start = 5

# Results saving params
file_name = "output.csv"

#############################################

data_list = []
alpha_list = []
trial_move = None

optimizer = lib.Optimizer(opt_args)

while not optimizer.converged:
	data = lib.MC_integration(E_local_f, prob_density, optimizer.alpha, dim,
							N_steps=N_steps, N_walkers=N_walkers, N_skip=N_skip, L_start=L_start, trial_move=trial_move)

	if abs(data[2] - 0.5) >= 0.05: #if acceptance ratio deviates from 0.5 by 0.05 we change the trial move
		trial_move = None
		print("The trial move will be updated")
	else:
		trial_move = data[3]

	alpha = optimizer.alpha

	data_list += [data]
	alpha_list += [alpha]
	
	if opt_args["method"] == "scan1D":
		optimizer.update_alpha()
		print("E({:0.7f}) = {:0.15f} | var(E) = {:0.15f}".format(*alpha, *data))

	if opt_args["method"] == "steepest_descent1D":
		optimizer.update_alpha(data)
		print("E({:0.7f}) = {:0.15f} | var(E) = {:0.15f}".format(*alpha, *data))

	if opt_args["method"] == "steepest_descent_ANY_D":
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