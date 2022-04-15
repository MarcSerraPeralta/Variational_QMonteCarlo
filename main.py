import numpy as np
from sympy import symbols, diff, exp, pi, Heaviside, sqrt

import lib as lib

#############################################

# Hamiltonian and trial wavefunction
	# Harmonic oscillator: 0.5*d^2 psi_t / dx^2 + 0.5*x^2*psi_t^2
x, alpha = symbols('x alpha', real=True)
psi_t = exp(-alpha*x**2) 
H_psi_t = -0.5 * diff(psi_t, x, 2) + 0.5 * x**2 * psi_t 
var = [x, alpha] 
dim = 1 # dimension of configuration space
init_alpha, final_alpha = np.array([0.25]), np.array([0.75])

	# Hydrogen atom (unoptimized)
x, y, z, alpha = symbols('x y z alpha', real=True)
psi_t = exp(-alpha*sqrt(x**2 + y**2 + z**2)) 
H_psi_t = -0.5 * (diff(psi_t, x, 2) + diff(psi_t, y, 2) + diff(psi_t, z, 2)) - 1/sqrt(x**2 + y**2 + z**2) * psi_t 
var = [x, y, z, alpha] 
dim = 3 # dimension of configuration space
init_alpha, final_alpha = np.array([0.75]), np.array([1.25])

# Monte Carlo integration params
N_steps = 25000
N_walkers = 250
N_skip = 1000
L_start = 5

# Alpha optimization params
opt_args = {"method":"scan1D", "init_alpha":init_alpha, "step":0.05, "final":final_alpha}

# Results saving params
file_name = "output.csv"

#############################################

E_local_f = lib.get_E_local_f(H_psi_t, psi_t, var)
prob_density = lib.get_prob_density(psi_t, var)

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