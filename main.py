import numpy as np
from sympy import symbols, diff, exp, pi

import lib as lib

#############################################

# Hamiltonian and trial wavefunction
x, alpha = symbols('x alpha', real=True)
dim = 1 #dimension of configuration space
psi_t = (2*alpha/pi)**0.25 * exp(-alpha*x**2) 
normalized = True # if psi_t is normalized
H = -0.5 * diff(psi_t, x, 2) + 0.5 * x**2 * psi_t # Harmonic oscillator: 0.5*d^2 psi_t / dx^2 + 0.5*x^2*psi_t^2
var = [x, alpha] 

# Monte Carlo integration params
N_steps = 5000
N_walkers = 250
N_skip = 100
L_start = 1

# Alpha optimization params
init_alpha = np.array([])
opt_args = {"method":None, "init_alpha": init_alpha}

# Results saving params
file_name = "output.csv"

#############################################

E_local_f = lib.get_E_local_f(H, psi_t, var)
prob_density = lib.get_prob_density(psi_t, var)

data_list = []
alpha_list = []

optimizer = lib.Optimizer(opt_args)

while not optimizer.converged:
	data = lib.MC_integration(E_local_f, prob_density, optimizer.alpha, dim,
							N_steps=N_steps, N_walkers=N_walkers, N_skip=N_skip, L_start=L_start, normalized=normalized)
	alpha = optimizer.alpha

	data_list += [data]
	alpha_list += [alpha]

	optimizer.update_alpha(data)

lib.save(file_name, alpha_list, data_list)