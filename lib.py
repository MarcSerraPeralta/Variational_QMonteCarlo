from sympy import simplify, lambdify, conjugate, integrate, oo, DiracDelta
from multiprocessing import Pool
import numpy as np
from scipy.optimize import brentq 

#######################################
#	   PROCESS INPUT PARAMETERS
#######################################

def get_E_local_f(H, psi_t, var):
	"""
	Returns the function local energy given a Hamiltonian an trial wave function.

	Parameters
	----------
	H : sympy expression
		Hamiltonian of the system
	psi_t : sympy expression
		Trial wavefunction
	var : sympy symbols
		Variables for position r and for parameters alpha

	Returns
	-------
	E_local_f : function(r, alpha)
		Local energy function depending on r and alpha
	"""

	E_local_f = H/psi_t
	E_local_f = simplify(E_local_f)
	E_local_f = lambdify(var, E_local_f, modules=['numpy'])

	return E_local_f


def get_prob_density(psi_t, var):
	"""
	Returns the probability density function for metropolis algorithm given a trial wave function.

	Parameters
	----------
	psi_t : sympy expression
		Trial wavefunction
	var : sympy symbols
		Variables for position r and for parameters alpha

	Returns
	-------
	prob_density : function(r, alpha)
		Probability density function for the specified trial wave function
	"""

	prob_density = conjugate(psi_t)*psi_t
	prob_density = simplify(prob_density)
	prob_density = lambdify(var, prob_density, modules=['numpy'])

	return prob_density


#######################################
#			   EXAMPLES
#######################################

def prob_density(r, alpha):
	"""
	Returns the value of the probability density function at the specified positions r and parameters alpha.

	dim = dimension of the integral space.

	Parameters
	----------
	r : np.ndarray(N*dim)
		Positions of N particles in the form: r1x, r1y, ..., r2x, r2y, ... 
	alpha : np.ndarray	
		Parameters of the trial wave function

	Returns
	-------
	prob : np.ndarray(N)
		Probability density function at the specified r and alpha
	"""

	prob = 0

	return prob


def E_local_f(r, alpha):
	"""
	Returns the value of the local energy at the specified positions r and parameters alpha.

	dim = dimension of the integral space.

	Parameters
	----------
	r : np.ndarray(N*dim)
		Positions of N particles in the form: r1x, r1y, ..., r2x, r2y, ... 
	alpha : np.ndarray
		Parameters of the trial wave function

	Returns
	-------
	E_local : float
		Local energy at the specified r and alpha
	"""

	E_local = 0

	return E_local


#######################################
#	    MONTE CARLO INTEGRATION
#######################################

def random_walker(prob_density, alpha, N_steps, init_point, trial_move):
	"""
	Returns steps of a random walker that follows a Markov chain given a probability
	density function using the Metropolis algorithm. 
	'dim' is obtained from init_point. 

	Parameters
	----------
	prob_density : function(r, alpha)
		Probability density function depending on position r and parameters alpha
	alpha : np.ndarray
		Parameters of the trial wave functon
	N_steps : int
		Number of steps that the random walker takes
	init_point : np.ndarray(dim)
		Starting point of the random walker
	trial_move : float
		Standard deviation that defines the trial move according to a normal distribution

	Returns
	-------
	steps : np.ndarray(N_steps, dim)
		Steps of the random walker
	"""

	dim = init_point.shape[0]
	steps = np.zeros((N_steps,dim))
	steps[0] = init_point

	for i in np.arange(1, N_steps):
		print("{}/{}\r".format(i+1, N_steps), end="")
		next_point = steps[i-1] + np.random.normal(0, trial_move, size=dim)
		if (np.random.rand(1) <= prob_density(*next_point, *alpha)/prob_density(*steps[i-1], *alpha)).all():
			steps[i] = next_point
		else:
			steps[i] = steps[i-1]
	
	return steps


def random_walkers(prob_density, alpha, N_steps, init_points, trial_move):
	"""
	Returns steps of N_walkers random walkers that follows a Markov chain given a probability
	density function using the Metropolis algorithm. 
	'N_walkers' and 'dim' are obtained from init_points.

	Parameters
	----------
	prob_density : function(r, alpha)
		Probability density function depending on position r and parameters alpha
	alpha : np.ndarray
		Parameters of the trial wave functon
	N_steps : int
		Number of steps that each random walker takes
	init_point : np.ndarray(N_walkers, dim)
		Starting points of all random walkers
	trial_move : float
		Standard deviation that defines the trial move according to a normal distribution

	Returns
	-------
	steps : np.ndarray(N_steps, N_walkers, dim)
		Steps of the random walker
	"""

	N_walkers, dim = init_points.shape
	steps = np.zeros((N_steps, N_walkers, dim))
	steps[0] = init_points

	for i in np.arange(1, N_steps):
		print("{}/{}\r".format(i+1, N_steps), end="")
		next_point = steps[i-1] + np.random.normal(0, trial_move, size=N_walkers*dim).reshape(N_walkers, dim)

		to_change = np.where(np.random.rand(N_walkers) <= prob_density(*next_point.T, *alpha)/prob_density(*steps[i-1].T, *alpha))

		steps[i] = steps[i-1]
		steps[i, to_change] = next_point[to_change]

	return steps


def rand_init_point(system_size, dim, N_points):
	"""
	Returns a random initial point for the random walkers given a typical size of the system according to a normal distribution. 

	Parameters
	----------
	system_size : float
		Typical size of the system, e.g. for the Hydrogen atom it could be 2*a_0
	dim : int
		Dimension of the configuration space, i.e. number of degrees of freedom in the system

	Returns
	-------
	init_point : np.ndarray(N_points, dim)
		Random initial point for the random walkers
	"""

	init_point = np.random.normal(scale=system_size, size=(N_points, dim))

	return init_point


def dev_av_rate(trial_move, prob_density, alpha, dim, N_av=100):
	"""
	Returns the deviation of the average acceptance ratio for a random walker giving N_av steps from 0.5

	Parameters
	----------
	prob_density : function(r, alpha)
		Probability density function depending on position r and parameters alpha
	alpha : np.ndarray
		Parameters of the trial wave function
	dim : int
		Dimension of the configuration space, i.e. number of degrees of freedom in the system
	trial_move : float
		Current initial trial move variance
	N_av : int
		Number of steps the walker takes to compute the acceptance ratio average
	Returns
	-------
	dev_av_ratio : float
		Average acceptance ratio for a random walker giving N_av steps
	"""

	av_ratio = 0
	steps = np.zeros((N_av, dim))
	steps[0] = np.zeros(dim)

	for i in np.arange(1, N_av):
		next_point = steps[i-1] + np.random.normal(scale=trial_move, size=(dim))
		ratio = min(prob_density(*next_point, *alpha)/prob_density(*steps[i-1], *alpha),1)
		if np.random.rand(1) <= ratio:
			steps[i] = next_point
		else:
			steps[i] = steps[i-1]
		av_ratio += ratio

	dev_av_ratio = av_ratio/N_av - 0.5

	return dev_av_ratio


def find_optimal_trial_move(prob_density, alpha, dim, trial_move_init, maxiter=5000, N_av=2000, tol=0.05):
	"""
	Returns trial_move such that the corresponding average accepting ratio is 0.5+-tol. 

	Parameters
	----------
	prob_density : function(r, alpha)
		Probability density function depending on position r and parameters alpha
	alpha : np.ndarray
		Parameters of the trial wave function
	dim : int
		Dimension of the configuration space, i.e. number of degrees of freedom in the system
	trial_move_init : float
		Initial guess of the initial trial move variance
	maxiter : int
		Maximum number of iterations 
	tol : float
		Tolerance for the deviation of the average acceptance ratio from 0.5

	Returns
	-------
	opt_trial_move: float
		trial_move such that the corresponding average accepting ratio is 0.5 +- tol
	"""

	arguments = (prob_density, alpha, dim, N_av)

	# Find a zero in dev_av_rate between 10*trial_move_init and trial_move_init/100,
	# the function must be of oposite signs at the two points
	opt_trial_move = brentq(dev_av_rate, trial_move_init/1000, 10*trial_move_init, args=arguments, maxiter=maxiter, xtol=tol) 
												
	return opt_trial_move


def MC_integration(E_local_f, prob_density, alpha, dim, N_steps=5000, N_walkers=250, N_skip=0, L_start=1):
	"""
	Returns expectation value of the energy E(alpha) averaged over N_walkers random walkers
	using Monte Carlo integration. 

	Parameters
	----------
	E_local_f : function(r, alpha)
		Local energy function depending on r and alpha
	prob_density : function(r, alpha)
		Probability density function depending on position r and parameters alpha
	alpha : np.ndarray
		Parameters of the trial wave function
	N_steps : int
		Number of steps that the random walker takes
	N_walkers : int
		Number of random walkers
	N_skip : int
		Number of initial steps to skip for the integration
	L_start : float
		Length of the box in which the random walkers are initialized randomly

	Returns
	E_alpha : float
		Expectation value of the energy for given parameters of the trial wave function
	E_alpha_std : float
		Standard deviation of E_alpha computed from E_alpha_walkers (E_alpha for each walker)
	-------
	...
	"""

	init_points = rand_init_point(L_start, dim, N_walkers)
	trial_move = find_optimal_trial_move(prob_density, alpha, dim, L_start) # I don't know if trial_move_init should be L_start
	print("Optimal trial_move is", trial_move)

	# initialization of variables and prepare the inputs
	inputs = [(prob_density, alpha, N_steps, dim, init_points[i], trial_move) for i in range(N_walkers)]

	# multiprocessing
	#pool = Pool() # uses maximum number of processors available
	#data_outputs = pool.map(random_walker, inputs)
	data_outputs = random_walkers(prob_density, alpha, N_steps, init_points, trial_move)

	# do stuff with data_outputs
	total_steps = np.array(data_outputs)[N_skip:, :, :]
	E_alpha, E_alpha_std = MC_sum(E_local_f, total_steps, alpha)

	return E_alpha, E_alpha_std


def MC_sum(E_local_f, steps, alpha):
	"""Computes expectation value of energy, E(alpha), given a distribution
	of steps

	Parameters
	----------
	E_local_f : function(r, alpha)
		Local energy function depending on r and alpha
	steps : np.ndarray(N_steps, N_walkers, dim)
		Array of points to be used in the computation of the integral.
		N_steps is the number of steps each walker takes.
		N_walkers is the number of walkers.
		dim is the dimension of the integral space.
	alpha : np.ndarray
		Parameters of the trial wave function

	Returns
	-------
	E_alpha : float
		Expectation value of the energy for given parameters of the trial wave function
	E_alpha_std : float
		Standard deviation of E_alpha computed from E_alpha_walkers (E_alpha for each walker)
	"""

	N_steps = steps.shape[0]

	E_local = E_local_f(*steps.T, *alpha)
	E_local_walkers = np.average(E_local, axis=1)
	E_alpha = np.average(E_local_walkers)
	E_alpha_std = np.std(E_local_walkers)

	return E_alpha, E_alpha_std


#######################################
#	     NUMERICAL OPTIMIZATION
#######################################

class Optimizer:
	"""
	Wrapper for numerical optimization of alpha.

	Initial Parameters
	------------------
	method : str
		Method for optimizing alpha
		Options: "scan", "steepest_descent"
	init_alpha : np.ndarray
		Initial values of alpha

	Functions
	---------
	update_alpha
		Updates alpha using the specified method and checks if the 
		optimization has converged (stored in Optimizer.converged)
	"""

	def __init__(self, opt_args):
		self.method = opt_args["method"]
		self.converged = False
		self.alpha = opt_args["init_alpha"]

		if self.method == "scan1D":
			self.step = opt_args["step"]
			self.final = opt_args["final"]
			self.alpha_range = np.arange(self.alpha[0], self.final + self.step, self.step).tolist()
			self.alpha = np.array([self.alpha_range.pop(0)]) # deletes also first element from list
		return

	def update_alpha(self, args):
		if self.method == "scan1D":
			if len(self.alpha_range) != 0:
				self.alpha = np.array([self.alpha_range.pop(0)]) # deletes also first element from list
			else:
				self.converged = True

		elif self.method == "steepest_descent":
			self.alpha = steepest_descent(self.alpha, self.args)
			if True:
				self.converged = True

		else: 
			self.converged = True
		
		return


def steepest_descent(alpha_old, args):
	"""
	Returns the new value of alpha using the method of steepest descent.

	Parameters
	----------
	alpha_new : np.ndarray
		Old value of alpha
	args : dict
		Information about the steepest descent method

	Returns
	-------
	alpha_new : np.ndarray
		New value of alpha
	"""

	alpha_new = 0

	return alpha_new


#######################################
#			  SAVE RESULTS
#######################################

def save(file_name, alpha_list, data_list, alpha_labels=None, data_labels=["E", "var(E)"]):
	"""
	Saves alpha and data results to file_name.

	Parameters
	----------
	file_name : str
		Name of the file to store the results
	alpha_list : list of np.ndarray
		List of alpha values
	alpha_list : list of list
		List of data values
	alpha_labels : list of str
		Labels for the elements of alpha_list[i]
		If None, it assumes labels: alpha1, alpha2, alpha3, ...
	data_labels : list of str
		Labels for the elements of data_list[i]
		If None, it assumes labels: data1, data2, data3, ...

	Returns
	-------
	None
	"""

	f = open(file_name, "w")
	header = ""
	if alpha_labels is None:
		header += ",".join(["alpha{}".format(i+1) for i in range(len(alpha_list[0]))]) + ","
	else: 
		header += ",".join(alpha_labels) + ","
	if data_labels is None:
		header += ",".join(["data{}".format(i+1) for i in range(len(data_list[0]))]) + ","
	else: 
		header += ",".join(data_labels) + ","
	header += "\n"
	f.write(header)

	for k in range(len(alpha_list)):
		s = ("{},"*len(alpha_list[k])).format(*alpha_list[k])
		s += ("{},"*len(data_list[k])).format(*data_list[k])
		s += "\n"
		f.write(s)

	f.close()

	return 