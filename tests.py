import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import lib as lib
import inputs as inputs


def check_random_walker_1D(tm_sigma, prob_density=None, alpha=None, N_steps=50000):
	"""
	Check for lib.random_walker() function with 1D input using visual inspection.

	Parameters
	----------
	tm_sigma : float
		Initial guess of the initial trial move variance for lib.random_walker()
	prob_density : function or None
		Probability density distribution for random_walker()
		If None, it uses a Gaussian probability density distribution 
	alpha : np.ndarray
		Parameters for prob_density()
		If None, alpha=np.array([1])
	N_steps : int
		Number of steps of the random walker

	Returns
	-------
	None
	"""

	if prob_density is None: prob_density = Gaussian
	if alpha is None: alpha = np.array([1])

	x_points, acceptance_probability, acceptance_ratio = lib.random_walker(prob_density, alpha, N_steps, np.zeros(1), tm_sigma)

	plt.hist(x_points, bins=40, density=True)
	xmin, xmax = np.min(x_points), np.max(x_points)
	x = np.linspace(xmin, xmax, 1000)
	plt.plot(x, prob_density(x, alpha), "--")
	plt.xlim(xmin, xmax)
	plt.xlabel("r")
	plt.ylabel("Probability density function")
	plt.tight_layout()
	plt.show()

	x1 = np.linspace(1, N_steps, N_steps)
	print(x.shape, acceptance_probability.shape)
	plt.scatter(x1, acceptance_probability, s=1)
	plt.plot(x1, acceptance_ratio*np.ones(N_steps), 'r')
	plt.ylim(0, 1)
	plt.xlabel("step")
	plt.ylabel("Acceptance probability")
	plt.tight_layout()
	plt.show()

	
	graph=sns.jointplot(x=x1, y=acceptance_probability, s=1, hue_norm=(0,1))
	graph.ax_joint.axhline(y=acceptance_ratio, c='r')
	plt.show()

	return


def check_random_walker_3D(tm_sigma, prob_density=None, alpha=None, N_steps=50000):
	"""
	Check for lib.random_walker() function with 3D input using visual inspection.

	Parameters
	----------
	tm_sigma : float
		Initial guess of the initial trial move variance for lib.random_walker()
	prob_density : function or None
		Probability density distribution for random_walker()
		If None, it uses a Gaussian probability density distribution 
	alpha : np.ndarray
		Parameters for prob_density()
		If None, alpha=np.array([1])
	N_steps : int
		Number of steps of the random walker

	Returns
	-------
	None
	"""

	if prob_density is None: prob_density = lambda x,y,z,std: Gaussian(x,std)*Gaussian(y,std)*Gaussian(z,std)
	if alpha is None: alpha = np.array([1])

	x_points = lib.random_walker(prob_density, alpha, N_steps, np.zeros(3), tm_sigma)

	plt.hist(x_points[:,0], bins=40, density=True, alpha=0.3, label="r1")
	plt.hist(x_points[:,1], bins=40, density=True, alpha=0.3, label="r2")
	plt.hist(x_points[:,2], bins=40, density=True, alpha=0.3, label="r3")
	plt.xlabel("r_i")
	plt.ylabel("Probability density function")
	plt.legend()
	plt.tight_layout()
	plt.show()

	return


def check_random_walkers_1D(tm_sigma, prob_density=None, alpha=None, N_steps=50000):
	"""
	Check for lib.random_walker() function with 3D input using visual inspection.

	Parameters
	----------
	tm_sigma : float
		Initial guess of the initial trial move variance for lib.random_walker()
	prob_density : function or None
		Probability density distribution for random_walker()
		If None, it uses a Gaussian probability density distribution 
	alpha : np.ndarray
		Parameters for prob_density()
		If None, alpha=np.array([1])
	N_steps : int
		Number of steps of the random walker

	Returns
	-------
	None
	"""

	if prob_density is None: prob_density = Gaussian
	if alpha is None: alpha = np.array([1])

	x_points = lib.random_walkers(prob_density, alpha, N_steps, np.zeros((3,1)), tm_sigma)

	plt.hist(x_points[:,0,0], bins=40, density=True, alpha=0.3, label="RW 1")
	plt.hist(x_points[:,1,0], bins=40, density=True, alpha=0.3, label="RW 2")
	plt.hist(x_points[:,2,0], bins=40, density=True, alpha=0.3, label="RW 3")
	xmin, xmax = np.min(x_points), np.max(x_points)
	x = np.linspace(xmin, xmax, 1000)
	plt.plot(x, prob_density(x, alpha), "--")
	plt.xlim(xmin, xmax)
	plt.xlabel("r_i")
	plt.ylabel("Probability density function")
	plt.legend()
	plt.tight_layout()
	plt.show()

	return


def monitor_ar_hydrogen(trial_move, N_walkers=10000, N_steps=10000):
	"""
	Plots a histogram of the acceptance ratio for hydrogen for all the random walkers.

	Parameters
	----------
	trial_move : float
		Standard deviation that defines the trial move according to a normal distribution
	N_walkers : int
		Number of random walkers
	N_steps : int
		Number of steps that each random walker takes

	Returns
	-------
	None
	"""

	L_start = 5
	dim = 3
	alpha = np.array([1])

	prob_density = inputs.prob_density_Hydrogen_atom
	init_points = lib.rand_init_point(L_start, dim, N_walkers)

	_, _, acceptance_ratio = lib.random_walkers(prob_density, alpha, N_steps, init_points, trial_move)
	plt.hist(acceptance_ratio, bins=40, density=True)
	plt.xlabel("acceptance ratio")
	plt.ylabel("density of walkers")
	plt.show()

	return 


def Gaussian(x, std):
	"""
	Returns the probability of x in a 1D Gaussian distribution of
	zero mean and standard deviation of std. 

	Parameters
	----------
	x : float
		Variable of the Gaussian distribution
	std : float
		Standard deviation of the aussian distribution

	Returns
	-------
	f : float
		Probability of x in a 1D Gaussian distribution N(0, std)
	"""

	f = np.exp(-x**2/(2*std**2))/(std*np.sqrt(2*np.pi))

	return f


if __name__ == '__main__':
	print("CHECK acceptance ratio Hydrogen (tm_sigma=opt)...")
	opt_sigma = lib.find_optimal_trial_move(inputs.prob_density_Hydrogen_atom,np.array([1]),3,2)
	print("The optimal trial move is: ", opt_sigma)
	monitor_ar_hydrogen(opt_sigma)
	print("DONE")

	print("CHECK GAUSSIAN 1D (tm_sigma=opt)...")
	opt_sigma = lib.find_optimal_trial_move(Gaussian,np.array([1]),1,1)
	print("The optimal trial move is: ", opt_sigma)
	check_random_walker_1D(opt_sigma)
	print("DONE")

	print("CHECK X^2*GAUSSIAN 1D (tm_sigma=opt)...")
	f = lambda x, std: std*x**2*Gaussian(x, std)
	opt_sigma = lib.find_optimal_trial_move(f,np.array([1]),1,1)
	print("The optimal trial move is: ", opt_sigma)
	check_random_walker_1D(opt_sigma, prob_density=f)
	print("DONE")

	print("CHECK GAUSSIAN 1D (tm_sigma=1)...")
	check_random_walker_1D(1)
	print("DONE")

	print("CHECK X^2*GAUSSIAN 1D (tm_sigma=1)...")
	f = lambda x, std: std*x**2*Gaussian(x, std)
	check_random_walker_1D(1, prob_density=f)
	print("DONE")

	print("CHECK X^2*GAUSSIAN 1D (tm_sigma=0.1)...")
	f = lambda x, std: std*x**2*Gaussian(x, std)
	check_random_walker_1D(0.1, prob_density=f)
	print("DONE")

	print("CHECK X^2*GAUSSIAN 1D (tm_sigma=0.01)...")
	f = lambda x, std: std*x**2*Gaussian(x, std)
	check_random_walker_1D(0.01, prob_density=f)
	print("DONE")

	print("CHECK GAUSSIAN 3D (tm_sigma=1)...")
	check_random_walker_3D(1)
	print("DONE")

	print("CHECK GAUSSIAN 1D WITH MULTIPLE RANDOM WALKERS (tm_sigma=1)...")
	check_random_walkers_1D(1)
	print("DONE")