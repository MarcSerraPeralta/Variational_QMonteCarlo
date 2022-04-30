import numpy as np
import matplotlib.pyplot as plt

import lib as lib


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

	x = np.linspace(1, N_steps, N_steps)
	print(x.shape, acceptance_probability.shape)
	plt.scatter(x,acceptance_probability, s=1)
	plt.plot(x,acceptance_ratio*np.ones(N_steps),'r')
	plt.ylim(0, 2)
	plt.xlabel("step")
	plt.ylabel("Acceptance probability")
	plt.tight_layout()
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
	"""
	return np.exp(-x**2/(2*std**2))/(std*np.sqrt(2*np.pi))


if __name__ == '__main__':
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