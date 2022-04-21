import numpy as np

################################
#     EXAMPLE OF FUNCTIONS     #
################################

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


################################
#     HARMONIC OSCILLATOR      #
################################

def E_local_Harmonic_Oscillator(r, alpha):
	"""
	Returns the value of the local energy for the 1D Harmonic Oscillator
	at the specified positions r and parameters alpha.

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

	E_local = alpha + (0.5 - 2*alpha**2)*r**2

	return E_local


def prob_density_Harmonic_Oscillator(r, alpha):
	"""
	Returns the value of the probability density function for the 1D Harmonic Oscillator
	at the specified positions r and parameters alpha.

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

	prob = np.exp(-2*alpha*r**2)

	return prob