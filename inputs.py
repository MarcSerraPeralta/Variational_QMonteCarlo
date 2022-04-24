import numpy as np

################################
#     EXAMPLE OF FUNCTIONS     #
################################

def prob_density(r1x, r1y, r2x, r2y, alpha1, alpha2):
	"""
	Returns the value of the probability density function at the specified positions r_ij and parameters alpha_k.

	Parameters
	----------
	r_ij : np.ndarray(N)
		Positions of N particles
	alpha_i : np.ndarray	
		Parameters of the trial wave function

	Returns
	-------
	prob : np.ndarray(N)
		Probability density function at the specified r and alpha
	"""

	prob = 0

	return prob


def E_local_f(r1x, r1y, r2x, r2y, alpha1, alpha2):
	"""
	Returns the value of the local energy at the specified positions r_ij and parameters alpha_k.

	Parameters
	----------
	r_ij : np.ndarray(N)
		Positions of N particles
	alpha_i : np.ndarray	
		Parameters of the trial wave function

	Returns
	-------
	E_local : np.ndarray(N)
		Local energy at the specified r and alpha
	"""

	E_local = 0

	return E_local


################################
#     HARMONIC OSCILLATOR      #
################################

# TRIAL WAVE FUNCTION = e^(-alpha*x^2)

def E_local_Harmonic_Oscillator(x, alpha):
	"""
	Returns the value of the local energy for the 1D Harmonic Oscillator
	at the specified positions r and parameters alpha.

	Parameters
	----------
	x : np.ndarray(N)
		Position of N particles
	alpha : np.ndarray
		Parameter of the trial wave function

	Returns
	-------
	E_local : np.ndarray(N)
		Local energy at the specified x and alpha
	"""

	E_local = alpha + (0.5 - 2*alpha**2)*x**2

	return E_local

def prob_density_Harmonic_Oscillator(x, alpha):
	"""
	Returns the value of the probability density function for the 1D Harmonic Oscillator
	at the specified positions x and parameters alpha.

	Parameters
	----------
	x : np.ndarray(N)
		Position of N particles
	alpha : np.ndarray	
		Parameter of the trial wave function

	Returns
	-------
	prob : np.ndarray(N)
		Probability density function at the specified x and alpha
	"""

	prob = np.exp(-2*alpha*x**2)

	return prob

def diff_wave_function_ratio_Harmonic_Oscillator(x, alpha):
	"""
	Returns the value Psi'(x,alpha)/Psi(x,alpha) where Psi' denotes the
	derivative of the wave function with respect to alpha for the 1D 
	Harmonic Oscillator	at the specified positions x and parameters alpha.
	The wave function Psi must be normalized.

	Parameters
	----------
	x : np.ndarray(N)
		Position of N particles
	alpha : np.ndarray	
		Parameter of the trial wave function

	Returns
	-------
	ratio : np.ndarray(N)
		Ratio Psi'(x,alpha)/Psi(x,alpha) at the specified x and alpha
	"""

	ratio = 1/(4*alpha)-x**2

	return ratio


################################
#        HYDROGEN ATOM         #
################################

# TRIAL WAVE FUNCTION = e^(-alpha*r) where r = sqrt(x^2 + y^2 + z^2)

def E_local_Hydrogen_atom(x, y, z, alpha):
	"""
	Returns the value of the local energy for the Hydrogen atom
	at the specified positions x,y,z and parameter alpha.

	Parameters
	----------
	x,y,z : np.ndarray(N)
		Position of N particles
	alpha : np.ndarray
		Parameter of the trial wave function

	Returns
	-------
	E_local : np.ndarray(N)
		Local energy at the specified x,y,z and alpha
	"""

	r = np.sqrt(x**2 + y**2 + z**2)
	E_local = -0.5*(alpha**2 - 2*alpha/r) - 1/r

	return E_local


def prob_density_Hydrogen_atom(x, y, z, alpha):
	"""
	Returns the value of the probability density function for the Hydrogen atom
	at the specified positions x,y,z and parameter alpha.

	Parameters
	----------
	x,y,z : np.ndarray(N)
		Position of N particles
	alpha : np.ndarray	
		Parameter of the trial wave function

	Returns
	-------
	prob : np.ndarray(N)
		Probability density function at the specified r and alpha
	"""

	r = np.sqrt(x**2 + y**2 + z**2)
	prob = np.exp(-2*alpha*r)

	return prob

def diff_wave_function_ratio_Hydrogen_atom(x, y, z, alpha):
	"""
	Returns the value Psi'(x,alpha)/Psi(x,alpha) where Psi' denotes the
	derivative of the wave function with respect to alpha for the Hydrogen 
	atom at the specified positions x and parameters alpha.
	The wave function Psi must be normalized.

	Parameters
	----------
	x,y,z : np.ndarray(N)
		Position of N particles
	alpha : np.ndarray	
		Parameter of the trial wave function

	Returns
	-------
	ratio : np.ndarray(N)
		Ratio Psi'(x,alpha)/Psi(x,alpha) at the specified x and alpha
	"""

	ratio = 3/(2*alpha)-np.sqrt(x**2+y**2+z**2)

	return ratio