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

# TRIAL WAVE FUNCTION = e**(-alpha*x**2)

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

# TRIAL WAVE FUNCTION = e**(-alpha*r) where r = sqrt(x**2 + y**2 + z**2)

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

################################
#        HELIUM ATOM         #
################################

# TRIAL WAVE FUNCTION = e**(-z*(r1+r2))e**(r12/2(1+alpha*r12)) where r1,2 = sqrt(x1,2**2 + y1,2**2 + z1,2**2) and r12 = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha):
	"""
	Returns the value of the wave function for the Helium atom
	at the specified positions x1,y1,z1,x2,y2,z2 and parameters z and alpha.

	Parameters
	----------
	x1,y1,z1,x2,y2,z2 : np.ndarray(N)
		Position of N particles
	z : np.darray
		Parameter of the trial wave function
	alpha : np.ndarray
		Parameter of the trial wave function
	Returns
	-------
	psi : np.ndarray(N)
		Wave function at the specified position, z and alpha
	"""
	r1 = np.sqrt(x1**2 + y1**2 + z1**2)
	r2 = np.sqrt(x2**2 + y2**2 + z2**2)
	r12 = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
	psi = np.exp(-z*(r1+r2))*np.exp(r12/(2*(1+alpha*r12)))
	return psi

def E_local_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha, h=0.001):
	"""
	Returns the value of the local energy for the Helium atom
	at the specified positions x1,y1,z1,x2,y2,z2 and parameters z and alpha.
	It uses numerical differentiation to compute the kinetic energy.

	Parameters
	----------
	x1,y1,z1,x2,y2,z2 : np.ndarray(N)
		Position of N particles
	z : np.darray
		Parameter of the trial wave function
	alpha : np.ndarray
		Parameter of the trial wave function
	h : float
		Step of the numerical derivative
	Returns
	-------
	E_local : np.ndarray(N)
		Local energy at the specified position, z and alpha
	"""
	d2x1 = (WF_Helium_atom(x1+h, y1, z1, x2, y2, z2, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1-h, y1, z1, x2, y2, z2, z, alpha))/h**2
	d2y1 = (WF_Helium_atom(x1, y1+h, z1, x2, y2, z2, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1, y1-h, z1, x2, y2, z2, z, alpha))/h**2
	d2z1 = (WF_Helium_atom(x1, y1, z1+h, x2, y2, z2, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1, y1, z1-h, x2, y2, z2, z, alpha))/h**2
	d2x2 = (WF_Helium_atom(x1, y1, z1, x2+h, y2, z2, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1, y1, z1, x2-h, y2, z2, z, alpha))/h**2
	d2y2 = (WF_Helium_atom(x1, y1, z1, x2, y2+h, z2, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1, y1, z1, x2, y2-h, z2, z, alpha))/h**2
	d2z2 = (WF_Helium_atom(x1, y1, z1, x2, y2, z2+h, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1, y1, z1, x2, y2, z2-h, z, alpha))/h**2
	E_kin = -0.5*(d2x1 + d2y1 + d2z1 + d2x2 + d2y2 + d2z2)

	r1 = np.sqrt(x1**2 + y1**2 + z1**2)
	r2 = np.sqrt(x2**2 + y2**2 + z2**2)
	r12 = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
	E_pot = (-2/r1 -2/r2 +1/r12)*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)

	return E_kin + E_pot

def prob_density_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha):
	"""
	Returns the value of the probability density function for the Helium atom
	at the specified positions x1,y1,z1,x2,y2,z2 and parameters z and alpha.

	Parameters
	----------
	x1,y1,z1,x2,y2,z2 : np.ndarray(N)
		Position of N particles
		z : np.darray
		Parameter of the trial wave function
	alpha : np.ndarray	
		Parameter of the trial wave function

	Returns
	-------
	prob : np.ndarray(N)
		Probability density function at the specified position, z and alpha
	"""
	r1 = np.sqrt(x1**2 + y1**2 + z1**2)
	r2 = np.sqrt(x2**2 + y2**2 + z2**2)
	r12 = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
	prob = np.exp(-2*z*(r1+r2))*np.exp(r12/(1+alpha*r12))

	return prob