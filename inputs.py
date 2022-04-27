import numpy as np

################################
#     EXAMPLE OF FUNCTIONS     #
################################

def prob_density(r, alpha):
	"""
	Returns the value of the probability density function at the specified positions r and parameters alpha.
	'dim' is the dimension of the integral space

	Parameters
	----------
	r : np.ndarray(N_walkers, dim)
		Positions of N_walkers walkers
	alpha : np.ndarray(N_parameters)	
		Parameters of the trial wave function

	Returns
	-------
	prob : np.ndarray(N_walkers)
		Probability density function of the different walkers at the specified r and alpha
	"""

	prob = 0

	return prob


def E_local_f(r, alpha):
	"""
	Returns the value of the local energy at the specified positions r and parameters alpha.
	'dim' is the dimension of the integral space

	Parameters
	----------
	r : np.ndarray(N_steps, N_walkers, dim)
		Positions of N_walkers walkers in different steps
	alpha : np.ndarray(N_parameters)	
		Parameters of the trial wave function

	Returns
	-------
	E_local : np.ndarray(N_steps, N_walkers)
		Local energy of the different walkers at the specified r and alpha
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
	x : np.ndarray(N_steps, N_walkers, 1) or np.ndarray(N_steps, N_walkers, 1)
		Position of N_walkers walkers
	alpha : np.ndarray(1) or float
		Parameter of the trial wave function

	Returns
	-------
	E_local : np.ndarray(N_steps, N_walkers)
		Local energy of the different walkers at the specified x and alpha
	"""

	x = x.reshape(x.shape[:-1]) # x.shape = N_steps, N_walkers
	E_local = alpha + (0.5 - 2*alpha**2)*x**2

	return E_local


def prob_density_Harmonic_Oscillator(x, alpha):
	"""
	Returns the value of the probability density function for the 1D Harmonic Oscillator
	at the specified positions x and parameters alpha.

	Parameters
	----------
	x : np.ndarray(N_walkers, 1)
		Position of N_walkers walkers
	alpha : np.ndarray(N_parameters)
		Parameter of the trial wave function

	Returns
	-------
	prob : np.ndarray(N_walkers)
		Probability density of the different walkers function at the specified x and alpha
	"""

	x = x.reshape(x.shape[:-1]) # x.shape = N_walkers
	prob = np.exp(-2*alpha*x**2)

	return prob


################################
#        HYDROGEN ATOM         #
################################

# TRIAL WAVE FUNCTION = e**(-alpha*r) where r = sqrt(x**2 + y**2 + z**2)

def E_local_Hydrogen_atom(r, alpha):
	"""
	Returns the value of the local energy for the Hydrogen atom
	at the specified positions r and parameter alpha.

	Parameters
	----------
	r : np.ndarray(Nsteps, N_walkers, 3)
		Position in 3D of N_walkers walkers (r = (x,y,z))
	alpha : np.ndarray(1) or float
		Parameter of the trial wave function

	Returns
	-------
	E_local : np.ndarray(N_walkers)
		Local energy of the different walkers at the specified r and alpha
	"""

	r = np.linalg.norm(r, axis=-1)
	E_local = -0.5*(alpha**2 - 2*alpha/r) - 1/r

	return E_local


def prob_density_Hydrogen_atom(r, alpha):
	"""
	Returns the value of the probability density function for the Hydrogen atom
	at the specified positions r and parameter alpha.

	Parameters
	----------
	r : np.ndarray(N_walkers, 3)
		Position in 3D of N_walkers walkers (r = (x,y,z))
	alpha : np.ndarray(1) or float
		Parameter of the trial wave function

	Returns
	-------
	prob : np.ndarray(N_walkers)
		Probability density function of the different walkers at the specified r and alpha
	"""

	r = np.linalg.norm(r, axis=-1)
	prob = np.exp(-2*alpha*r)

	return prob


################################
#        HELIUM ATOM         #
################################

# TRIAL WAVE FUNCTION = e**(-z*(r1+r2))e**(r12/2(1+alpha*r12)) where r1,2 = sqrt(x1,2**2 + y1,2**2 + z1,2**2) and r12 = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def WF_Helium_atom(r, alpha):
	"""
	Returns the value of the wave function for the Helium atom
	at the specified positions r and parameters alpha.

	Parameters
	----------
	r : np.ndarray(N_steps, N_walkers, 6)
		Position in 3D of N_walkers walkers (from WF: r = (x1,y1,z1,x2,y2,z2))
	alpha : np.ndarray(2)
		Parameters of the trial wave function (from WF: alpha = (z, alpha))

	Returns
	-------
	psi : np.ndarray(N_steps, N_walkers)
		Wave function at the specified position and alpha
	"""

	r1 = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
	r2 = np.sqrt(r[3]**2 + r[4]**2 + r[5]**2)
	r12 = np.sqrt((r[0]-r[3])**2 + (r[1]-r[4])**2 + (r[2]-r[5])**2)
	psi = np.exp(-alpha[0]*(r1+r2))*np.exp(r12/(2*(1+alpha[1]*r12)))

	return psi


def E_local_Helium_atom(r, alpha, h=0.001):
	"""
	Returns the value of the local energy for the Helium atom
	at the specified positions r and parameters alpha.
	It uses numerical differentiation to compute the kinetic energy.

	Parameters
	----------
	r : np.ndarray(N_steps, N_walkers, 6)
		Position in 3D of N_walkers walkers (from WF: r = (x1,y1,z1,x2,y2,z2))
	alpha : np.ndarray(2)
		Parameters of the trial wave function (from WF: alpha = (z, alpha))
	h : float
		Step for the numerical derivative

	Returns
	-------
	E_local : np.ndarray(N_steps, N_walkers)
		Local energy at the specified position and alpha
	"""

	x1, y1, z1, x2, y2, z2 = r.T # x_i = N_walkers, N_steps
	z, alpha = alpha

	# kinetic energy / WF
	d2x1 = (WF_Helium_atom(x1+h, y1, z1, x2, y2, z2, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1-h, y1, z1, x2, y2, z2, z, alpha))/h**2
	d2y1 = (WF_Helium_atom(x1, y1+h, z1, x2, y2, z2, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1, y1-h, z1, x2, y2, z2, z, alpha))/h**2
	d2z1 = (WF_Helium_atom(x1, y1, z1+h, x2, y2, z2, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1, y1, z1-h, x2, y2, z2, z, alpha))/h**2
	d2x2 = (WF_Helium_atom(x1, y1, z1, x2+h, y2, z2, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1, y1, z1, x2-h, y2, z2, z, alpha))/h**2
	d2y2 = (WF_Helium_atom(x1, y1, z1, x2, y2+h, z2, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1, y1, z1, x2, y2-h, z2, z, alpha))/h**2
	d2z2 = (WF_Helium_atom(x1, y1, z1, x2, y2, z2+h, z, alpha)-2*WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha)+WF_Helium_atom(x1, y1, z1, x2, y2, z2-h, z, alpha))/h**2
	E_kin = -0.5*(d2x1 + d2y1 + d2z1 + d2x2 + d2y2 + d2z2).T / WF_Helium_atom(x1, y1, z1, x2, y2, z2, z, alpha) # E_kin = N_steps, N_walkers

	# potential energy 
	r1 = np.sqrt(x1**2 + y1**2 + z1**2)
	r2 = np.sqrt(x2**2 + y2**2 + z2**2)
	r12 = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
	E_pot = (-2/r1 -2/r2 +1/r12).T # E_kin = N_steps, N_walkers

	return E_kin + E_pot


def prob_density_Helium_atom(r, alpha):
	"""
	Returns the value of the probability density function for the Helium atom
	at the specified positions r and parameters z and alpha.

	Parameters
	----------
	r : np.ndarray(N_walkers, 6)
		Position in 3D of N_walkers walkers (from WF: r = (x1,y1,z1,x2,y2,z2))
	alpha : np.ndarray(2)
		Parameters of the trial wave function (from WF: alpha = (z, alpha))

	Returns
	-------
	prob : np.ndarray(N_walkers)
		Probability density function at the specified position, z and alpha
	"""

	x1, y1, z1, x2, y2, z2 = r.T # x_i = N_walkers
	z, alpha = alpha

	r1 = np.sqrt(x1**2 + y1**2 + z1**2)
	r2 = np.sqrt(x2**2 + y2**2 + z2**2)
	r12 = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
	prob = np.exp(-2*z*(r1+r2))*np.exp(r12/(1+alpha*r12))

	return prob