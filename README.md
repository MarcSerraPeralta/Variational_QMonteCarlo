# Project 2: Variational Quantum Monte Carlo

Variational Quantum Monte Carlo to approximate the ground state of a quantum system given a trial wavefunction with a set of free parameters. 

The best approximation is found by minimizing the expectation value of the energy $\langle E\rangle (\vec{\alpha})$ with respect to $\vec{\alpha}$, where $\vec{\alpha}$ is the set of free parameters. The expectation value of the energy is computed using Monte Carlo integration withm importance sampling method, where the sampling from the probability density function is done using the Metropolis algorithm. To find the minimum of the expectation value can be done performing a scan of a region of parameter space or using the steepest descent method.

The necessary problem-related inputs are the probability density function $|\Psi_\text{T}(\vec{r},\vec{\alpha})|^2$ and the local energy function $E_\text{loc}(\vec{r}, \vec{\alpha}) = \frac{H\Psi_\text{T}(\vec{r},\vec{\alpha})}{\Psi_\text{T}(\vec{r},\vec{\alpha})}$. The included examples of problem systems are:

- 1D Harmonic Oscillator
- Hydrogen atom
- Helium atom (ground state, first excited state, second excited state, third excited state, fourth excited state)

The calculation for the expectation value of the energy $\langle E\rangle (\vec{\alpha})$ is parallelized using the _multiprocessing_ Python library. 

## Setup

Clone this repo and run `pip install -r requirements.txt` to install its dependencies.

**Important note**: For Windows users only, `if __name__ == '__main__':` has to be added to `main.py` right after the import of libraries. Therefore, the code below this statement has to be indented one tab to the right. This solves some issues coming from the _multiprocessing_ library.

## Usage

Open `main.py` and specify the input parameters:

- Hamiltonian and trial wave function
    - `E_local_f`: defined as $E_{local} = (\hat{H}\Psi)/\Psi$, specific of each problem and trial wavefunction $\Psi$. The functions for the examples are stored in `inputs.py`. 
    - `prob_density`: probability density function of the trial wavefunction $p_{density} = |\Psi|^2$. The functions for the examples are stored in `inputs.py`. 
    - `dim`: number of degrees of freedom in the system
    - `opt_args`: defines the optimization algorithm to be used for finding the minimum expectation value of the energy and the parameters associated with it (stored in a dictionary object). The available algorithms and its parameters are defined below.
- Monte Carlo integration parameters
    - `N_steps`: number of points to extract from the probability density function for each walker.
    - `N_walkers`: number of walkers. Each walker has a randomly generated initial point according to a normal distribution of standard deviation `L_start`.
    - `N_skip`: number of initial steps to remove from each walker to avoid bias of the initial point.
    - `L_start`: standard deviation of the normal distribution with which each walker is initialized.
    - `N_cores`: number of cores used for the computation of the expectation value of the energy.
- Results saving parameters
    - `file_name`: name of the CSV file in which to store the computed expectation values of the energy.

## Algorithms for finding the minimum of energy

The following algorithms are included to find the minimum of energy:

- `scan1D`: scans the parameter space (1 free parameter) from an initial alpha to a final alpha with a given step. The parameters to specify are: initial alpha (`init_alpha`), step (`step`), final alpha (`final_alpha`).
- `steepest_descent1D`: searches for the minimum expectation value of the energy using gradient descent (1 free parameter). The parameters to specify are: initial alpha (`init_alpha`), step to be used for the first computation of the derivative (`init_step`), damping parameter parameter for $\alpha$, i.e. $\alpha_{\text{new}}=\alpha_{\text{old}}-\gamma \frac{dE}{d\alpha}$, (`gamma`), desired precision with which to find the local minimum, i.e. $|(\alpha_{\text{new}} - \alpha_{\text{old}}) / \alpha_{\text{old}}| < \epsilon$ (`precision`). 
- `steepest_descent_ANY_D`: searches for the minimum expectation value of the energy using gradient descent (arbitrary number of free parameters). Same parameters to specify as `steepest_descent1D`. 

## Examples 

See beginning of `inputs.py` for the full description of the examples. To implement them in `main.py`, copy the following blocks of code

```
# HARMONIC OSCILLATOR (SCAN)
E_local_f = E_local_Harmonic_Oscillator
prob_density = prob_density_Harmonic_Oscillator
dim = 1 # dimension of configuration space
opt_args = {"method":"scan1D", "init_alpha":np.array([0.25]), "step":0.05, "final":np.array([0.75])}

# HYDROGEN ATOM (SCAN, STEEPEST DESCENT IN 1D AND ANY D)
E_local_f = E_local_Hydrogen_atom
prob_density = prob_density_Hydrogen_atom
dim = 3 # dimension of configuration space
opt_args = {"method":"scan1D", "init_alpha":np.array([0.75]), "step":0.05, "final":np.array([1.25])}
opt_args = {"method":"steepest_descent1D", "init_alpha":np.array([0.75]), "init_step":0.05, "gamma":0.25, "precision":1E-3}
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([0.75]), "init_step":np.array([0.05]), "gamma":0.25, "precision":1E-3}

# HELIUM ATOM: GROUND STATE 1 PARAMETER (SCAN, STEEPEST DESCENT IN 1D AND ANY D)
E_local_f = E_local_Helium_atom_GS_1param_numeric
E_local_f = E_local_Helium_atom_GS_1param_analytic
prob_density = prob_density_Helium_atom_GS_1param
dim = 6
opt_args = {"method":"scan1D", "init_alpha":np.array([0.15]), "step":0.005, "final":np.array([0.17])}
opt_args = {"method":"steepest_descent1D", "init_alpha":np.array([0.15]), "init_step":0.01, "gamma":0.03, "precision":1E-5}
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([0.15]), "init_step":np.array([0.01]), "gamma":0.03, "precision":1E-5}

# HELIUM ATOM: GROUND STATE 2 PARAMETERS (STEEPEST DESCENT)
E_local_f = E_local_Helium_atom_GS
prob_density = prob_density_Helium_atom_GS
dim = 6
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([2,0.5]), "init_step":np.array([0.1,-0.1]), "gamma":0.5, "precision":1E-5}

# HELIUM ATOM: FIRST EXCITED 3 PARAMETERS (STEEPEST DESCENT)
E_local_f = E_local_Helium_atom_1E
prob_density = prob_density_Helium_atom_1E
dim = 6
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([2,2,0.5]), "init_step":np.array([-0.1,-0.1,-0.1]), "gamma":0.5, "precision":1E-5}

# HELIUM ATOM: SECOND EXCITED 3 PARAMETERS (STEEPEST DESCENT)
E_local_f = E_local_Helium_atom_2E
prob_density = prob_density_Helium_atom_2E
dim = 6
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([2,2,0.5]), "init_step":np.array([-0.1,-0.1,-0.1]), "gamma":0.5, "precision":1E-5}

# HELIUM ATOM: THIRD EXCITED 3 PARAMETERS (STEEPEST DESCENT)
E_local_f = E_local_Helium_atom_3E
prob_density = prob_density_Helium_atom_3E
dim = 6
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([2,2,0.5]), "init_step":np.array([-0.1,-0.1,-0.1]), "gamma":0.5, "precision":1E-5}

# HELIUM ATOM: FOURTH EXCITED 3 PARAMETERS (STEEPEST DESCENT)
E_local_f = E_local_Helium_atom_4E
prob_density = prob_density_Helium_atom_4E
dim = 6
opt_args = {"method":"steepest_descent_ANY_D", "init_alpha":np.array([2,2,0.5]), "init_step":np.array([-0.1,-0.1,-0.1]), "gamma":0.5, "precision":1E-5}
```

## Results and discusion

The whole project is described in detail in [report.pdf](https://github.com/abermejillo/variational_QMonteCarlo/blob/40cf9892ab1150e591f16fd2642e5f0dcf75af05/report.pdf). The results of the spectra of Helium are shown and discussed.

## Authors 
- Álvaro Bermejillo Seco
- Daniel Bedialauneta Rodríguez
- Marc Serra Peralta
