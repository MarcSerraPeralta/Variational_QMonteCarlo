# Project 2: Variational Quantum Monte Carlo

Variational Quantum Monte Carlo to approximate the ground state of a quantum system given a trial wavefunction with a set of free parameters. 

The best approximation is found by minimizing the expectation value of the energy $`\langle E\rangle (\bm{\alpha})`$ with respect to $`\bm{\alpha}`$, where $`\bm{\alpha}`$ is the set of free parameters. The expectation value of the energy is computed using Monte Carlo integration withm importance sampling method, where the sampling from the probability density function is done using the Metropolis algorithm. To find the minimum of the expectation value can be done performing a scan of a region of parameter space or using the steepest descent method.

The necessary problem-related inputs are the probability density function $`|\Psi_\text{T}(\bm{r},\bm{\alpha})|^2`$ and the local energy function $`E_\text{loc}(\bm{r}, \bm{\alpha}) = \frac{H\Psi_\text{T}(\bm{r},\bm{\alpha})}{\Psi_\text{T}(\bm{r},\bm{\alpha})}`$. The included examples of problem systems are:

- 1D Harmonic Oscillator
- Hydrogen atom
- Helium atom (ground state, first excited state, second excited state, third excited state, fourth excited state)

The calculation for the expectation value of the energy $`\langle E\rangle (\bm{\alpha})`$ is parallelized using the _multiprocessing_ Python library. 

## Setup

Clone this repo and run `pip install -r requirements.txt` to install its dependencies.

**Important note**: For Windows users only, `if __name__ == '__main__':` has to be added to `main.py` right after the import of libraries. Therefore, the code below this statement has to be indented one tab to the right. This solves some issues coming from the _multiprocessing_ library.

## Usage

Open `main.py` and specify the input parameters:

- Hamiltonian and trial wave function
    - `E_local_f`: defined as $`E_{local} = (\hat{H}\Psi)/\Psi`$, specific of each problem and trial wavefunction $`\Psi`$. The functions for the examples are stored in `inputs.py`. 
    - `prob_density`: probability density function of the trial wavefunction $`p_{density} = |\Psi|^2`$. The functions for the examples are stored in `inputs.py`. 
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
- `steepest_descent1D`: searches for the minimum expectation value of the energy using gradient descent (1 free parameter). The parameters to specify are: initial alpha (`init_alpha`), step to be used for the first computation of the derivative (`init_step`), damping parameter parameter for $`\alpha`$, i.e. $`\alpha_{\text{new}}=\alpha_{\text{old}}-\gamma \frac{dE}{d\alpha}`$, (`gamma`), desired precision with which to find the local minimum, i.e. $`|(\alpha_{\text{new}} - \alpha_{\text{old}}) / \alpha_{\text{old}}| < \epsilon`$ (`precision`). 
- `steepest_descent_ANY_D`: searches for the minimum expectation value of the energy using gradient descent (arbitrary number of free parameters). Same parameters to specify as `steepest_descent1D`. 

## Examples 

See beginning of `inputs.py` for the full description of the examples. 

## Authors 
- Álvaro Bermejillo Seco
- Daniel Bedialauneta Rodríguez
- Marc Serra Peralta
