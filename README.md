# Project 2: Variational Quantum Monte Carlo



## Setup

Clone this repo and run `pip install -r requirements.txt` to install its dependencies.

**Important note**: For Windows users only, `if __name__ == '__main__':` has to be added to `main.py` right after the import of libraries. Therefore, the code below has to be indented one tab to the right. This is so that the multiprocessing of different walkers works correctly.

## Usage

Open `main.py` and specify the input parameters:

- Hamiltonian and trial wave function
    - `E_local_f` and `prob_density`: $`E_{local} = \Psi^{-1} \cdot \hat{H}\Psi`$ and $`p_{density} = |\Psi|^2`$, specific of each problem and trial wavefunction. The functions are stored in `inputs.py`. 
    - `dim`: number of degrees of freedom in the system
    - `opt_args`: defines the optimization algorithm to be used for finding the minimum expectation value of the energy and the parameters associated in a dictionary object. The available algorithms and its parameters are defined below.
- Monte Carlo integration parameters
    - `N_steps`: number of points to extract from the probability density function for each walker.
    - `N_walkers`: number of walkers. Each walker has a randomly generated initial point according to a normal distribution of standard deviation `L_start`.
    - `N_skip`: number of initial steps to remove from each walker to avoid bias of the initial point.
    - `L_start`: standard deviation of the normal distribution with which each walker is initialized.
- Results saving parameters
    - `file_name`: name of the file in which to store the computed expectation values of the energy.

## Optimization algorithms
- scan1D: Scans the parameter space from `init_alpha` to `final_alpha` with a given step. The parameters to specify are:
    - Initial alpha: this will be given by `init_alpha`.
    - Step: to be specified by user.
    - Final alpha: this will be given by `final_alpha`

## Authors 
- Álvaro Bermejillo Seco
- Daniel Bedialauneta Rodríguez
- Marc Serra Peralta
