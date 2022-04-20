# Project 2: Variational Quantum Monte Carlo

The README.md file serves as a reference for other users visiting your repository.
It should contain a brief description of your project, and document the steps others need to take to get your application up and running.
In addition, it should list the authors of the project.

## Setup

Clone this repo and run `pip install -r requirements.txt` to install its dependencies.

## Usage

Open `main.py` and specify the input parameters:

- Hamiltonian and trial wave function (with parameter $`\alpha`$ to be optimized) defined with symbols from simpy
    - `x` and `alpha`: define the degrees of freedom in the system and the parameter space with sympy symbols.
    - `psi_t`: trial wave function in terms of `x` and `alpha`. Does not need to be normalized.
    - `H_psi_t`: the hamiltonian acting on `psi_t`. It can contain the `diff` operator.
    - `var`: list with all simpy defined variables (`x` and `alpha`)
    - `dim`: number of degrees of freedom in the system
    - `init_alpha` and `final_alpha`: range of parameter space where to look for the minimum expected energy.
- Monte Carlo integration parameters
    - `N_steps`: number of points to extract from the probability density function for each walker.
    - `N_walkers`: number of walkers. Each walker has a randomly generated initial point according to a normal distribution of standard deviation `L_start`.
    - `N_skip`: number of initial steps to remove from each walker to avoid bias of the initial point.
    - `L_start`: standard deviation of the normal distribution with which each walker is initialized.
- Alpha optimization parameters
    - `opt_args`: defines the optimization algorithm to be used for finding the minimum expectation value of the energy and the parameters associated in a dictionary object. The available algorithms and its parameters are defined below.
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
