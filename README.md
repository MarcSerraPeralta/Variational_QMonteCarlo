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
- scan1D: Scans the parameter space from an initial alpha to a final alpha with a given step. The parameters to specify are:
    - Initial alpha: this will be given in `"init_alpha"`.
    - Step: this will be given in `"step"`
    - Final alpha: this will be given in `"final_alpha"`.
- steepest_descent1D: Searches for the minimum expectation value of the energy using gradient descent. The parameters to specify are:
    - Initial alpha: this will be given in `"init_alpha"`.
    - Initial step: The derivative of the energy with respect to alpha is computed numerically. The initial step indicates the displacement in alpha to be used for the first computation of the derivative. This will be given in `"init_step"`.
    - Gamma: The parameter $`\alpha`$ is updated in the following way: $`\alpha_{\text{new}}=\alpha_{\text{old}}-\gamma \frac{dE}{d\alpha}`$. The $`\gamma`$ factor serves the purpose of damping the gradient descent. This will be given in `"gamma"`.
    - Precision: Desired precision with which to find the local minimum. Let $`\epsilon`$ be this precision, then the minimization algorithm stops when $`|(\alpha_{\text{new}}-\alpha_{\text{old}})/\alpha_{\text{old}}|<\epsilon`$. This will be given in `"precision"`.

## Authors 
- Álvaro Bermejillo Seco
- Daniel Bedialauneta Rodríguez
- Marc Serra Peralta
