# Single Cell Phylodynamic Inference
 
Implementation of the algorithm published in [[Dieselhorst and Berg 2024](https://doi.org/10.1101/2024.12.05.627005 )] for phylodynamic inference from mutations accumulated under cellular reproduction. All code is written in the [Julia programming language](https://julialang.org/) and tested on Julia version 1.11.3 (but can also be executed with earlier Julia versions).

- `project.toml`and `manifest.toml` list the dependencies and package versions on which the code was developed.
- `parameter_inference.jl` contains the algorithm. See below for a detailed explanation of the important functions.
- `simulate_population.jl` is a simple implementation to simulate a population and reconstruct the phylogenetic tree.
- `simulate_and_infer.jl` gives an example for simulation based assessment of the inference scheme. A population is simulated, the phylogenetic tree is reconstructed, the parameters are inferred and results are saved in a CSV file. The underlying parameters and discretizations can be set in the first few lines of the script. Running this code usually takes around three minutes (4.4GHz clock speed).

## Data formatting
In order to use this code for inference from empirical data, the phylogenetic tree must be defined as a [SimpleWeightedGraph](https://juliagraphs.org/SimpleWeightedGraphs.jl/dev/) with integer weights corresponding to the number of mutations. Note that branches of length zero cannot be added after initializing the graph. We give an example workflow with the following components

- `read_data` contains example data (`example_data.tree` from simulations) and a Python script (`newick_to_csv.py`) to create the corresponding adjacency matrix (`example_data_adj_matrix_p1.csv`) from a Newick tree. In order to account for branches with length zero, the adjacency matrix is saved with one additional mutation per branch, such that the minimum branch length is one.
- `inference_from_adj_matrix.jl` shows how to perform inference based on a tree saved as as an adjacency matrix in CSV format (with branch lengths increased by one mutation). Inference results are saved in a CSV file. Running this code usually takes around three minutes (4.4GHz clock speed).

## Main functions
Our implementations provide docstrings for all functions. Here, we introduce the two main functions, i.e. for the likelihood computation and the parameter inference, which can be found in `parameter_inference.jl`.
### Likelihood computation
```julia
function compute_param_likelihood(tree::SimpleWeightedGraph{T, T}, δ::Real, μ::Real;
   β::Real=1.0, ρ::Real=1.0, μ_e::Real=0, N_τ::Integer=1000, τ_max::Real=0.0, τ_max_fact::Real=2.0, i_max::Integer=100,
   norm::Real=20.0,integral::Function=trapezoidal_integral, return_without_integrating::Bool=false) where T<:Integer
```
Computes the parameter log-likelihood for a given tree.

#### Arguments:
- `tree`: Phylogenetic tree (as a `SimpleWeightedGraph`).
- `δ`: Death rate.
- `μ`: Mean number of mutations per birth event.
- `β`: Birth rate (default is 1.0).
- `ρ`: Sampling probability (default is 1.0).
- `μ_e`: Mean amplification error (default is 0)
- `N_τ`: Number of time intervals for integration (default is 1000).
- `τ_max`: Maximum time for integration (computed as τ_max_fact*(log(length(leaves)/ρ)) / (β-δ) if not provided).
- `τ_max_fact`: Factor to scale maximum time (default is 2.0).
- `i_max`: Maximum number of generations for calculation (default is 100).
- `norm`: Normalization factor for p_1(i|τ_s,τ_e) to avoid vanishing probabilities (default is 20.0).
- `integral`: Integration function (default is `trapezoidal_integral`).
- `return_without_integrating`: do not integrate over the time of the MRCA (needed if the tree is a subtree of a heterogeneous tree)

#### Returns:
- The log-likelihood value of the parameters for the given tree.

### Parameter inference
```julia
function infer_parameters(tree::SimpleWeightedGraph{T, T}, ρ::Real=1.0, μ_e::Real=0;
   q_0::Real=rand(), μ_0::Real=2*rand(), ρ_0::Real=rand(), μ_e_0::Real=2*rand(),
   infer_ρ::Bool=false, infer_μ_e::Bool=false, return_likelihood::Bool=false, show_trace::Bool=false, norm::Real=20.0, L_tol::Real=1e-8, kwargs...) where T<:Integer
```
Infers most likely parameters of a given tree using Nelder-Mead optimization.

#### Arguments:
- `tree`: Phylogenetic tree.
- `ρ`: Sampling probability (default is 1.0).
- `μ_e`: Mean of the Poisson distributed amplification error (default is 0)
- `q_0`: Initial relative death rate (default is random).
- `μ_0`: Initial mean number of mutations per birth event. (default is random).
- `ρ_0`: Initial sampling probability, if inferred (default is random).
- `μ_e_0`: Initial amplification error mean (default is random).
- `infer_ρ`: Whether to infer sampling probability `ρ` (default is false).
- `infer_μ_e`: Whether to infer the mean amplification error `ρ` (default is false).
- `return_likelihood`: Return likelihood alongside parameters (default is false).
- `show_trace`: Show optimization trace (default is false).
- `norm`: Initial normalization factor (default is 20.0).
- `L_tol`: Tolerance for likelihood convergence (default is 1e-8).
- `kwargs`: Additional arguments passed to likelihood computation.

#### Returns:
- The inferred parameters (and optionally their likelihood).
