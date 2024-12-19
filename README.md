# Single Cell Phylodynamic Inference
 
Implementation of the algorithm published in [[Dieselhorst and Berg 2024](https://doi.org/10.1101/2024.12.05.627005 )] for phylodynamic inference from mutations accumulated under cellular reproduction.

- `parameter_inference.jl` contains the algorithm.
- `simulate_population.jl` is a simple implementation to simulate a population and reconstruct the phylogenetic tree.
- `simulate_and_infer.jl` gives an example on how to infer the parameters of a simulated population.
