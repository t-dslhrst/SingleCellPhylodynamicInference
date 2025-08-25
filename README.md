# Single Cell Phylodynamic Inference
 
Implementation of the algorithm published in [[Dieselhorst and Berg 2024](https://doi.org/10.1101/2024.12.05.627005 )] for phylodynamic inference from mutations accumulated under cellular reproduction. All code is written in the [Julia programming language](https://julialang.org/) and tested on Julia version 1.11.3 (but can also be executed with earlier Julia versions).

- `project.toml`and `manifest.toml` list the dependencies and package versions on which the code was developed.
- `parameter_inference.jl` contains the algorithm.
- `simulate_population.jl` is a simple implementation to simulate a population and reconstruct the phylogenetic tree.
- `simulate_and_infer.jl` gives an example for simulation based assessment of the inference scheme. A population is simulated, the phylogenetic tree is reconstructed, the parameters are inferred and results are saved in a CSV file. The underlying parameters and discretizations can be set in the first few lines of the script. Running this code usually takes around three minutes (4.4GHz clock speed).


In order to use this code for inference from empirical data, the phylogenetic tree must be defined as a [SimpleWeightedGraph](https://juliagraphs.org/SimpleWeightedGraphs.jl/dev/) with integer weights corresponding to the number of mutations. Note that branches of length zero cannot be added after initializing the graph. We give an example workflow with the following components

- `read_data` contains example data (`example_data.tree` from simulations) and a Python script (`newick_to_csv.py`) to create the corresponding adjacency matrix (`example_data_adj_matrix_p1.csv`) from a Newick tree. In order to account for branches with length zero, the adjacency matrix is saved with one additional mutation per branch, such that the minimum branch length is one.
- `inference_from_adj_matrix.jl` shows how to perform inference based on a tree saved as as an adjacency matrix in CSV format (with branch lengths increased by one mutation). Inference results are saved in a CSV file. Running this code usually takes around three minutes (4.4GHz clock speed).
