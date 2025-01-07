# Single Cell Phylodynamic Inference
 
Implementation of the algorithm published in [[Dieselhorst and Berg 2024](https://doi.org/10.1101/2024.12.05.627005 )] for phylodynamic inference from mutations accumulated under cellular reproduction. All code is written in the [Julia programming language](https://julialang.org/) and was tested on Julia version 1.11.1 (but can also be executed with earlier Julia versions).

- `project.toml`and `manifest.toml` list the dependencies and package versions on which the code was developed.
- `parameter_inference.jl` contains the algorithm.
- `simulate_population.jl` is a simple implementation to simulate a population and reconstruct the phylogenetic tree.
- `simulate_and_infer.jl` gives an example on how to infer the parameters of a simulated population. The inference results are saved in a CSV file. Running this code usually takes around ten minutes (4GHz clock speed).

In order to use this code for inference from empirical data, the phylogenetic tree must be defined as a [SimpleWeightedGraph](https://juliagraphs.org/SimpleWeightedGraphs.jl/dev/) with integer weights corresponding to the number of mutations. Note that branches of length zero cannot be added after initializing the graph.
