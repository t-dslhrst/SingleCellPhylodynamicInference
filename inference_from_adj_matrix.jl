include("parameter_inference.jl")
using Tables, CSV, StatsBase, DataFrames

N_sample = 300
ρ = 1.0
N_τ = 1000
i_max = 50
start_norm = 20.0

function run()
    # Read the adjacency matrix where all branch lengths are increased by one do avoid branch lengths of zero.
    weights_p1 = Matrix(CSV.read("read_data/example_data_adj_matrix_p1.csv", DataFrame, header=false))
    # Initialize a sparse matrix as adjacency matrix of the tree with Float64 weights
    # Subtract a number just below one from all branch lengths from the csv file and add them to the sparse matrix
    weights = convert(SparseMatrixCSC, weights_p1 .- (weights_p1 .> 0) .* 0.9999999999)
    weights_p1 = nothing
    # Change all elements to integer numbers (still in Float64 format)
    for (i, weight) in enumerate(weights)
        if weight < 1
            if weight > 0
                weights[i] = 0.0
            end
        end
    end
    weights[weights .>= 1] .= floor.(weights[weights .>= 1])
    # create a tree with Int32 weights
    tree = SimpleWeightedGraph(convert(SparseMatrixCSC{Int32, Int32}, weights))

    # infer parameters
    params_inf = infer_parameters(tree, ρ, N_τ=N_τ, norm=start_norm, i_max=i_max)

    filename = "example_data_inference_results.csv"
    if !isfile(filename)
        CSV.write(filename, Tables.table(params_inf'), header=["q", "mu"])
    else
        CSV.write(filename, Tables.table(params_inf'), append=true)
    end
end

run()