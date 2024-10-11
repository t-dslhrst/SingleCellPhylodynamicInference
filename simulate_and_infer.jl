include("parameter_inference.jl")
using Tables, CSV, StatsBase

N_sample = 300
ρ = 1.0
q = 0.75
μ = 1.0
N_τ = 1000
i_max = 50
start_norm = 20.0

function run()
    living, dead, _, _ = grow_population(round(Int64, N_sample/ρ), 1.0, q, μ, 0.001)
    tree = build_tree(vcat(living, dead))
    pop_sample = sample(living[:,1], round(Int64, N_sample), replace=false)
    living, dead = nothing, nothing

    rtree, _ = reduce_tree(tree, pop_sample)
    tree, pop_sample = nothing, nothing


    params_inf = infer_parameters(rtree, ρ, N_τ=N_τ, norm=start_norm, i_max=i_max)
    
    filename = "rho$(ρ)_N$(N_sample)_q$(q)_mu$(μ).csv"
    if !isfile(filename)
        CSV.write(filename, Tables.table(params_inf'), header=["q", "mu"])
    else
        CSV.write(filename, Tables.table(params_inf'), append=true)
    end
end

run()