include("simulate_population.jl")
using LinearAlgebra, Optim

# ChatGPT was used to improve docstrings

"""
Compute integral using the trapezoidal rule.

# Arguments:
- `f_x`: Vector of function values to integrate.
- `dx`: The spacing between points in `f_x` (assumed to be uniform).

# Returns:
- The approximate integral.
"""
function trapezoidal_integral(f_x::Vector{T}, dx::T) where T<:Real
    if length(f_x) == 1
        return 0.0
    else
        temp = f_x[1]
        temp += 2*sum(f_x[2:length(f_x)-1])
        temp += f_x[end]
        return 0.5 * temp * dx
    end
end


"""
Calculate the mean number of generations on a branch going from time `τ_s` to `τ_e`.

# Arguments:
- `τ_s`: Start time.
- `τ_e`: End time.
- `δ`: Death rate.
- `ρ`: Sampling probability.
- `β`: Birth rate (default is 1.0).

# Returns:
- Mean number of generations.
"""
function mean_generations(τ_s::Real, τ_e::Real, δ::Real, ρ::Real, β::Real=1.0)
    int_s = 2*δ*τ_s - 2*log( (β-δ-β*ρ)*exp(-(β-δ)*τ_s) + β*ρ )
    int_e = 2*δ*τ_e - 2*log( (β-δ-β*ρ)*exp(-(β-δ)*τ_e) + β*ρ )
    return int_s - int_e
end


# precompute the logarithm of factorials
const log_facts = cumsum(log.(1:1000000))


"""
Computed the joint probability p_1(i|τ_s, τ_e) that a lineage starting at time `τ_s` and runs until time `τ_e` with `i` generations in between (eq. (1) & (7)).

# Arguments:
- `i`: Number of generations.
- `τ_s`: Start time.
- `τ_e`: End time.
- `δ`: Death rate.
- `ρ`: Sampling probability.
- `β`: Birth rate (default is 1.0).

# Returns:
- The probability of having `i` generations.
"""
function p_lineage_and_gens(i::Integer, τ_s::Real, τ_e::Real, δ::Real, ρ::Real, β::Real=1.0)
    if i==0
        return exp(-(β+δ)*(τ_s-τ_e))
    elseif i <= 20
        expo = -(β+δ)*(τ_s-τ_e)
        # if the exponent is too large, compute the logarithm of the terms
        if expo <= - 600.0
            log_result = expo + i*log(mean_generations(τ_s, τ_e, δ, ρ, β)) - log_facts[i]
            return exp(log_result)
        else
            return exp(expo) * mean_generations(τ_s, τ_e, δ, ρ, β)^i / factorial(i)
        end
    else
        # if i is too large, compute the logarithm of the terms
        log_result = -(β+δ)*(τ_s-τ_e) + i*log(mean_generations(τ_s, τ_e, δ, ρ, β)) - log_facts[i]
        return exp(log_result)
    end
end


"""
Poisson distribution for the number of mutations `m` given the number of generations `i` (eq. (2)).

# Arguments:
- `m`: Number of mutations.
- `i`: Number of generations.
- `μ`: Mean mutations per birth event.

# Returns:
- The probability of observing `m` mutations.
"""
function p_m_given_i(m::T, i::Integer, μ::Real) where T<:Integer
    if m <= 20
        return exp(-(i+1)*μ) * ((i+1)*μ)^m / factorial(m)
    else
        log_result = -(i+1)*μ + m*log((i+1)*μ) - log_facts[m]
        return exp(log_result)
    end
end



"""
Compute the parameter log-likelihood for a given tree.

# Arguments:
- `tree`: Phylogenetic tree (as a `SimpleWeightedGraph`).
- `δ`: Death rate.
- `μ`: Mutation rate.
- `β`: Birth rate (default is 1.0).
- `ρ`: Sampling probability (default is 1.0).
- `N_τ`: Number of time intervals for integration (default is 1000).
- `τ_max`: Maximum time for integration (computed as τ_max_fact*(log(length(leaves)/ρ)) / (β-δ) if not provided).
- `τ_max_fact`: Factor to scale maximum time (default is 2.0).
- `norm`: Normalization factor for p_1(i|τ_s,τ_e) to avoid vanishing probabilities (default is 20.0).
- `i_max`: Maximum number of generations for calculation (default is 100).
- `integral`: Integration function (default is `trapezoidal_integral`).

# Returns:
- The log-likelihood value of the parameters for the given tree.
"""
function compute_param_likelihood(tree::SimpleWeightedGraph{T, T}, δ::Real, μ::Real; β::Real=1.0, ρ::Real=1.0, N_τ::Integer=1000, τ_max::Real=0.0, τ_max_fact::Real=2.0, norm::Real=20.0, i_max::Integer=100, integral::Function=trapezoidal_integral) where T<:Integer
    # Find the leaves
    leaves = findall(x -> x==1, degree(tree))
    # =Find the root and its direct descendants
    root = findfirst(x -> x==2, degree(tree))
    root_desc = neighbors(tree, root)

    # Set the maximum time for integration if not provided
    if τ_max==0.0
        τ_max = τ_max_fact * (log(length(leaves)/ρ)) / (β-δ)
    end
    # Create an equidistant vector for the time
    dτ = τ_max/(N_τ-1)
    τ_vec = collect(0.0:dτ:τ_max+dτ/2)
    # Create a vector for the number of generations
    i_vec = collect(0:i_max)

    # Matrix for p_1(i|τ_s, τ_e), i.e. the second term in P(ν|τ_s) (eqs. (3) and (4))
    P_lineage_and_gens = zeros(Float64, length(i_vec), N_τ, N_τ)
    for (n_i, i) in enumerate(i_vec) # Iterate over all i
        for (n_e, τ_e) in enumerate(τ_vec) # Iterate over τ_e
            @inbounds P_lineage_and_gens[n_i,n_e:end,n_e] .= p_lineage_and_gens.(i, τ_vec[n_e:end], τ_e, δ, ρ, β) * norm
        end
    end
    
    # Identify all numbers of mutations along branches in the tree
    ms = sort(unique(tree.weights))
    mp1_to_ms = zeros(Int64, ms[end]+1) # Vector to connect the number of mutations m to the corresponding element in vector ms. Element mp1_to_ms[m+1] will correspond to the index for m mutations
    for (j, m) in enumerate(ms)
        mp1_to_ms[m+1] = j
    end    

    # Matrix for the sum over P(m|i)*p_1(i|τ_s, τ_e), i.e. the first two terms in P(ν|τ_s) (eqs. (3) and (4))
    P_branch_only =  zeros(Float64, length(ms), N_τ, N_τ)
    for (j, m) in enumerate(ms) # Iterate over all relevant values of m
        p_m_given_i_vec = p_m_given_i.(m, i_vec, μ) # Store all values of P(m|i) for the given m
        for n_τ in 1:N_τ # Iterate over all discrete τ_e
            @inbounds P_branch_only[j,n_τ:end,n_τ] .= sum(view(P_lineage_and_gens, :, n_τ:N_τ, n_τ) .* p_m_given_i_vec, dims=1)[1,:]
        end
    end

    P_lineage_and_gens = nothing # Free memory
    
    # Function to compute P(ν|τ_s) for clades starting with an interior branch (eq. (4))
    function p_interior_branch!(P_nu, m::T, p_below_1::Vector{Float64}, p_below_2::Vector{Float64})
        for n_τ in 1:N_τ
            @inbounds P_nu[n_τ] = 2*β * integral(P_branch_only[mp1_to_ms[m+1], n_τ, 1:n_τ] .* p_below_1[1:n_τ] .* p_below_2[1:n_τ], dτ)
        end
    end

    # Matrix to store P(ν|τ_s) for all clades. Row n corresponds to the clade starting with the branch above n, the columns correspond to the discrete times of τ_s
    P_clade = zeros(nv(tree),N_τ)
    remaining = trues(nv(tree)) # Track unprocessed nodes
    
    # compute the P(ν|τ_s) for leaves (eq. (3))
    for leaf in leaves
        m = get_weight(tree, leaf, neighbors(tree, leaf)[1])
        P_clade[leaf,:] .= P_branch_only[mp1_to_ms[m+1], :, 1] .* ρ   # factor of ρ for extant individuals to be sampled
        remaining[leaf] = false
    end

    # Compute the P(ν|τ_s) for all nodes for which both descendant clades have already been computed.
    # End when the entire tree up to the two descendants of the root has been computed
    while all(remaining[root_desc].==false)==false
        for node in 1:nv(tree)
            if remaining[node] == true && node != root
                # Collect the three neighbors of the interior node
                nb_computed = Array{Int64}(undef,0)
                nb_remaining = Array{Int64}(undef,0)
                for i in neighbors(tree, node)
                    if remaining[i]
                        push!(nb_remaining, i)
                    else
                        push!(nb_computed, i)
                    end
                end
                # Compute P(ν|τ_s) if two of the neighbors have already been covered. The third remaining neighbor is the predecessor of the current node
                if length(nb_computed) == 2
                    p_interior_branch!(view(P_clade, node,:), get_weight(tree, node, nb_remaining[1]), P_clade[nb_computed[1],:], P_clade[nb_computed[2],:])
                    remaining[node] = false
                end
            end
        end
    end

    # Integrate over the two clades descending from the root (eq. (5))
    # Take the logarithm and renormalize to 1
    return log(2*β * integral(P_clade[root_desc[1],:] .* P_clade[root_desc[2],:], dτ)) - (nv(tree)-1)*log(norm)
end


"""
Infer most likely parameters of a given tree using Nelder-Mead optimization.

# Arguments:
- `tree`: Phylogenetic tree.
- `ρ`: Sampling probability (default is 1.0).
- `q_0`: Initial death rate (default is random).
- `μ_0`: Initial mutation rate (default is random).
- `ρ_0`: Initial sampling probability (default is random).
- `return_likelihood`: Return likelihood alongside parameters (default is false).
- `show_trace`: Show optimization trace (default is false).
- `infer_ρ`: Whether to infer sampling probability `ρ` (default is false).
- `norm`: Initial normalization factor (default is 20.0).
- `L_tol`: Tolerance for likelihood convergence (default is 1e-8).
- `kwargs`: Additional arguments passed to likelihood computation.

# Returns:
- The inferred parameters (and optionally the likelihood).
"""
function infer_parameters(tree::SimpleWeightedGraph{T, T}, ρ::Real=1.0; q_0::Real=rand(), μ_0::Real=2*rand(), ρ_0::Real=rand(), return_likelihood::Bool=false, show_trace::Bool=false, infer_ρ::Bool=false, norm::Real=20.0, L_tol::Real=1e-8, kwargs...) where T<:Integer
    # Define negative log-likelihood function for the optimizer
    function log_likelihood(params::Vector{T}) where T<:Real
        # Unpack the parameters
        if infer_ρ
            q_, μ_, ρ_ = params
        else
            q_, μ_ = params
            ρ_ = ρ
        end
        # Return a high value if any parameters are negative or either q or ρ > 1
        if any(params.<0.0) || any([q_, ρ_].>1.0)
            return 10000.0 - sum((params.<0.0).*params)*1000.0 + sum(([q_, ρ_].>1.0).*[q_, ρ_])*1000
        else
        
            # Compute the log-likelihood
            L = compute_param_likelihood(tree, q_, μ_; ρ=ρ_, norm=norm, kwargs...)

            # Check if `norm` needs to be adapted. If so, repeat with the new value.
            while isnan(L) || abs(L)==Inf
                if L==-Inf
                    println("L=$(L), choose larger normalization!")
                    norm *= 1.08
                else
                    println("L=$(L), choose smaller normalization!")
                    norm *= 0.9
                end
                println("norm changed to $(norm)")
                L = compute_param_likelihood(tree, q_, μ_; ρ=ρ_, norm=norm, kwargs...)
            end
            return - L
        end
    end

    # Initialize parameters
    params_0 = [q_0, μ_0]
    if infer_ρ
        push!(params_0, ρ_0)
    end
    
    # Optimize using Nelder-Mead
    result = optimize(params -> log_likelihood(params), params_0, NelderMead(), Optim.Options(g_tol=L_tol, trace_simplex=show_trace, show_trace=show_trace))
    params_inf = Optim.minimizer(result)
    
    # Return the parameters of the maximum (and the log-likelihood value, if desired)
    if return_likelihood
        return params_inf, - Optim.minimum(result)
    else
        return params_inf
    end
end