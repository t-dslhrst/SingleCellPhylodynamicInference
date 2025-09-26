using Distributions
using SparseArrays
using Graphs, SimpleWeightedGraphs

# ChatGPT was used to improve docstrings

"""
Performs a single time step in a birth-death process.

# Arguments:
- `living`: 2D matrix of alive individuals. Each row corresponds to an individual, with the columns representing:
    1. Individual's ID
    2. Parent's ID
    3. Number of mutations since the origin of the population
    4. (Optional) Last time step the individual was alive
- `dead`: 2D matrix of dead individuals, with the same format as `living`.
- `β`: Birth rate.
- `δ`: Death rate.
- `μ`: Mutation parameter.
- `dt`: Time step size.
- `birth_weights`: If `true`, one mutation occurs in every birth event.
- `mol_clock`: If `true`, mutations are generated at rate `μ`. Otherwise, mutations follow a Poisson distribution with mean `μ`.

# Optional Arguments:
- `last_born`: The maximum number used for labeling individuals. Defaults to 0, which computes it from the matrices.
    
# Returns:
- Updated `living` matrix after the time step.
- Updated `dead` matrix after the time step.
- Number of births that occurred.
- Number of deaths that occurred.
- Number of mutations that occurred.

"""
function timestep(living::Matrix{T}, dead::Matrix{T}, β::Real, δ::Real, μ::Real, dt::Real, birth_weights::Bool, mol_clock::Bool; last_born::Integer=0) where T <: Integer
    # Initialize last_born if not provided
    if last_born==0
        last_born = maximum([maximum(living[:,1]),maximum(dead[:,1], init=0)]) # number of the most recently born individual
    end
    # Initialize new matrix for individuals still alive after this time step
    living_new = zeros(T, 2*length(living[:,1]),length(living[1,:]))

    # initialize counters
    births = 0
    deaths = 0
    total_mutations = 0
    living_new_counter = 0
    dead_counter = length(dead[:,1])

    # Expand dead matrix to accommodate newly dead individuals
    dead = vcat(dead, zeros(Int32, length(living[:,1]),length(living[1,:])))

    # Iterate over all living individuals and determine their fate
    for i in 1:length(living[:,1])
        x = rand()
        # Death
        if x < dt*δ
            dead_counter += 1
            dead[dead_counter,:] = living[i,:]
            deaths += 1

        # Birth event: Two offspring (with mutations) are created, the parent is considered dead
        elseif x < dt*(β+δ)
            # add parent to dead individuals
            dead_counter += 1
            dead[dead_counter,:] = living[i,:]
            # create offspring
            for _ in 1:2
                last_born += 1
                if birth_weights
                    mutations = 1
                elseif mol_clock
                    mutations = Int64(rand() < μ*dt)
                else
                    mutations = rand(Poisson(μ))
                end
                total_mutations += mutations

                if length(living[1,:])==4
                    living_new_counter += 1
                    living_new[living_new_counter,:] .= [last_born, living[i,1], living[i,3]+mutations, living[i,4]]
                else
                    living_new_counter += 1
                    living_new[living_new_counter,:] = [last_born, living[i,1], living[i,3]+mutations]
                end
            end
            births += 1
        
        # If neither birth nor death, individual survives and is added to the living_new matrix
        else
            if mol_clock
                # Nn case of a molecular clock, add mutations with rate μ
                mutation = Int64(rand() < μ*dt)
                if length(living[1,:])==4
                    living_new_counter += 1
                    living_new[living_new_counter,:] = [living[i,1], living[i,2], living[i,3]+mutation, living[i,4]]
                else
                    living_new_counter += 1
                    living_new[living_new_counter,:] = [living[i,1], living[i,2], living[i,3]+mutation]
                end
            else
                living_new_counter += 1
                living_new[living_new_counter,:] = living[i,:]
            end
        end
    end
    return living_new[1:living_new_counter,:], dead[1:dead_counter,:], births, deaths, total_mutations
end


"""
Grows a population to a specified size (`N_final`) or time (`T_final`), whichever is reached first.
It also computes empirical values of q (deaths/births) and μ (mutations/birth events).

# Arguments:
- `N_final`: Final population size.
- `β`: Birth rate.
- `δ`: Death rate.
- `μ`: Mutation parameter.
- `dt`: Time step size.

# Optional Arguments:
- `T_final`: Final time. Defaults to infinity.
- `birth_weights`: If `true`, one mutation occurs per birth event.
- `timestep_weights`: If `true`, track the last time step until which individuals are alive.
- `return_T`: If `true`, return the final time `T`.
- `mol_clock`: If `true`, use `μ` as mutation rate. Otherwise, mutations occur at birth events following a Poisson distribution with mean `μ`.
- `IntType`: Data type for individual IDs and mutations. Defaults to `Int32`.

# Returns:
- `living`: 2D matrix of alive individuals. Each row corresponds to an individual, with the columns representing:
    1. Individual's ID
    2. Parent's ID
    3. Number of mutations since the origin of the population
    4. (Optional) Last time step the individual was alive
- `dead`: 2D matrix of dead individuals, with the same format as `living`.
- Empirical death rate (q).
- Empirical mutation rate (μ).
- (Optional) Final time `T` if `return_T` is `true`.

"""
function grow_population(N_final::Integer,  β::Real, δ::Real, μ::Real, dt::Real; T_final::Real=Inf, birth_weights::Bool=false, timestep_weights::Bool=false, return_T::Bool=false, mol_clock::Bool=false, IntType::DataType=Int32)
    # Initialize matrices for alive and dead individuals
    if timestep_weights
        living = ones(IntType,1,4)
        living[1,4] = 0
        dead = Matrix{IntType}(undef, 0,4)
    else 
        living = ones(IntType,1,3)
        dead = Matrix{IntType}(undef, 0,3)
    end
    # set the number of mutations of the first individual to 0
    living[1,3] = 0

    # initialize some counters
    births = 0
    deaths = 0
    mutations = 0
    T = 0

    # Perform time steps until population size or time limit is reached
    while length(living[:,1]) < N_final && T < T_final
        # Update last time step for living individuals if tracking time
        if timestep_weights
            living[:,4] .+= 1
        end

        # Perform a time step
        living, dead, new_births, new_deaths, new_mutations = timestep(living, dead, β, δ, μ, dt, birth_weights, mol_clock)

        # Restart if the population has died out
        if length(living[:,1]) == 0
            if timestep_weights
                living = ones(IntType,1,4)
                living[1,4] = 0
                dead = Matrix{IntType}(undef, 0,4)
            else 
                living = ones(IntType,1,3)
                dead = Matrix{IntType}(undef, 0,3)
            end
            living[1,3] = 0
            births = 0
            deaths = 0
            mutations = 0
            T = 0
        else
            # Update counters
            births += new_births
            deaths += new_deaths
            mutations += new_mutations
            if length(living[:,1]) == 1
                T = 0
            else
                T += dt
            end
        end
    end
    # Return population matrices and empirical values.
    # If desired, also return the age T of the population.
    if return_T
        return living, dead, deaths/births, mutations/(2*births), T
    else
        return living, dead, deaths/births, mutations/(2*births)
    end
end


"""
Builds a phylogenetic tree (as `SimpleWeightedGraph`) from a population matrix, 
with edge weights representing the number of mutations between individuals.

# Arguments:
- `complete_population`: 2D matrix of individuals (both alive and dead). The columns represent:
    1. Individual's ID
    2. Parent's ID
    3. Number of mutations

# Returns:
- A `SimpleWeightedGraph` where nodes represent individuals and edges represent parent-child relationships. 
  Edge weights represent the number of mutations.

"""
function build_tree(complete_population::Matrix{T}) where T <: Integer
    # Compute mutations between each individual and its parent
    mutation_diffs = complete_population[:,3]
    for i in 1:length(complete_population[:,1])
        mutation_diffs[i] -= complete_population[findfirst(y->y==complete_population[i,2], complete_population[:,1]),3]
    end

    # Create a graph with edges weighted by mutation differences
    graph = SimpleWeightedGraph(complete_population[:,1], complete_population[:,2], mutation_diffs)
    rem_edge!(graph, 1, 1) # Remove self-loop on the root node
    return graph
end



"""
Reduces a complete tree to only include the extant individuals (leaves).
This reduces the tree to the phylogenetic relationships between living individuals.

# Arguments:
- `graph_`: The full `SimpleWeightedGraph` containing all individuals.
- `leaves_`: Array of IDs of extant individuals (leaves).

# Returns:
- A reduced `SimpleWeightedGraph` containing only the extant individuals and their ancestors.
- Updated array of `leaves` corresponding to the new graph.

"""
function reduce_tree(graph_::SimpleWeightedGraph{T, T}, leaves_::Array{T}) where T <: Integer
    # Unfortunately, new edges with weight 0 cannot be added.
    # Copy graph and convert weights to Float to allow adding small weights (instead of 0)
    graph = SimpleWeightedGraph(convert(SparseMatrixCSC{Float32, T}, copy(graph_.weights)))
    leaves = copy(leaves_)

    # Remove all "dead ends" (degree=1 and not a leaf)
    i = 2
    while i <= nv(graph)
        if !(i in leaves)
            nb = neighbors(graph,i)
            if length(nb)==1
                i_temp = i
                # Continue removing neighbors that become dead ends
                while length(nb)==1
                    # Remove dead end i_temp. The n-th node will become the (n-1)-th node, if n>i_temp
                    rem_vertex!(graph, i_temp)
                    # Descrease the index of all nodes of interest if their number was larger than i_temp
                    i -= (i>=i_temp)
                    leaves .-= leaves .> i_temp
                    # Now look at the neighbor of the former dead end
                    i_temp = nb[1]
                    nb = neighbors(graph,i_temp)
                end
            end
        end
        # Go on
        i += 1
    end

    # Remove intermediate generations (degree=2) and connect their neighbors (add weights together)
    # Start at node 2, as 1 is the root node (which always has two neighbors)
    i = 2
    while i <= nv(graph)
        nb = neighbors(graph,i)
        if length(nb)==2 && i>1
            # Determine the weight of the new edge
            weight_temp = get_weight(graph, nb[1], i) + get_weight(graph, nb[2], i)
            if weight_temp==0.0 # Add a small weight if it is zero
                weight_temp += 1e-9
            end
            # Add an edge between the two neighbors
            add_edge!(graph, nb[1], nb[2], weight_temp)
            # Remove the intermediate node
            rem_vertex!(graph, i)
            # Adjust the indices of the leaves
            leaves .-= leaves .> i     
        else
            # Go on, if not an intermediate generation
            i += 1          
        end
    end

    # Remove root if it no longer splits up directly
    while degree(graph)[1]==1
        rem_vertex!(graph, 1)
        leaves .-= leaves .> 1
    end

    # Create a new tree with integer weights
    new_weights = graph.weights
    for (i, weight) in enumerate(new_weights)
        if weight < 1
            if weight > 0
                new_weights[i] = 0.0
            end
        end
    end
    new_weights[new_weights .>= 1] .= floor.(new_weights[new_weights .>= 1])
    
    return SimpleWeightedGraph(convert(SparseMatrixCSC{T, T}, new_weights)), leaves
end



"""
Finds the direct descendants (children) of a given node in a rooted tree.

# Arguments:
- `tree`: A `SimpleWeightedGraph` representing the phylogenetic tree.
- `node`: The node (vertex) whose descendants are to be found.

# Returns:
- A vector of node IDs representing the descendants of the given node.
  Assumes that descendants have higher node IDs than their parents.
"""
function find_descendants(tree, node::Integer)
    nb = neighbors(tree, node)
    return nb[nb .> node]
end


"""
Finds all nodes in the subtree rooted at a given node.

# Arguments:
- `tree`: A `SimpleWeightedGraph` representing the phylogenetic tree.
- `root_subtree`: The node ID of the root of the desired subtree.

# Returns:
- A vector of node IDs in the subtree rooted at `root_subtree`, including the root itself.
"""
function find_subtree(tree, root_subtree)
    vlist = Vector{Int64}(undef,0)
    add = [root_subtree]
    while length(add)>0
        node = add[1]
        add = vcat(add, find_descendants(tree,node))
        push!(vlist, add[1])
        popfirst!(add)
    end
    return vlist
end


"""
Finds the Most Recent Common Ancestor (MRCA) of a set of nodes in a rooted tree.

# Arguments:
- `tree`: A `SimpleWeightedGraph` representing the phylogenetic tree.
- `node_list`: A vector of node IDs for which the MRCA is to be found.

# Returns:
- The node ID of the MRCA of the nodes in `node_list`.
"""
function find_MRCA(tree::SimpleWeightedGraph{T, T}, node_list) where T <: Integer
    # choose first node of the list, set it as preliminary MRCA and find its path to the root
    node = node_list[1]
    MRCA = node
    path_1 = [node]
    while node > 1
        node = neighbors(tree, node)[1]
        push!(path_1, node)
    end
    # iterate over all other nodes in the list and find the first intersection of their path to the root with the path_1
    for i in node_list[2:end]
        node = i
        while !(node in path_1)
            node = neighbors(tree, node)[1]
        end
        if node < MRCA
            MRCA = node
        end
    end
    return MRCA
end


"""
Simulates allelic dropout in a phylogenetic tree based on a probabilistic model.

Each mutation in each sample can be lost with probability `σ`. If lost, the mutation
can either disappear (if not detected in any sample) or be reassigned tor the branch above the MRCA
of all samples where the mutation is still detected.
If a mutation is only detected in one sample, it is assigned to the pendant branch of that sample.

# Arguments:
- `tree`: A `SimpleWeightedGraph` representing the original phylogenetic tree.
- `σ`: Dropout probability (between 0 and 1).

# Returns:
- A new `SimpleWeightedGraph` with updated mutation weights after allelic dropout.
"""
function simulate_allelic_dropout(tree::SimpleWeightedGraph{T, T}, σ::Real) where T <: Integer
    # initialize tree with weights of zero
    weights0 = convert(SparseMatrixCSC{Float16, T}, copy(tree.weights))
    weights0[weights0 .> 0] .= 0.1
    weights0[weights0 .>= 0] .= floor.(weights0[weights0 .>= 0])
    tree_post_ad = SimpleWeightedGraph(convert(SparseMatrixCSC{T, T}, weights0))
    weights0 = nothing
    # iterate over all branches and all mutations along them
    for node in 2:nv(tree)
        # identify predecessor and descendants of the node
        predecessor = neighbors(tree, node)[1]
        node_is_sample = degree(tree)[node] == 1
        mutations_along_branch = tree.weights[node, predecessor]
        if node_is_sample
            for _ in 1:mutations_along_branch
                # add mutation to the pendant branch with probability 1-σ
                if rand() < 1-σ
                    tree_post_ad.weights[node, predecessor] += 1
                    tree_post_ad.weights[predecessor, node] += 1
                end
            end
        else
            # make list of descendants that inherit the mutation
            descendants = find_descendants(tree, node)
            for _ in 1:mutations_along_branch
                # make list of descendants where the mutation is detected
                descendants_with_mutation = Vector{T}(undef,0)
                for desc in descendants
                    if rand() < 1-σ
                        push!(descendants_with_mutation, desc)
                    end
                end
                # if only found in one sample, add mutation to the pendant branch
                if length(descendants_with_mutation) == 1
                    tree_post_ad.weights[descendants_with_mutation[1], neighbors(tree, descendants_with_mutation[1])[1]] += 1
                    tree_post_ad.weights[neighbors(tree, descendants_with_mutation[1])[1], descendants_with_mutation[1]] += 1
                # if found in multiple samples, add mutation to the branch above their MRCA
                elseif length(descendants_with_mutation) > 1
                    mrca = find_MRCA(tree, descendants_with_mutation)
                    tree_post_ad.weights[mrca, neighbors(tree, mrca)[1]] += 1
                    tree_post_ad.weights[neighbors(tree, mrca)[1], mrca] += 1
                end
            end
        end
    end
    return tree_post_ad
end