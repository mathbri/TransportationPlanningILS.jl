# Perturbation object to define the neighborhoods used for milps in LNS

struct Perturbation
    # Type of neighborhood 
    type::Symbol
    # Bundle idxs concerned
    bundleIdxs::Vector{Int}
    # Nodes used for the source and destination of the bundles for two shared node neighborhoods 
    src::Int
    dst::Int
    # New paths to use for attract and reduce neighborhoods
    oldPaths::Vector{Vector{Int}}
    newPaths::Vector{Vector{Int}}
    # Transport units completion through time 
    loads::SparseMatrixCSC{Int,Int}
end

function Perturbation(
    type::Symbol,
    bundleIdxs::Vector{Int},
    oldPaths::Vector{Vector{Int}},
    loads::SparseMatrixCSC{Int,Int},
)
    return Perturbation(type, bundleIdxs, 0, 0, oldPaths, Vector{Int}[], loads)
end

function Perturbation(
    type::Symbol,
    bundleIdxs::Vector{Int},
    oldPaths::Vector{Vector{Int}},
    newPaths::Vector{Vector{Int}},
    loads::SparseMatrixCSC{Int,Int},
)
    return Perturbation(type, bundleIdxs, 0, 0, oldPaths, newPaths, loads)
end

function Perturbation(
    type::Symbol,
    bundleIdxs::Vector{Int},
    src::Int,
    dst::Int,
    oldPaths::Vector{Vector{Int}},
    loads::SparseMatrixCSC{Int,Int},
)
    return Perturbation(type, bundleIdxs, src, dst, oldPaths, Vector{Int}[], loads)
end

function is_attract_reduce(perturbation::Perturbation)
    return perturbation.type == :attract_reduce
end

function is_two_shared_node(perturbation::Perturbation)
    return perturbation.type == :two_shared_node
end

function number_of_variables(perturbation::Perturbation, instance::Instance)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    nPathVar = if is_attract_reduce(perturbation)
        length(instance.bundles)
    else
        sum(b -> length(TTGraph.bundleArcs[b]), perturbation.bundleIdxs)
    end
    nPackVar = length(TSGraph.commonArcs)
    return nPathVar + nPackVar
end

function perturbation_variables(type::Symbol, bundleIdxs::Vector{Int}, instance::Instance)
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    nPathVar = if type == :attract_reduce
        length(instance.bundles)
    else
        sum(b -> length(TTGraph.bundleArcs[b]), bundleIdxs)
    end
    nPackVar = length(TSGraph.commonArcs)
    return nPathVar + nPackVar
end

function lower_bound_milp_variables(instance::Instance)
    return perturbation_variables(:single_plant, idx(instance.bundles), instance)
end