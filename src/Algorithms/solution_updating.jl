# Updating functions for the solution

function is_path_partial(TTGraph::TravelTimeGraph, bundle::Bundle, path::Vector{Int};)
    bundle.supplier != TTGraph.networkNodes[path[1]] && return true
    bundle.customer != TTGraph.networkNodes[path[end]] && return true
    return false
end

function add_bundle!(
    solution::Solution,
    instance::Instance,
    bundle::Bundle,
    path::Vector{Int};
    sorted::Bool=false,
    skipFill::Bool=false,
    verbose::Bool=false,
)
    verbose && println("Adding bundle $(bundle) on path $path")
    # If nothing to do, returns nothing
    length(path) == 0 && return 0.0
    TSGraph, TTGraph = instance.timeSpaceGraph, instance.travelTimeGraph
    # Adding the bundle to the solution
    remove_shortcuts!(path, TTGraph)
    verbose && println("Shortcut removed : $path")
    verbose && println("Is path partial : $(is_path_partial(TTGraph, bundle, path))")
    verbose && println("Path nodes : $(TTGraph.networkNodes[path])")
    add_path!(solution, bundle, path; partial=is_path_partial(TTGraph, bundle, path))
    # Updating the bins
    verbose && println("Skip fill : $skipFill")
    skipFill && return 0.0
    verbose && println("Updating bins on path $path for bundle $bundle")
    updateCost = update_bins!(
        solution, TSGraph, TTGraph, bundle, path; sorted=sorted, verbose=verbose
    )
    verbose && println("Cost added for bundle : $updateCost")
    return updateCost
end

# Combine all bundles paths in arguments into a sparse matrix indicating the arcs to work with
function get_bins_updated(
    TSGraph::TimeSpaceGraph,
    TTGraph::TravelTimeGraph,
    bundles::Vector{Bundle},
    paths::Vector{Vector{Int}},
)
    I, J = Int[], Int[]
    # For every bundle path and every order in the bundle, adding the timed nodes in the matrix indices
    for (bundle, path) in zip(bundles, paths)
        for order in bundle.orders
            timedPath = time_space_projector(TTGraph, TSGraph, path, order)
            # Without checking overlapping as the combine function will take care of it
            append!(I, timedPath[1:(end-1)])
            append!(J, timedPath[2:end])
        end
    end
    # Garbage colecting is here
    V = ones(Bool, length(I))
    # Combine function for bools is | by default
    return sparse(I, J, V)
end

function refill_bins!(
    bins::Vector{Bin}, fullCapacity::Int, ALL_COMMODITIES::Vector{Commodity}
)
    # Bound filtering on the recomputation : if the bins already attain the lower bound, no need to optimize the storage
    length(bins) == 0 && return 0
    binsBefore = length(bins)
    if sum(bin.load for bin in bins) == 0
        empty!(bins)
        return -binsBefore
    end
    ceil(sum(bin.load for bin in bins) / fullCapacity) == length(bins) && return 0
    allCommodities = get_all_commodities(bins, ALL_COMMODITIES)
    empty!(bins)
    # Filling it back again
    first_fit_decreasing!(bins, fullCapacity, allCommodities; sorted=false)
    return length(bins) - binsBefore
end

# Refill bins on the working arcs, to be used after bundle removal
function refill_bins!(
    solution::Solution,
    TSGraph::TimeSpaceGraph,
    workingArcs::SparseMatrixCSC{Bool,Int},
    ALL_COMMODITIES::Vector{Commodity},
)
    costAdded = 0.0
    # Efficient iteration over sparse matrices
    rows = rowvals(workingArcs)
    for tDst in 1:Base.size(workingArcs, 2)
        for srcIdx in nzrange(workingArcs, tDst)
            tSrc = rows[srcIdx]
            arcData = TSGraph.networkArcs[tSrc, tDst]
            # No need to refill bins on linear arcs
            arcData.isLinear && continue
            # Adding new bins cost
            addedBins = refill_bins!(
                solution.bins[tSrc, tDst], arcData.capacity, ALL_COMMODITIES
            )
            costAdded += addedBins * arcData.unitCost
        end
    end
    return costAdded
end

# Refill bins on the path given, to be used after single bundle removal
function refill_bins!(
    solution::Solution,
    TTGraph::TravelTimeGraph,
    TSGraph::TimeSpaceGraph,
    bundle::Bundle,
    path::Vector{Int},
    ALL_COMMODITIES::Vector{Commodity};
    verbose::Bool=false,
)
    costAdded = 0.0
    verbose && println("Refilling bins on path $path for bundle $bundle")
    # Projecting path for every order
    for order in bundle.orders
        timedPath = time_space_projector(TTGraph, TSGraph, path, order)
        # Refilling all bins on this path
        for (tSrc, tDst) in partition(timedPath, 2, 1)
            arcData = TSGraph.networkArcs[tSrc, tDst]
            verbose && println("Arc $tSrc -> $tDst : $(arcData)")
            # No need to refill bins on linear arcs
            arcData.isLinear && continue
            # Adding new bins cost
            arcBins = solution.bins[tSrc, tDst]
            verbose &&
                println("Bins : $(length(arcBins)) -> $([bin.load for bin in arcBins]) m3")
            addedBins = refill_bins!(arcBins, arcData.capacity, ALL_COMMODITIES)
            verbose && println("Refilling")
            verbose &&
                println("Bins : $(length(arcBins)) -> $([bin.load for bin in arcBins]) m3")
            verbose && println("Bins added : $addedBins")
            costAdded += addedBins * arcData.unitCost
            verbose && println("Total cost added : $costAdded")
        end
    end
    return costAdded
end

# Remove the bundle only on the path portion provided
function remove_bundle!(
    solution::Solution,
    instance::Instance,
    bundle::Bundle,
    path::Vector{Int};
    verbose::Bool=false,
)
    TSGraph, TTGraph = instance.timeSpaceGraph, instance.travelTimeGraph
    oldPart = Int[]
    if length(path) == 0 || !is_path_partial(TTGraph, bundle, path)
        verbose && println(
            "\nRemoving full path for bundle $bundle (path = $path and sol path = $(solution.bundlePaths[bundle.idx]))",
        )
        oldPart = remove_path!(solution, bundle)
    else
        verbose && println(
            "\nRemoving partial path for bundle $bundle (path = $path and sol path = $(solution.bundlePaths[bundle.idx]))",
        )
        oldPart = remove_path!(solution, bundle; src=path[1], dst=path[end])
    end
    verbose && println("Old part : $oldPart")
    costRemoved = update_bins!(
        solution, TSGraph, TTGraph, bundle, oldPart; remove=true, verbose=verbose
    )
    verbose && println("\nCost removed (without refilling) : $costRemoved")
    return (costRemoved, oldPart)
end

# Update the current solution for several bundles
# Providing a path with remove option means to remove the bundle only on the path portion provided
function update_solution!(
    solution::Solution,
    instance::Instance,
    bundles::Vector{Bundle},
    paths::Vector{Vector{Int}}=[Int[] for _ in 1:length(bundles)],
    ALL_COMMODITIES::Vector{Commodity}=Commodity[];
    remove::Bool=false,
    sorted::Bool=false,
    skipRefill::Bool=false,
)
    costAdded = 0.0
    if !remove
        # If remove = false, adding the bundle to the solution
        for (bundle, path) in zip(bundles, paths)
            costAdded += add_bundle!(
                solution, instance, bundle, path; sorted=sorted, skipFill=skipRefill
            )
        end
    else
        pathsToUpdate = Vector{Vector{Int}}()
        # If remove = true, removing the bundle from the solution
        for (bundle, path) in zip(bundles, paths)
            costRemoved, oldPart = remove_bundle!(solution, instance, bundle, path)
            costAdded += costRemoved
            push!(pathsToUpdate, oldPart)
        end
        # If skipRefill than no recomputation
        skipRefill && return costAdded
        # Than refilling the bins
        binsUpdated = get_bins_updated(
            instance.timeSpaceGraph, instance.travelTimeGraph, bundles, pathsToUpdate
        )
        costAdded += refill_bins!(
            solution, instance.timeSpaceGraph, binsUpdated, ALL_COMMODITIES
        )
    end
    return costAdded
end

# Update the current solution (optimized for a unique bundle)
# Providing a path with remove option means to remove the bundle only on the path portion provided
function update_solution!(
    solution::Solution,
    instance::Instance,
    bundle::Bundle,
    path::Vector{Int}=Int[],
    ALL_COMMODITIES::Vector{Commodity}=Commodity[];
    remove::Bool=false,
    sorted::Bool=false,
    skipRefill::Bool=false,
    verbose::Bool=false,
)
    # If remove = false, adding the bundle to the solution
    !remove && return add_bundle!(
        solution,
        instance,
        bundle,
        path;
        sorted=sorted,
        skipFill=skipRefill,
        verbose=verbose,
    )
    verbose && println("\nRemoving bundle $bundle of path $path")
    # Otherwise, removing it from the solution
    costRemoved, oldPart = remove_bundle!(solution, instance, bundle, path; verbose=verbose)
    verbose && println("\nCost removed (without refilling) : $costRemoved")
    # If skipRefill than no recomputation
    verbose && println("Skip refill : $skipRefill")
    skipRefill && return costRemoved
    # Than refilling the bins
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    verbose && println("Refilling bins on old part $oldPart")
    costRemoved += refill_bins!(
        solution, TTGraph, TSGraph, bundle, oldPart, ALL_COMMODITIES; verbose=verbose
    )
    verbose && println("\nCost removed (with refilling) : $costRemoved")
    return costRemoved
end

# Removing all empty bins from the linear arcs (to be used before extraction)
function clean_empty_bins!(solution::Solution, instance::Instance)
    TSGraph = instance.timeSpaceGraph
    for (i, arc) in enumerate(edges(TSGraph.graph))
        filter!(bin -> bin.load > 0, solution.bins[src(arc), dst(arc)])
        i % 1000 == 0 && print("|")
    end
    return println()
end
