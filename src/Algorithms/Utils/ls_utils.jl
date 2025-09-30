# Utils function for local search neighborhoods

function is_bin_candidate(bins::Vector{Bin}, arcData::NetworkArc; skipLinear::Bool)
    # If there is no bins, one bin or the arc is linear, skipping arc
    length(bins) <= 1 && return false
    skipLinear && arcData.isLinear && return false
    # If there is no gap with the lower bound, skipping arc
    arcVolume = sum(bin.load for bin in bins)
    ceil(arcVolume / arcData.capacity) == length(bins) && return false
    return true
end

function compute_new_bins(
    arcData::NetworkArc, allCommodities::Vector{Commodity}; sorted::Bool
)
    newBins = first_fit_decreasing(Bin[], arcData.capacity, allCommodities; sorted=sorted)
    bfdBins = best_fit_decreasing(Bin[], arcData.capacity, allCommodities; sorted=sorted)
    length(newBins) > length(bfdBins) && (newBins = bfdBins)
    return newBins
end

function compute_new_bins(
    arcData::NetworkArc,
    allCommodities::SubArray{
        Commodity,1,Vector{Commodity},Tuple{UnitRange{Int64}},true
    };
    sorted::Bool,
)
    newBins = first_fit_decreasing(Bin[], arcData.capacity, allCommodities; sorted=sorted)
    bfdBins = best_fit_decreasing(Bin[], arcData.capacity, allCommodities; sorted=sorted)
    length(newBins) > length(bfdBins) && (newBins = bfdBins)
    return newBins
end

function tentative_first_fit(
    arcData::NetworkArc,
    commodities::SubArray{
        Commodity,1,Vector{Commodity},Tuple{UnitRange{Int64}},true
    },
    CAPACITIES::Vector{Int},
)
    return tentative_first_fit(
        Bin[], arcData.capacity, commodities, CAPACITIES; sorted=false
    )
end

function tentative_best_fit(
    arcData::NetworkArc,
    commodities::SubArray{
        Commodity,1,Vector{Commodity},Tuple{UnitRange{Int64}},true
    },
    CAPACITIES::Vector{Int},
)
    capa = arcData.capacity
    return tentative_best_fit([Bin(capa)], capa, commodities, CAPACITIES; sorted=false) + 1
end

# Store previous bins before removing commodities from them
function save_previous_bins(solution::Solution, workingArcs::SparseMatrixCSC{Bool,Int})
    I, J, _ = findnz(workingArcs)
    oldBins = Vector{Vector{Bin}}(undef, length(I))
    # Efficient iteration over sparse matrices
    rows = rowvals(workingArcs)
    for timedDst in 1:Base.size(workingArcs, 2)
        for srcIdx in nzrange(workingArcs, timedDst)
            timedSrc = rows[srcIdx]
            # Storing old bins
            oldBins[srcIdx] = my_deepcopy(solution.bins[timedSrc, timedDst])
        end
    end
    return sparse(I, J, oldBins)
end

# Revert the bin loading the the vector of bins given
function revert_bins!(solution::Solution, previousBins::SparseMatrixCSC{Vector{Bin},Int})
    # Efficient iteration over sparse matrices
    rows = rowvals(previousBins)
    for timedDst in 1:Base.size(previousBins, 2)
        for srcIdx in nzrange(previousBins, timedDst)
            timedSrc = rows[srcIdx]
            # Reverting to previous bins
            solution.bins[timedSrc, timedDst] = previousBins[timedSrc, timedDst]
        end
    end
end

# Checking if two nodes are candidates for two node neighborhhods
function are_nodes_candidate(TTGraph::TravelTimeGraph, src::Int, dst::Int)
    src == dst && return false
    # Cannot go back in time
    TTGraph.stepToDel[src] < TTGraph.stepToDel[dst] && return false
    # Not authorizing going to the same node at different time steps
    TTGraph.networkNodes[src] == TTGraph.networkNodes[dst] && return false
    # Need to have existing path
    return has_path(TTGraph.graph, src, dst)
end

# Checks whether the node is before node in the path
function is_node1_before_node2(path::Vector{Int}, node1::Int, node2::Int)
    # look for the first node encountered
    idx = findfirst(node -> (node == node1) || (node == node2), path)
    return path[idx] == node1
end

# Return a vector of bundles to update 
# If node 1 and node 2 are given : bundles that flow from 1 to 2
# If only node 1 : bundle that have node 1 for destination
function get_bundles_to_update(
    TTGraph::TravelTimeGraph, solution::Solution, node1::Int, node2::Int=-1
)
    if node2 == -1
        if node1 == -1
            return Int[]
        else
            return solution.bundlesOnNode[node1]
        end
    end
    # Adaptation for split bundles (we will want to do a two node neighborhood with all bundles from one supplier to one plant)
    if TTGraph.networkNodes[node1].type == :supplier
        # Keeping only the bundles flowing from one supplier to one plant
        return filter(b -> TTGraph.bundleSrc[b] == node1, solution.bundlesOnNode[node2])
    end
    # Classical case
    twoNodeBundleIdxs = intersect(
        get(solution.bundlesOnNode, node1, Int[]), get(solution.bundlesOnNode, node2, Int[])
    )
    return filter(
        b -> is_node1_before_node2(solution.bundlePaths[b], node1, node2), twoNodeBundleIdxs
    )
end

# Return a vector of path portion to update for bundles that flow from 1 to 2
function get_paths_to_update(
    solution::Solution, bundles::Vector{Bundle}, node1::Int, node2::Int
)
    paths = Vector{Vector{Int}}(undef, length(bundles))
    for (idx, bundle) in enumerate(bundles)
        oldPath = solution.bundlePaths[bundle.idx]
        srcIdx = findfirst(node -> node == node1, oldPath)
        dstIdx = findlast(node -> node == node2, oldPath)
        paths[idx] = oldPath[srcIdx:dstIdx]
    end
    return paths
end

# computes an (over)estimated number of bins removed 
function estimated_transport_units(order::Order, bins::Vector{Bin})
    n, vol = 0, 0
    for binIdx in sortperm(bins; by=bin -> bin.load)
        vol += bins[binIdx].load
        vol > order.volume && return n
        n += 1
    end
    return n
end

# computes an estimated removal cost of a bundle
function bundle_estimated_removal_cost(
    bundle::Bundle, path::Vector{Int}, instance::Instance, solution::Solution
)
    cost = 0.0
    TTGraph, TSGraph = instance.travelTimeGraph, instance.timeSpaceGraph
    for (i, j) in partition(path, 2, 1), order in bundle.orders
        cost += volume_stock_cost(TTGraph, i, j, order)
        arcData = TTGraph.networkArcs[i, j]
        ti, tj = time_space_projector(TTGraph, TSGraph, i, j, order)
        transportUnits = if arcData.isLinear
            get_transport_units(order, arcData)
        else
            estimated_transport_units(order, solution.bins[ti, tj])
        end
        cost += transportUnits * arcData.unitCost
    end
    return cost
end

# Fusing all bundles into one having the same supplier and customer as the first one
function merge_bundles(instance::Instance, bundles::Vector{Bundle}, CAPACITIES::Vector{Int})
    # Putting one order for each delivery date to fuse them together
    newOrders = [Order(UInt(0), i) for i in 1:(instance.timeHorizon)]
    for bundle in bundles
        # Fusing all orders
        for order in bundle.orders
            append!(newOrders[order.deliveryDate].content, order.content)
        end
    end
    filter!(o -> length(o.content) > 0, newOrders)
    for order in newOrders
        sort!(order.content; rev=true)
    end
    newOrders = [
        add_properties(order, tentative_first_fit, CAPACITIES) for order in newOrders
    ]
    supp, cust = bundles[1].supplier, bundles[1].customer
    maxDelTime, maxDelIdx = findmax(b -> b.maxDelTime, bundles)
    maxPackSize = maximum(b -> b.maxPackSize, bundles)
    fusedBunIdx = bundles[maxDelIdx].idx
    return Bundle(supp, cust, newOrders, fusedBunIdx, UInt(0), maxPackSize, maxDelTime)
end

# Increase direct costs a hundred fold and reintroduce directs 
function forbid_directs!(instance::Instance)
    # Multiplying by a constant factor the number of trucks precomputed for directs
    for bundle in instance.bundles
        for order in bundle.orders
            order.bpUnits[:direct] *= 100
        end
    end
end

# Revert direct cost to normal values and reintroduce all bundles
function allow_directs!(instance::Instance)
    for bundle in instance.bundles
        for order in bundle.orders
            order.bpUnits[:direct] /= 100
        end
    end
end

