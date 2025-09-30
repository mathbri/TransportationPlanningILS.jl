# Utils function only used in greedy

# With the pre-processing, we only have to check the shortcut type now
# Check whether the arc is fit for a cost update
function is_update_candidate(TTGraph::TravelTimeGraph, src::Int, dst::Int, bundle::Bundle)
    arcData = TTGraph.networkArcs[src, dst]
    # If it is a shortcut leg, cost alredy set to EPS
    arcData.type == :shortcut && return false
    return true
end

# Check whether the arc is forbidden for the bundle
function is_forbidden(TTGraph::TravelTimeGraph, src::Int, dst::Int, bundle::Bundle)
    # If it is an inland bundle, I want to avoid ports
    inlandBundle = (bundle.customer.continent == bundle.supplier.continent)
    return (inlandBundle && (is_port(TTGraph, src) || is_port(TTGraph, dst)))
end

# Computes volume and lead time costs for an order
function volume_stock_cost(
    TTGraph::TravelTimeGraph, src::Int, dst::Int, order::Order; verbose::Bool=false
)
    dstData, arcData = TTGraph.networkNodes[dst], TTGraph.networkArcs[src, dst]
    verbose && println(
        "Dst $dstData cost : $(dstData.volumeCost) * $(order.volume) / $VOLUME_FACTOR = $(dstData.volumeCost * order.volume / VOLUME_FACTOR)",
    )
    verbose && println(
        "Arc $arcData cost : $(arcData.carbonCost) * $(order.volume) / $(arcData.capacity) = $(arcData.carbonCost * order.volume / arcData.capacity)",
    )
    verbose && println("Order : $order $(order.volume) $(order.stockCost)")
    # Node volume cost + Arc carbon cost + Commodity stock cost
    return dstData.volumeCost * order.volume / VOLUME_FACTOR +
           arcData.carbonCost * order.volume / arcData.capacity +
           arcData.distance * order.stockCost
end

# Computes transport units for an order
function transport_units(
    solution::Solution,
    TSGraph::TimeSpaceGraph,
    timedSrc::Int,
    timedDst::Int,
    order::Order,
    CAPACITIES::Vector{Int};
    sorted::Bool,
    use_bins::Bool,
    verbose::Bool=false,
)
    arcData = TSGraph.networkArcs[timedSrc, timedDst]
    # Transport cost 
    orderTrucks = get_transport_units(order, arcData)
    verbose && println("Empty transport units for order : $orderTrucks")
    # If we take into account the current solution
    if use_bins && !arcData.isLinear
        bins = solution.bins[timedSrc, timedDst]
        verbose && println(
            "Bins on $timedSrc -> $timedDst : $(length(bins)) -> $([bin.load for bin in bins]) m3",
        )
        # If the arc is not empty, computing a tentative first fit 
        if length(bins) > 0
            orderTrucks = tentative_first_fit(
                bins, arcData, order, CAPACITIES; sorted=sorted
            )
            verbose && println("Tentative first fit for order : $orderTrucks")
        end
    end
    return orderTrucks
end

# Returns the corresponding transport costs to be used
function transport_cost(
    TSGraph::TimeSpaceGraph, timedSrc::Int, timedDst::Int; current_cost::Bool=false
)
    current_cost && return TSGraph.currentCost[timedSrc, timedDst]
    return TSGraph.networkArcs[timedSrc, timedDst].unitCost
end

# Compute the arc update cost to be used in the path computation
function arc_update_cost(
    sol::Solution,
    TTGraph::TravelTimeGraph,
    TSGraph::TimeSpaceGraph,
    bundle::Bundle,
    src::Int,
    dst::Int,
    CAPACITIES::Vector{Int};
    sorted::Bool=false,
    use_bins::Bool=true,
    opening_factor::Float64=1.0,
    current_cost::Bool=false,
    verbose::Bool=false,
)
    # Otherwise, computing the new cost
    arcBundleCost = EPS
    verbose && println(
        "\nComputing tentative cost of bundle $bundle for arc $src -> $dst (sorted = $sorted, use_bins = $use_bins, current_cost = $current_cost, opening_factor = $opening_factor)",
    )
    verbose && println("arcBundleCost = $arcBundleCost")
    for order in bundle.orders
        costBefore = arcBundleCost
        # Getting time space projection
        tSrc, tDst = time_space_projector(TTGraph, TSGraph, src, dst, order)
        verbose && println("\nComputing for order $order")
        verbose && println("Timed projection $tSrc -> $tDst")
        # Node volume cost 
        arcBundleCost += volume_stock_cost(TTGraph, src, dst, order; verbose=verbose)
        verbose && println(
            "Volume cost added for order : $(volume_stock_cost(TTGraph, src, dst, order))",
        )
        # Arc transport cost 
        units = transport_units(
            sol,
            TSGraph,
            tSrc,
            tDst,
            order,
            CAPACITIES;
            sorted=sorted,
            use_bins=use_bins,
            verbose=verbose,
        )
        verbose && println("Transport units added for order : $units")
        verbose && println("Transport cost : $(transport_cost(TSGraph, tSrc, tDst))")
        verbose && println(
            "Transport cost added : $(units * transport_cost(TSGraph, tSrc, tDst) * opening_factor)",
        )
        arcBundleCost +=
            units *
            transport_cost(TSGraph, tSrc, tDst; current_cost=current_cost) *
            opening_factor
        verbose && println("Cost added for order : $(arcBundleCost - costBefore)")
        verbose && println("\narcBundleCost = $arcBundleCost")
    end
    return arcBundleCost
end

function find_other_src_node(travelTimeGraph::TravelTimeGraph, src::Int)::Int
    otherSrcIdx = findfirst(
        dst -> travelTimeGraph.networkArcs[src, dst].type == :shortcut,
        outneighbors(travelTimeGraph.graph, src),
    )
    otherSrcIdx === nothing && return -1
    return outneighbors(travelTimeGraph.graph, src)[otherSrcIdx::Int]
end

# Creating start node vector
function get_all_start_nodes(travelTimeGraph::TravelTimeGraph, bundle::Bundle)
    src = travelTimeGraph.bundleSrc[bundle.idx]
    startNodes = Int[src]
    # Iterating through outneighbors of the start node
    otherSrc = find_other_src_node(travelTimeGraph, src)
    # Iterating through outneighbors of the other start node 
    while otherSrc != -1
        push!(startNodes, otherSrc)
        src = otherSrc
        otherSrc = find_other_src_node(travelTimeGraph, src)
    end
    return startNodes
end

# Updating cost matrix on the travel time graph for a specific bundle 
function update_cost_matrix!(
    solution::Solution,
    TTGraph::TravelTimeGraph,
    TSGraph::TimeSpaceGraph,
    bundle::Bundle,
    CAPACITIES::Vector{Int};
    sorted::Bool=false,
    use_bins::Bool=true,
    opening_factor::Float64=1.0,
    current_cost::Bool=false,
    findSources::Bool=true,
)
    # Iterating through outneighbors of the start nodes and common nodes
    startNodes = findSources ? get_all_start_nodes(TTGraph, bundle) : Int[]
    for src in vcat(startNodes, TTGraph.commonNodes)
        for dst in outneighbors(TTGraph.graph, src)
            # If the arc doesn't need an update, skipping
            is_update_candidate(TTGraph, src, dst, bundle) || continue
            # Otherwise, computing the new cost
            TTGraph.costMatrix[src, dst] = arc_update_cost(
                solution,
                TTGraph,
                TSGraph,
                bundle,
                src,
                dst,
                CAPACITIES;
                sorted=sorted,
                use_bins=use_bins,
                opening_factor=opening_factor,
                current_cost=current_cost,
            )
        end
    end
end

# Check whether the path of the bundle needs to be recomputed
function is_path_admissible(travelTimeGraph::TravelTimeGraph, path::Vector{Int})
    # Checking elementarity on network
    return is_path_elementary(travelTimeGraph, path)
    # Too long ? Too many node ? To be difined
end

function path_cost(path::Vector{Int}, costMatrix::SparseMatrixCSC{Float64,Int})
    cost = 0.0
    for (i, j) in partition(path, 2, 1)
        cost += costMatrix[i, j]
    end
    return cost
end

######################################################
###           Parallel Implementation              ###
######################################################

# Creating channel of CAPACITIES to limit memory footprint
function create_channel(::Type{Vector{Int}}; n::Int=Threads.nthreads())
    chnl = Channel{Vector{Int}}(n)
    foreach(1:n) do _
        put!(chnl, Vector{Int}(undef, 0))
    end
    return chnl
end

function create_filled_channel(; n::Int=Threads.nthreads(), bufferSize::Int=200)
    chnl = Channel{Vector{Int}}(n)
    foreach(1:n) do _
        put!(chnl, [0 for _ in 1:bufferSize])
    end
    return chnl
end

function parallel_update_cost_matrix!(
    solution::Solution,
    TTGraph::TravelTimeGraph,
    TSGraph::TimeSpaceGraph,
    bundle::Bundle,
    CHANNEL::Channel{Vector{Int}},
    sorted::Bool=true,
    use_bins::Bool=true,
    opening_factor::Float64=1.0,
)
    # Iterating in parallel (thanks to @tasks) through the bundle arcs
    @tasks for (src, dst) in TTGraph.bundleArcs[bundle.idx]
        CAPACITIES = take!(CHANNEL)
        TTGraph.costMatrix[src, dst] = if is_update_candidate(TTGraph, src, dst, bundle)
            arc_update_cost(
                solution,
                TTGraph,
                TSGraph,
                bundle,
                src,
                dst,
                CAPACITIES;
                sorted=sorted,
                use_bins=use_bins,
                opening_factor=opening_factor,
            )
        else
            EPS
        end
        put!(CHANNEL, CAPACITIES)
    end
end
