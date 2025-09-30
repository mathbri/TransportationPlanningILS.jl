# Updating functions for the bins

function compute_new_cost(
    arcData::NetworkArc, dstData::NetworkNode, newBins::Int, commodities::Vector{Commodity}
)
    volume = sum(com.size for com in commodities)
    leadTimeCost = sum(com.stockCost for com in commodities)
    # Node cost 
    cost =
        dstData.volumeCost * volume / VOLUME_FACTOR +
        arcData.carbonCost * volume / arcData.capacity
    # Transport cost 
    addedBins = arcData.isLinear ? (volume / arcData.capacity) : newBins
    cost += addedBins * arcData.unitCost
    # Commodity cost
    return cost += arcData.distance * leadTimeCost
end

# Add order content to solution truck loads with packing function
function add_order!(
    solution::Solution,
    TSGraph::TimeSpaceGraph,
    timedPath::Vector{Int},
    order::Order;
    sorted::Bool=false,
    verbose::Bool=false,
)
    if verbose
        return add_order_verbosed!(solution, TSGraph, timedPath, order; sorted=sorted)
    end
    costAdded = 0.0
    for (timedSrc, timedDst) in partition(timedPath, 2, 1)
        bins = solution.bins[timedSrc, timedDst]
        dstData = TSGraph.networkNodes[timedDst]
        arcData = TSGraph.networkArcs[timedSrc, timedDst]
        # Updating bins
        addedBins = first_fit_decreasing!(bins, arcData, order; sorted=sorted)
        # Updating cost
        costAdded += compute_new_cost(arcData, dstData, addedBins, order.content)
    end
    return costAdded
end

function add_order_verbosed!(
    solution::Solution,
    TSGraph::TimeSpaceGraph,
    timedPath::Vector{Int},
    order::Order;
    sorted::Bool=false,
)
    costAdded = 0.0
    println("\nAdding order $(order) on timed path $timedPath")
    for (timedSrc, timedDst) in partition(timedPath, 2, 1)
        println("\nAdding on arc $timedSrc -> $timedDst")
        bins = solution.bins[timedSrc, timedDst]
        dstData = TSGraph.networkNodes[timedDst]
        arcData = TSGraph.networkArcs[timedSrc, timedDst]
        println("Dst node $timedDst : $(dstData)")
        println("Arc $timedSrc -> $timedDst : $(arcData)")
        previousLoads = [bin.load for bin in bins]
        println("Bins : $(length(bins)) -> $([bin.load for bin in bins]) m3")
        # Updating bins
        addedBins = first_fit_decreasing!(bins, arcData, order; sorted=sorted)
        println("Order added : $addedBins new bins")
        if length(previousLoads) < length(bins)
            append!(previousLoads, zeros(length(bins) - length(previousLoads)))
        end
        println(
            "Bins : $(length(bins)) -> $([bin.load for bin in bins] .- previousLoads) m3"
        )
        println("Bins : $(length(bins)) -> $([bin.load for bin in bins]) m3")
        # Updating cost
        costAddedForOrder = compute_new_cost(arcData, dstData, addedBins, order.content)
        println("Cost added for order : $costAddedForOrder")
        costAdded += compute_new_cost(arcData, dstData, addedBins, order.content)
    end
    return costAdded
end

# Remove order content from solution truck loads, does not refill bins
function remove_order!(
    solution::Solution,
    TSGraph::TimeSpaceGraph,
    timedPath::Vector{Int},
    order::Order;
    verbose::Bool=false,
)
    if verbose
        return remove_order_verbosed!(solution, TSGraph, timedPath, order)
    end
    costAdded, orderUniqueCom = 0.0, unique(order.content)
    # For all arcs in the path, updating the right bins
    for (timedSrc, timedDst) in partition(timedPath, 2, 1)
        for bin in solution.bins[timedSrc, timedDst]
            remove!(bin, orderUniqueCom)
        end
        dstData = TSGraph.networkNodes[timedDst]
        arcData = TSGraph.networkArcs[timedSrc, timedDst]
        costAdded -= compute_new_cost(arcData, dstData, 0, order.content)
    end
    return costAdded
end

function remove_order_verbosed!(
    solution::Solution, TSGraph::TimeSpaceGraph, timedPath::Vector{Int}, order::Order
)
    costAdded, orderUniqueCom = 0.0, unique(order.content)
    println("\nRemoving order $(order) on timed path $timedPath")
    # For all arcs in the path, updating the right bins
    for (timedSrc, timedDst) in partition(timedPath, 2, 1)
        println("\nRemoving on arc $timedSrc -> $timedDst")
        arcBins = solution.bins[timedSrc, timedDst]
        previousLoads = [bin.load for bin in arcBins]
        println("Bins : $(length(arcBins)) -> $([bin.load for bin in arcBins]) m3")
        for bin in solution.bins[timedSrc, timedDst]
            remove!(bin, orderUniqueCom)
        end
        println("Order removed")
        println(
            "Bins : $(length(arcBins)) -> $([bin.load for bin in arcBins] .- previousLoads) m3",
        )
        println("Bins : $(length(arcBins)) -> $([bin.load for bin in arcBins]) m3")
        dstData = TSGraph.networkNodes[timedDst]
        println("Dst Node $timedDst : $(dstData)")
        arcData = TSGraph.networkArcs[timedSrc, timedDst]
        println("Arc $timedSrc -> $timedDst : $(arcData)")
        costRemoved = compute_new_cost(arcData, dstData, 0, order.content)
        println("Cost removed : $costRemoved")
        costAdded -= compute_new_cost(arcData, dstData, 0, order.content)
        println("Total cost added : $costAdded")
    end
    return costAdded
end

function update_bins!(
    solution::Solution,
    TSGraph::TimeSpaceGraph,
    TTGraph::TravelTimeGraph,
    bundle::Bundle,
    path::Vector{Int};
    sorted::Bool=false,
    remove::Bool=false,
    verbose::Bool=false,
)
    costAdded = 0.0
    verbose && println(
        "\nUpdating bins for bundle $bundle (bunH = $(bundle.hash) on path travel time path $path",
    )
    for order in bundle.orders
        verbose && println("Updating order $order")
        # Projecting path
        timedPath = time_space_projector(TTGraph, TSGraph, path, order)
        verbose && println("Timed path : $timedPath")
        if -1 in timedPath
            bundleSrcDst = (TTGraph.bundleSrc[bundle.idx], TTGraph.bundleDst[bundle.idx])
            pathStr = join(path, ", ")
            pathInfo = join(string.(TTGraph.networkNodes[path]), ", ")
            pathSteps = join(string.(TTGraph.stepToDel[path]), ", ")
            timedPathStr = join(timedPath, ", ")
            @error "At least one node was not projected in bin updating" :bundle = bundle :bundleSrcDst =
                bundleSrcDst :order = order :path = pathStr :pathInfo = pathInfo :pathSteps =
                pathSteps :timedPath = timedPathStr
        end
        # Add or Remove order
        if remove
            verbose && println("\nRemoving order $(order)")
            costAdded += remove_order!(solution, TSGraph, timedPath, order; verbose=verbose)
            verbose && println("Cost removed : $costAdded")
        else
            verbose && println("\nAdding order $(order)")
            costAdded += add_order!(
                solution, TSGraph, timedPath, order; sorted=sorted, verbose=verbose
            )
            verbose && println("Cost added for order : $costAdded")
        end
    end
    return costAdded
end
